"""Run Spark on the Terra task VM and collect checked WDL outputs."""
import argparse
import json
import os
import subprocess
import tempfile
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from .contracts import (GLOW_COMMIT, LOG, SPARK_ARROW_PROPERTIES, InputArtifact,
                        add_scientific_arguments, options_from_args, require_local_file,
                        scientific_arguments, setup_logging)
from .storage import collect_outputs, input_identity

@dataclass(frozen=True)
class LocalOptions:
    cpu: int = 64
    memory_gb: int = 256
    driver_memory_gb: int = 64
    parallelism: int = 100
    shuffle_partitions: int = 100

    def __post_init__(self):
        if min(self.cpu,self.memory_gb,self.driver_memory_gb,self.parallelism,self.shuffle_partitions)<1:
            raise ValueError('Spark CPU, memory, and partition settings must be positive')
        if self.driver_memory_gb>=self.memory_gb:
            raise ValueError('Spark driver memory must be less than VM memory; reserve space for Python and the operating system')

def spark_command(options,scratch):
    properties={**SPARK_ARROW_PROPERTIES,
                'spark.default.parallelism':str(options.parallelism),
                'spark.sql.shuffle.partitions':str(options.shuffle_partitions),
                'spark.local.dir':str(scratch),
                'spark.dynamicAllocation.enabled':'false',
                'spark.ui.enabled':'false',
                'spark.executorEnv.OPENBLAS_NUM_THREADS':'1',
                'spark.executorEnv.OMP_NUM_THREADS':'1'}
    return ['spark-submit','--master',f'local[{options.cpu}]','--driver-memory',f'{options.driver_memory_gb}g',
            *[item for key,value in properties.items() for item in ('--conf',key+'='+value)],
            '--jars','/opt/glow/jars/glow-assembly.jar']

def read_file_list(path):
    paths=[require_local_file(Path(line)).resolve() for line in require_local_file(path).read_text().splitlines()]
    if not paths: raise ValueError('Empty association input list')
    return paths

def run_local(options,science,inputs,out_dir):
    out_dir=Path(out_dir).resolve(); out_dir.mkdir(parents=True,exist_ok=True)
    inputs=[InputArtifact(item.name,require_local_file(item.path).resolve()) for item in inputs]
    identity=input_identity(inputs,{'science':asdict(science),'runtime':asdict(options),'glow_commit':GLOW_COMMIT,
                                    'project_revision':os.environ.get('PROJECT_REVISION','unknown')})
    by_name={item.name:item.path for item in inputs}
    report={'state':'RUNNING','execution_mode':'single_vm','spark_master':f'local[{options.cpu}]',
            'resources':asdict(options),'identity_sha256':identity['sha256']}
    # Keep work on the task disk. Cromwell collects only the declared result files.
    with tempfile.TemporaryDirectory(prefix='glow-local-',dir=out_dir.parent) as directory:
        work=Path(directory); scratch=work/'scratch'; scratch.mkdir()
        # Stable names prevent Hadoop glob characters in localized filenames.
        for name,path in list(by_name.items()):
            alias=work/'inputs'/name; alias.parent.mkdir(parents=True,exist_ok=True)
            alias.symlink_to(path); by_name[name]=alias
        run_uri=(work/'run').as_uri()
        command=spark_command(options,scratch)+[str(Path(__file__).with_name('spark_main.py'))]
        for name,flag in (('model/bgen','model-bgen'),('model/sample','model-sample'),('model/variants','model-variants'),
                          ('sample_map','sample-map'),('phenotypes','phenotypes'),('covariates','covariates'),('phenotype_batches','phenotype-batches')):
            command+=['--'+flag,by_name[name].as_uri()]
        indices=sorted({int(name.split('/')[1]) for name in by_name if name.startswith('association/')})
        if not indices or indices!=list(range(len(indices))): raise ValueError('Association sources must be nonempty and contiguous')
        for index in indices:
            for kind in ('bgen','sample','variants'): command+=['--association-'+kind,by_name[f'association/{index}/{kind}'].as_uri()]
        command+=['--output-prefix',run_uri,'--identity-sha256',identity['sha256']]+scientific_arguments(science)
        env=dict(os.environ,SPARK_LOCAL_DIRS=str(scratch),SPARK_LOCAL_IP='127.0.0.1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        try:
            LOG.info('Start Glow on one VM: %s, driver heap %d GB',report['spark_master'],options.driver_memory_gb)
            # Inherit stdout/stderr so the WDL task log includes Spark failures.
            subprocess.run(command,env=env,check=True)
            paths=collect_outputs(None,run_uri,identity['sha256'],out_dir)
            # Local Spark URIs expire when this temporary directory is removed.
            # Keep manifest paths as relative metadata for the collected WDL Files.
            manifest=json.loads((out_dir/'output_manifest.json').read_text())
            for row in manifest['outputs']: row.pop('uri',None)
            (out_dir/'output_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
            report['state']='SUCCEEDED'
            LOG.info('Glow local analysis completed; %d files collected',len(paths))
            return paths
        except Exception as error:
            report.update(state='FAILED',error=str(error)); raise
        finally:
            (out_dir/'execution.json').write_text(json.dumps(report,indent=2)+'\n')

def main(argv=None):
    setup_logging(); parser=argparse.ArgumentParser(description=__doc__)
    for field in fields(LocalOptions): parser.add_argument('--'+field.name.replace('_','-'),type=int,default=field.default)
    for name in ('model-bgen','model-sample','model-variants','sample-map','phenotypes','covariates','phenotype-batches','out-dir'):
        parser.add_argument('--'+name,type=Path,required=True)
    for kind in ('bgen','sample','variants'): parser.add_argument('--association-'+kind+'-list',type=Path,required=True)
    add_scientific_arguments(parser); args=parser.parse_args(argv)
    options=LocalOptions(**{field.name:getattr(args,field.name) for field in fields(LocalOptions)})
    inputs=[InputArtifact(name,require_local_file(getattr(args,flag))) for name,flag in (
        ('model/bgen','model_bgen'),('model/sample','model_sample'),('model/variants','model_variants'),('sample_map','sample_map'),
        ('phenotypes','phenotypes'),('covariates','covariates'),('phenotype_batches','phenotype_batches'))]
    groups=[read_file_list(getattr(args,'association_'+kind+'_list')) for kind in ('bgen','sample','variants')]
    if len({len(group) for group in groups})!=1: raise ValueError('Association input lists differ in length')
    for index,triplet in enumerate(zip(*groups)):
        inputs.extend(InputArtifact(f'association/{index}/{kind}',path) for kind,path in zip(('bgen','sample','variants'),triplet))
    inputs.append(InputArtifact('application/spark_main.py',Path(__file__).with_name('spark_main.py')))
    run_local(options,options_from_args(args),inputs,args.out_dir)
    return 0
if __name__=='__main__': raise SystemExit(main())
