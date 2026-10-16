FROM mambaorg/micromamba:2.3.3 AS glow-builder
USER root
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
ENV PATH="/opt/conda/bin:${PATH}"
COPY envs/glow/build-environment.yml /tmp/build-environment.yml
RUN micromamba install --yes --strict-channel-priority --name base --file /tmp/build-environment.yml \
    && micromamba clean --all --yes
RUN git init /opt/glow-source \
    && git -C /opt/glow-source remote add origin https://github.com/projectglow/glow.git \
    && git -C /opt/glow-source fetch --depth 1 origin 72637ad17ff267ca6bbc63cfc77c78db619f38b5 \
    && git -C /opt/glow-source checkout --detach FETCH_HEAD
RUN mkdir -p /opt/sbt \
    && curl --fail --location --retry 3 https://repo.maven.apache.org/maven2/org/scala-sbt/sbt-launch/1.9.9/sbt-launch-1.9.9.jar -o /opt/sbt/sbt-launch.jar
COPY scripts/glow/ /opt/regenie-spark/scripts/glow/
ENV PYTHONPATH="/opt/regenie-spark"
RUN python -m scripts.glow.build_artifacts --source-dir /opt/glow-source --out-dir /opt/artifacts --sbt-launcher /opt/sbt/sbt-launch.jar

FROM mambaorg/micromamba:2.3.3 AS production
USER root
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
ARG PROJECT_REVISION=unknown
LABEL org.opencontainers.image.source="https://github.com/AoU-Multiomics-Analysis/regenie"
LABEL org.opencontainers.image.revision="${PROJECT_REVISION}"
ENV PROJECT_REVISION="${PROJECT_REVISION}" \
    PATH="/opt/conda/bin:${PATH}" \
    PYTHONPATH="/opt/regenie-spark" \
    PYSPARK_PYTHON="/opt/conda/bin/python" \
    SPARK_EXTRA_CLASSPATH="/opt/glow/jars/*" \
    PYTHONUNBUFFERED=1 \
    OPENBLAS_NUM_THREADS=1 \
    OMP_NUM_THREADS=1 \
    HOME="/home/spark"
COPY envs/glow/environment.yml /tmp/environment.yml
RUN micromamba install --yes --strict-channel-priority --name base --file /tmp/environment.yml \
    && micromamba clean --all --yes
COPY --from=glow-builder /opt/artifacts/ /tmp/glow-artifacts/
RUN python -m pip install --no-deps /tmp/glow-artifacts/*.whl \
    && python -m pip check \
    && mkdir -p /opt/glow/jars /work \
    && cp /tmp/glow-artifacts/glow-assembly.jar /opt/glow/jars/ \
    && cp /tmp/glow-artifacts/build.json /opt/glow/ \
    && rm -rf /tmp/glow-artifacts \
    && groupadd --gid 1099 spark \
    && useradd --uid 1099 --gid 1099 --create-home --shell /bin/bash spark \
    && chown -R 1099:1099 /work /home/spark
COPY scripts/glow/ /opt/regenie-spark/scripts/glow/
RUN chmod -R a+rX /opt/regenie-spark /opt/glow \
    && python -c "import importlib.util, shutil; assert importlib.util.find_spec('pyspark') is None; assert shutil.which('java') is None" \
    && micromamba list --explicit --name base > /opt/glow/environment.explicit.txt \
    && python -m pip freeze > /opt/glow/python-packages.txt
WORKDIR /work
USER 1099:1099
CMD ["python", "-m", "scripts.glow.controller", "--help"]

FROM production AS ci-test
USER root
RUN micromamba install --yes --strict-channel-priority --name base --channel conda-forge pyspark=3.5.1 openjdk=17.0.15 \
    && micromamba clean --all --yes
ENV JAVA_HOME="/opt/conda/lib/jvm" \
    SPARK_HOME="/opt/conda/lib/python3.11/site-packages/pyspark" \
    SPARK_LOCAL_IP="127.0.0.1"
USER 1099:1099
CMD ["python", "-m", "scripts.glow.smoke_test", "--work-dir", "/work/smoke"]
