FROM mambaorg/micromamba:2.3.3

LABEL org.opencontainers.image.source="https://github.com/AoU-Multiomics-Analysis/regenie"

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/environment.yml /tmp/environment.yml
RUN micromamba install --yes --strict-channel-priority --name base --file /tmp/environment.yml \
    && micromamba clean --all --yes

USER root
COPY scripts/*.py /opt/regenie/scripts/
RUN chmod -R a+rX /opt/regenie

# Terra can start a shell without the micromamba entrypoint.
ENV PATH="/opt/conda/bin:${PATH}"
ENV PYTHONUNBUFFERED=1
WORKDIR /work
CMD ["python3", "/opt/regenie/scripts/run_regenie.py", "--help"]
