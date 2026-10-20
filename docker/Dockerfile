FROM --platform=linux/amd64 condaforge/miniforge3:latest

ENV DEBIAN_FRONTEND=noninteractive \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/envs/mfix-26.1.2/bin:/opt/conda/bin:$PATH \
    MFIX_ENV=mfix-26.1.2 \
    OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OMPI_MCA_hwloc_base_binding_policy=core

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        ca-certificates curl less procps unzip \
 && rm -rf /var/lib/apt/lists/* \
 && curl -fsSL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o /tmp/awscliv2.zip \
 && unzip -q /tmp/awscliv2.zip -d /tmp \
 && /tmp/aws/install \
 && rm -rf /tmp/aws /tmp/awscliv2.zip

RUN mamba create -y -n mfix-26.1.2 \
        mfix=26.1.2 mfix-solver \
        -c conda-forge \
        -c https://mfix.netl.doe.gov/s3/4a4141af/f7bb3848a7d3c9b90a9d380e687cc523//conda/dist \
 && mamba clean -afy

WORKDIR /work

COPY entrypoint.sh /opt/entrypoint.sh
RUN chmod +x /opt/entrypoint.sh

ENTRYPOINT ["/opt/entrypoint.sh"]
