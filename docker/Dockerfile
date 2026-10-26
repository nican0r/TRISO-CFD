FROM --platform=linux/amd64 condaforge/miniforge3:latest

ENV DEBIAN_FRONTEND=noninteractive \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/envs/mfix-23.2/bin:/opt/conda/bin:$PATH \
    MFIX_ENV=mfix-23.2 \
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

# MFiX 23.2 as NETL's Linux install command gives it, plus a toolchain for
# the DMP build below.  The stock 23.2 mfix-solver is serial only (bin/
# mfixsolver; no MPI, no OpenMP; openmpi first appears as a dependency in
# 23.3.2).  gfortran 13 matches the libgfortran5 >= 13.1 the stock 23.2
# solver was built against; binutils supplies the unprefixed ar/ranlib CMake
# archives with (gcc brings only x86_64-conda-linux-gnu-ar and the base image
# has none).  The channel URL is the per-download link from
# NETL's registration e-mail and expires (the 26.1.2 one now returns 401).
# conda, not mamba: mamba 2.x's sharded conda-forge index fetches only the
# shards of the named packages, so the MFiX channel's conda-forge
# dependencies (flask, matplotlib, ...) "do not exist" (2026-10-24 build).
RUN conda create -y -n mfix-23.2 \
        mfix=23.2 mfix-gui=23.2 mfix-src=23.2 mfix-doc=23.2 mfix-solver=23.2 \
        gcc=13 gfortran=13 binutils openmpi=4.1 cmake make \
        -c conda-forge \
        -c https://mfix.netl.doe.gov/s3/4a4141af/e98cde21edc4d275274d8084577accfa//conda/dist \
 && conda clean -afy

# Batch DMP solver from the 23.2 source, via MFiX's own build tool
# (mfixgui/build_mfixsolver.py: --batch = plain Fortran executable, no
# pymfix; --dmp = ENABLE_MPI), run as a module (same code as the
# bin/build_mfixsolver launcher); CONDA_PREFIX lets it find share/mfix/src.
# Installed as mfixsolver_dmp, the name the entrypoint uses.
ARG BUILD_JOBS=4
RUN mkdir /tmp/mfix-dmp && cd /tmp/mfix-dmp \
 && CONDA_PREFIX=/opt/conda/envs/mfix-23.2 \
        python -m mfixgui.build_mfixsolver --batch --dmp -j "$BUILD_JOBS" \
 && install -m 755 mfixsolver /opt/conda/envs/mfix-23.2/bin/mfixsolver_dmp \
 && cd / && rm -rf /tmp/mfix-dmp

# Open MPI 4.1's mpirun selects its rsh launcher even on one node and aborts
# when neither ssh nor rsh exists (none here); "isolated" launches every
# rank locally, which is all a single-instance Batch job needs.  (26.1.2's
# image had Open MPI 5, which did not need this.)
ENV OMPI_MCA_plm=isolated

WORKDIR /work

COPY entrypoint.sh /opt/entrypoint.sh
RUN chmod +x /opt/entrypoint.sh

ENTRYPOINT ["/opt/entrypoint.sh"]
