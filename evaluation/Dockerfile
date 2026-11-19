# Reproducible Evaluation Image for OpenBench (Polymer + Protein-Ligand)
#
# Uses a single unified OpenStructure 2.12.0 runtime (`/usr/local/bin/ost`) for both
# Polymer tasks and Runs N' Poses Protein-Ligand evaluation, refreshes the wwPDB CCD
# library (`compounds.chemlib`), and installs RDKit, PoseBusters, and Gemmi into an
# isolated `/opt/eval-venv` so OpenStructure's system Python/C++ ABI remains untouched.
#
# Build (from OpenBench/ root):
#   docker build -t openbench-eval -f evaluation/Dockerfile evaluation/

FROM registry.scicore.unibas.ch/schwede/openstructure:2.12.0

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    UV_NO_CACHE=1 \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1

ENTRYPOINT []
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates curl procps \
    && mkdir -p /tmp/ccd /usr/local/share/openstructure \
    && curl -fsSL https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz \
        -o /tmp/ccd/components.cif.gz \
    && /usr/local/bin/chemdict_tool create \
        /tmp/ccd/components.cif.gz \
        /usr/local/share/openstructure/compounds.chemlib \
        pdb -i \
    && curl -LsSf https://astral.sh/uv/install.sh -o /tmp/uv-install.sh \
    && sh /tmp/uv-install.sh \
    && install -m 0755 /root/.local/bin/uv /usr/local/bin/uv \
    && uv venv /opt/eval-venv --python /usr/bin/python3 \
    && uv pip install --python /opt/eval-venv/bin/python --no-cache \
        rdkit \
        gemmi \
        posebusters \
        pandas \
        numpy \
        scipy \
        biopython \
        tqdm \
    && apt-get purge -y --auto-remove curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/ccd /tmp/uv-install.sh /root/.local /root/.cache

ENV OST_COMPOUNDS_CHEMLIB=/usr/local/share/openstructure/compounds.chemlib
ENV PATH=/opt/eval-venv/bin:/usr/local/bin:$PATH

COPY . /opt/openbench/evaluation

# Smoke test unified OpenStructure (`ost`) and Python evaluation dependencies
RUN ost --version \
    && ost compare-structures --help >/dev/null \
    && ost compare-ligand-structures --help >/dev/null \
    && python3 -c "import rdkit, gemmi, pandas, numpy, scipy, posebusters; print('OpenBench unified eval environment ready')"

WORKDIR /workspace
CMD ["/bin/bash"]
