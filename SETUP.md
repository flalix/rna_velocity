# Environment setup — Neural Crest
#
# Every command below was executed against pyproject_corrected.toml and the
# resulting venv was import-checked (Python 3.12.3; yaml 6.0.3, numpy 2.5.3,
# scanpy 1.12.4, squidpy 1.8.3).

# ---------------------------------------------------------------- quick start

# 1. Pin the interpreter. Writes `.python-version`; commit it so everyone
#    resolves against the same minor version.
uv python pin 3.12

# 2. Create the venv and install. `uv sync` creates `.venv` itself using the
#    pin, so a separate `uv venv` step is unnecessary.
uv sync --group dev

# 3. acitvate

source .venv/bin/activate

# 4. Verify.
uv run python -c "import yaml, numpy, scanpy, squidpy; print('ok')"


# ------------------------------------------------------- what changed and why

# The snippet this replaces was:
#
#   uv python pin 3.12
#   uv venv                       # uses the pinned 3.12
#   uv sync --extra dev           # add --extra llm if you want the LangChain stack
#
# Three problems:
#
# `uv venv` is redundant. `uv sync` creates the project environment itself and
# honours `.python-version`. Running `uv venv` first is harmless but does no
# work that the next command would not do. It also misleads: if you run
# `uv venv` with no pin present, you get whatever interpreter uv discovers, and
# `uv sync` may then replace that environment.
#
# `--extra dev` fails. Fails with:
#     error: Extra `dev` is not defined in the `optional-dependencies` table
#   `--extra` reads PEP 621 `[project.optional-dependencies]`. `dev` is a PEP 735
#   entry in `[dependency-groups]`, which is `--group`. The two flags are not
#   interchangeable. This project defines no extras at all, so `--extra` is never
#   the right flag here.
#
# `--extra llm` does not exist. There is no `llm` extra and no LangChain
#   dependency anywhere in this project — not in the original file, not in the
#   corrected one. Nothing in the HDCA heart pipeline uses an LLM stack; the
#   dependencies are scanpy/squidpy/spatialdata for spatial transcriptomics.
#   The comment appears to have been carried over from an unrelated template.


# ---------------------------------------------------------- available groups

# analysis   core stack — scanpy, squidpy, spatialdata, numpy, pyyaml (default)
# viz        napari-spatialdata — interactive viewer, pulls the Qt GUI stack
# image      centrosome, cp-measure, dask-image — image featurization
# test       analysis + pytest, coverage, pytest-xdist
# docs       analysis + sphinx, myst-nb
# dev        test + viz + pre-commit, ruff

uv sync                        # analysis only (default-groups)
uv sync --group dev            # full development environment
uv sync --group test           # CI
uv sync --group analysis --group image   # headless cluster job with image features

# Nothing here needs `--extra`; this project defines no extras.


# ------------------------------------------------------------------ with pixi

# The same environments are also declared under [tool.pixi], which manages the
# conda-side dependencies the original repo relied on:

pixi install -e dev-py312
pixi run -e dev-py312 test

# Use one tool or the other per checkout, not both — they maintain separate
# environments and separate lockfiles.
