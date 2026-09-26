# 1. Pin the interpreter for this project (writes .python-version)
uv python install 3.12        # downloads a managed CPython if you don't have one
uv python pin 3.12

# 2. Create the venv and resolve/install everything
uv venv                       # uses the pinned 3.12
conda actiate renv
uv sync --group dev

# 3. Sanity check
uv run python -c "import yaml, numpy, scanpy, squidpy; print('ok')"