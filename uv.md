`default-groups = ["analysis"]` in your `[tool.uv]` is the cause — `analysis` is always on, so `--group velocity` enables *both*, which is exactly what `conflicts` forbids. You need to disable the default, not just add the group. I reproduced your error verbatim on uv 0.12.19 and verified both fixes:

```bash
# either
uv sync --only-group velocity --python 3.12

# or
uv sync --no-default-groups --group velocity --python 3.12
```

Both land: Python 3.12.14, numpy 1.26.4, scanpy 1.11.5, scvelo 0.3.4, kb_python 0.30.2.

**One thing to know before you use it:** both forms write into the project's `.venv`, so syncing the velocity group *replaces* your analysis environment, and a later plain `uv sync` swaps it back. If you want both on disk simultaneously, give the velocity group its own path:

```bash
UV_PROJECT_ENVIRONMENT=.venv-velo uv sync --only-group velocity --python 3.12
uv sync --python 3.12          # .venv stays on the analysis stack
```

Verified side by side just now:

```
.venv-velo   numpy 1.26.4   scvelo 0.3.4
.venv        numpy 2.5.3    scanpy 1.11.5
```

Then run the pipeline against the pinned interpreter:

```bash
VELO_PY=$PWD/.venv-velo/bin/python KBREF=$WSref/kb_nac_human FQDIR=$WSuv/fastq \
  bash $WSsrc/run_all.sh
```

Worth adding `.venv-velo` to `.gitignore` alongside `.venv`, and a one-line comment on the `velocity` group recording *why* it's forked (scvelo 0.3.4 is not numpy-2 safe) so the fork can be deleted when that's fixed upstream rather than becoming permanent furniture.