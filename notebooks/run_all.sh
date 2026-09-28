#!/usr/bin/env bash
# GSE148190 / patient K409: counting -> QC -> PCA -> clusters -> velocity.
# Run from the directory where the outputs should land (count_*, *.h5ad, *.csv).
#
#   KBREF=/path/to/kb_nac_human FQDIR=/path/to/fastq bash run_all.sh
#
# Environments: kb count needs kb-python; everything downstream needs numpy<2
# (scvelo 0.3.4 is not numpy-2 compatible). Set VELO_PY to that interpreter.
set -euo pipefail
SRC=$(cd "$(dirname "$0")" && pwd)
KBREF=${KBREF:-kbref}
FQDIR=${FQDIR:-fastq}
VELO_PY=${VELO_PY:-python}
export KBREF FQDIR

echo "== 1/6 kb count (three libraries, ~11 min each at 8 threads) =="
bash "$SRC/count_libs.sh" tumor_primary:SRR11492025 lymph_node:SRR11492017 blood_PBMC:SRR11492013

echo "== 2/6 author cell barcodes from GEO =="
[ -s geo_barcodes/author_cells.json ] || "$VELO_PY" "$SRC/fetch_author_barcodes.py"

echo "== 3/6 load counts, call cells =="
"$VELO_PY" "$SRC/load_qc.py"

echo "== 4/6 QC metrics, doublets, filtering =="
"$VELO_PY" "$SRC/step_qc.py"

echo "== 5/6 normalise, HVGs, PCA diagnostics, clusters =="
"$VELO_PY" "$SRC/step_pca.py"
"$VELO_PY" "$SRC/step_cluster.py"

echo "== 6/6 scVelo (stochastic + dynamical; ~20 min single-process) =="
"$VELO_PY" "$SRC/step_velocity.py"

echo "done. key outputs:"
ls -1 k409_*.h5ad *.csv *.json 2>/dev/null | sed 's/^/  /'
