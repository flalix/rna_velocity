#!/usr/bin/env bash
# GSE148190 / patient K409: counting -> QC -> PCA -> clusters -> velocity.
# Run from the directory where the outputs should land (count_*, *.h5ad, *.csv).
#
#   KBREF=/path/to/kb_nac_human FQDIR=/path/to/fastq bash run_all.sh
#
# Environments: kb count needs kb-python; everything downstream needs numpy<2
# (scvelo 0.3.4 is not numpy-2 compatible). Set VELO_PY to that interpreter.
set -euo pipefail

SRC=/home/flavio/uv/rna_velocity/src/scripts
KBREF=${KBREF:-kbref}
FQDIR=${FQDIR:-fastq}
VELO_PY=${VELO_PY:-python}
COUNT_ROOT=${COUNT_ROOT:-.}   # kb count writes count_<library>/ here; load_qc.py reads from here


echo "SRC: $SRC"
echo "KBREF: $KBREF"
echo "FQDIR: $FQDIR"
echo "VELO_PY: $VELO_PY"
echo "COUNT_ROOT: $COUNT_ROOT"
