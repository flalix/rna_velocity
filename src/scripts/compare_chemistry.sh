#!/usr/bin/env bash
# Compare kb count settings (10x chemistry x strand) on one pair of FASTQs.
#
# Usage:
#   bash compare_chemistry.sh <R1.fastq.gz> <R2.fastq.gz> <index_dir> [outdir] [threads]
#
# <index_dir> must contain index.idx, t2g.txt, cdna_t2c.txt, nascent_t2c.txt
# (build once with:  kb ref -d human --workflow nac -i index.idx -g t2g.txt \
#                      -c1 cdna_t2c.txt -c2 nascent_t2c.txt )
set -euo pipefail

R1=${1:?need R1 fastq.gz}
R2=${2:?need R2 fastq.gz}
IDX=${3:?need index dir}
OUT=${4:-kb_compare}
THREADS=${5:-8}

mkdir -p "$OUT"
for spec in "v2_reverse:10XV2:reverse" \
            "v2_unstranded:10XV2:unstranded" \
            "v3_reverse:10xv3:reverse" \
            "v3_unstranded:10xv3:unstranded"; do
  lab=${spec%%:*}; rest=${spec#*:}; tech=${rest%%:*}; strand=${rest##*:}
  dir="$OUT/kb_$lab"
  if [ -s "$dir/counts_unfiltered/cells_x_genes.mature.mtx" ]; then
    echo "[skip] $lab already counted"; continue
  fi
  echo "[run ] $lab  (-x $tech --strand $strand)"
  kb count --workflow nac -x "$tech" --strand "$strand" -t "$THREADS" -m 8G -o "$dir" \
    -i "$IDX/index.idx" -g "$IDX/t2g.txt" \
    -c1 "$IDX/cdna_t2c.txt" -c2 "$IDX/nascent_t2c.txt" \
    "$R1" "$R2" > "$OUT/kb_$lab.log" 2>&1
  python - "$dir" <<'PY'
import json, sys
d = sys.argv[1]
ri = json.load(open(f"{d}/run_info.json"))
print("       %d reads, %.1f%% pseudoaligned" % (ri["n_processed"], ri["p_pseudoaligned"]))
PY
done

python $WSsrc/compare_layers.py \
  kb_compare/kb_v2_reverse kb_compare/kb_v2_unstranded \
  kb_compare/kb_v3_reverse kb_compare/kb_v3_unstranded