#!/usr/bin/env bash
# Verify md5 (if a targets.tsv is present) then kb count (nac) each library.
#
# Usage:   bash count_libs.sh tumor_primary:SRR11492025 lymph_node:SRR11492017 ...
#          several lanes of one library: tumor_primary:SRR11492025,SRR11492026
# Paths:   KBREF=<dir with index.idx,t2g.txt,cdna_t2c.txt,nascent_t2c.txt>  (default ./kbref)
#          FQDIR=<dir with <SRR>_1.fastq.gz / _2.fastq.gz>                  (default ./fastq)
#          TARGETS=<tsv: filename \t bytes \t url \t md5>                   (default $FQDIR/targets.tsv)
#          COUNT_ROOT=<dir for count_<library>/ outputs>                   (default .)
#          THREADS (default 8), MEM (default 16G)
#
# Chemistry is fixed to the GSE148190 5' libraries: -x 10XV2 (16 bp barcode + 10 bp UMI)
# and --strand reverse. Both were verified on the data, not taken from the kit spec.
set -u
KBREF=${KBREF:-kbref}
FQDIR=${FQDIR:-fastq}
TARGETS=${TARGETS:-$FQDIR/targets.tsv}
COUNT_ROOT=${COUNT_ROOT:-.}      # where count_<library>/ dirs are written (load_qc.py reads the same var)
THREADS=${THREADS:-8}
MEM=${MEM:-16G}

for f in index.idx t2g.txt cdna_t2c.txt nascent_t2c.txt; do
  [ -s "$KBREF/$f" ] || { echo "missing $KBREF/$f (build with: kb ref -d human --workflow nac ...)"; exit 1; }
done

for spec in "$@"; do
  lab=${spec%%:*}; srrs=${spec##*:}
  IFS=',' read -r -a runs <<< "$srrs"
  FQ_ARGS=()
  for srr in "${runs[@]}"; do
  for r in 1 2; do
    fq="$FQDIR/${srr}_${r}.fastq.gz"
    [ -s "$fq" ] || { echo "missing $fq"; exit 1; }
    if [ -s "$TARGETS" ]; then
      exp=$(awk -F'\t' -v f="${srr}_${r}.fastq.gz" '$1==f {print $4}' "$TARGETS")
      if [ -n "${exp:-}" ]; then
        obs=$(md5sum "$fq" | cut -d' ' -f1)
        [ "$obs" = "$exp" ] || { echo "MD5 MISMATCH ${srr}_${r}: $obs != $exp"; exit 1; }
        echo "[$(date -u +%H:%M:%S)] md5 ok ${srr}_${r}"
      fi
    fi
  done
  FQ_ARGS+=("$FQDIR/${srr}_1.fastq.gz" "$FQDIR/${srr}_2.fastq.gz")
  done
  if [ -s "$COUNT_ROOT/count_${lab}/counts_unfiltered/cells_x_genes.mature.mtx" ]; then
    echo "[$(date -u +%H:%M:%S)] $lab already counted"; continue
  fi
  echo "[$(date -u +%H:%M:%S)] kb count $lab (${#runs[@]} lane(s): $srrs)"
  kb count -i "$KBREF/index.idx" -g "$KBREF/t2g.txt" \
    -c1 "$KBREF/cdna_t2c.txt" -c2 "$KBREF/nascent_t2c.txt" \
    -x 10XV2 --workflow nac --strand reverse -t "$THREADS" -m "$MEM" -o "$COUNT_ROOT/count_${lab}" \
    "${FQ_ARGS[@]}" > "$COUNT_ROOT/count_${lab}.log" 2>&1
  echo "[$(date -u +%H:%M:%S)] $lab kb exit=$?"
  python - "$COUNT_ROOT/count_${lab}" <<'PY'
import json, sys
ri = json.load(open(f"{sys.argv[1]}/run_info.json"))
print("  processed %d reads, pseudoaligned %.1f%%" % (ri["n_processed"], ri["p_pseudoaligned"]))
PY
done
