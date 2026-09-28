#!/usr/bin/env bash
# Download ONE file from ENA with N parallel byte-range requests, then verify md5.
# Single-stream throughput from ENA is ~0.7-1.5 MB/s; 8 ranges give ~8-9 MB/s.
#
# Usage:   bash range_fetch.sh <filename> [parts]
# Env:     MANIFEST  tsv of  filename \t bytes \t url \t md5   (default ./targets.tsv)
#          DRYRUN=1  print the byte ranges and exit without transferring
#
# Runs in the CURRENT directory (cd to the FASTQ dir first) -- it does not write
# next to the script.
set -u
fn=${1:?filename from the manifest}
nparts=${2:-8}
MANIFEST=${MANIFEST:-targets.tsv}
[ -s "$MANIFEST" ] || { echo "manifest not found: $MANIFEST"; exit 2; }

read -r exp url md5 < <(awk -F'\t' -v f="$fn" '$1==f {print $2, $3, $4}' "$MANIFEST")
[ -n "${exp:-}" ] || { echo "$fn not listed in $MANIFEST"; exit 2; }

chunk=$(( (exp + nparts - 1) / nparts ))
echo "[$(date -u +%H:%M:%S)] $fn: $exp bytes in $nparts parts of $chunk"

if [ "${DRYRUN:-0}" = "1" ]; then
  tot=0
  for i in $(seq 0 $((nparts-1))); do
    s=$(( i * chunk )); e=$(( s + chunk - 1 )); [ "$e" -ge "$exp" ] && e=$(( exp - 1 ))
    [ "$s" -ge "$exp" ] && continue
    echo "  part $i: bytes $s-$e  ($(( e - s + 1 )) B)"
    tot=$(( tot + e - s + 1 ))
  done
  echo "  covered $tot of $exp bytes"
  [ "$tot" -eq "$exp" ] || { echo "  RANGE ARITHMETIC ERROR"; exit 1; }
  exit 0
fi

rm -f "$fn"
for i in $(seq 0 $((nparts-1))); do
  s=$(( i * chunk )); e=$(( s + chunk - 1 )); [ "$e" -ge "$exp" ] && e=$(( exp - 1 ))
  [ "$s" -ge "$exp" ] && continue
  (
    part="${fn}.part${i}"; want=$(( e - s + 1 ))
    for attempt in $(seq 1 100); do
      have=0; [ -f "$part" ] && have=$(stat -c %s "$part")
      [ "$have" -ge "$want" ] && break
      curl -sS -o "$part" --range "$(( s + have ))-${e}" --max-time 1800 --connect-timeout 30 \
           $( [ "$have" -gt 0 ] && echo --append ) "$url" || true
    done
  ) &
done
wait
cat $(for i in $(seq 0 $((nparts-1))); do p="${fn}.part${i}"; [ -f "$p" ] && echo "$p"; done) > "$fn"
obs=$(md5sum "$fn" | cut -d' ' -f1)
if [ "$obs" = "$md5" ]; then
  echo "[$(date -u +%H:%M:%S)] $fn OK ($(stat -c %s "$fn") bytes)"; rm -f "${fn}".part*; exit 0
else
  echo "[$(date -u +%H:%M:%S)] $fn BAD observed=$obs expected=$md5 (parts kept)"; exit 1
fi
