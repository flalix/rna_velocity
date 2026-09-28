#!/usr/bin/env bash
# Fetch every file in a manifest, skipping the ones whose md5 already matches.
# Files are fetched one at a time; range_fetch.sh already uses N parallel
# connections per file, and more concurrent files does not add throughput.
#
# Usage:  cd <fastq dir> && bash /path/to/fetch_all.sh <manifest.tsv> [parts]
# The manifest is  filename \t bytes \t url \t md5  (subset it to the lanes you want).
set -u
SRC=$(cd "$(dirname "$0")" && pwd)
MANIFEST=$(cd "$(dirname "${1:?manifest tsv}")" && pwd)/$(basename "$1")
PARTS=${2:-8}
ok=0; got=0; bad=0
while IFS=$'\t' read -r fn bytes url md5; do
  [ -z "${fn:-}" ] && continue
  if [ -f "$fn" ] && [ "$(stat -c %s "$fn")" = "$bytes" ] \
     && [ "$(md5sum "$fn" | cut -d' ' -f1)" = "$md5" ]; then
    echo "[ok]   $fn"; ok=$((ok+1)); continue
  fi
  echo "[get]  $fn ($(( bytes / 1000000 )) MB)"
  if MANIFEST="$MANIFEST" bash "$SRC/range_fetch.sh" "$fn" "$PARTS"; then got=$((got+1)); else bad=$((bad+1)); fi
done < "$MANIFEST"
echo "verified=$ok fetched=$got failed=$bad"
[ "$bad" -eq 0 ]
