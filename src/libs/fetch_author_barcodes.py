"""Download the authors' CellRanger cell barcodes for the K409 libraries from GEO.

Writes geo_barcodes/author_cells.json -> {library: [barcode, ...]} with the -1 suffix
stripped, which is what load_qc.py expects. Run once before load_qc.py.

Usage:  python fetch_author_barcodes.py [outdir]
"""
import gzip
import json
import os
import re
import sys
import urllib.request

GSM = {"tumor_primary": "GSM4455935",   # K409 cutaneous primary tumour
       "lymph_node": "GSM4455933",      # K409 involved regional lymph node
       "blood_PBMC": "GSM4455932"}      # K409 peripheral blood


def fetch(gsm):
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/"
    with urllib.request.urlopen(base, timeout=60) as r:
        listing = r.read().decode("utf8", "replace")
    files = sorted(set(re.findall(r'href="([^"]+barcodes[^"]*)"', listing)))
    if not files:
        raise RuntimeError(f"no barcodes file listed for {gsm} at {base}")
    with urllib.request.urlopen(base + files[0], timeout=300) as r:
        raw = r.read()
    barcodes = gzip.decompress(raw).decode().split()
    return files[0], [b.split("-")[0] for b in barcodes]


def main(outdir="geo_barcodes"):
    os.makedirs(outdir, exist_ok=True)
    out = {}
    for lib, gsm in GSM.items():
        fn, bcs = fetch(gsm)
        out[lib] = sorted(set(bcs))
        print(f"{lib:14s} {gsm}  {fn}  {len(bcs)} cells")
    path = os.path.join(outdir, "author_cells.json")
    json.dump(out, open(path, "w"))
    print("wrote", path)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "geo_barcodes")
