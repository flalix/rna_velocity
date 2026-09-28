"""Summarise kb-python `nac` count directories: mature / nascent / ambiguous totals.

Usage:  python compare_layers.py kb_v2_reverse kb_v2_unstranded kb_v3_reverse ...

Each argument is a `kb count -o` output directory (must contain
counts_unfiltered/cells_x_genes.{mature,nascent,ambiguous}.mtx). Percentages are
reported against the directory whose name contains 'v2_reverse' if present,
otherwise against the first one given. Writes kb_chemistry_comparison.csv.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import scipy.io as sio


def summarise(d):
    c = os.path.join(d, "counts_unfiltered")
    layers = {k: sio.mmread(os.path.join(c, f"cells_x_genes.{k}.mtx")).tocsr()
              for k in ("mature", "nascent", "ambiguous")}
    umi = np.asarray((layers["mature"] + layers["ambiguous"]).sum(1)).ravel()
    total = sum(float(v.sum()) for v in layers.values())
    row = {"setting": os.path.basename(d.rstrip("/")),
           "mature": int(layers["mature"].sum()),
           "nascent": int(layers["nascent"].sum()),
           "ambiguous": int(layers["ambiguous"].sum()),
           "unspliced_fraction": round(float(layers["nascent"].sum()) / total, 4) if total else float("nan"),
           "barcodes_ge500_umi": int((umi >= 500).sum())}
    ri = os.path.join(d, "run_info.json")
    if os.path.exists(ri):
        info = json.load(open(ri))
        row["reads"] = int(info.get("n_processed", 0))
        row["pct_pseudoaligned"] = round(float(info.get("p_pseudoaligned", float("nan"))), 1)
    return row


def main(dirs, out_csv="kb_chemistry_comparison.csv"):
    df = pd.DataFrame([summarise(d) for d in dirs])
    ref_rows = df[df.setting.str.contains("v2_reverse")]
    ref = ref_rows.iloc[0] if len(ref_rows) else df.iloc[0]
    for col in ("mature", "nascent"):
        df[f"{col}_vs_ref_pct"] = (100 * (df[col] / ref[col] - 1)).round(2)
    df = df.sort_values("mature", ascending=False)
    df.to_csv(out_csv, index=False)
    print(df.to_string(index=False))
    print(f"\nreference = {ref['setting']}  ->  {out_csv}")
    return df


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
