"""Convert a kb-python `nac` count directory to an scVelo-ready .h5ad.

Usage:  python nac_to_h5ad.py <kb_count_dir> [out.h5ad] [--filtered]

Writes an AnnData with
    X                = mature + ambiguous   (cytoplasmic convention, ~CellRanger)
    layers.spliced   = mature
    layers.unspliced = nascent              <- the layer scVelo needs
    layers.ambiguous = ambiguous
Naming matters: kb's own `--h5ad` output keeps the kb layer names
(mature/nascent/ambiguous), which scv.pp.moments does not recognise.
"""
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio


def _read_mtx(path):
    try:
        return sio.mmread(path, spmatrix=False)
    except TypeError:
        return sio.mmread(path)


def convert(d, out=None, filtered=False):
    sub = "counts_filtered" if filtered else "counts_unfiltered"
    c = os.path.join(d, sub)
    if not os.path.isdir(c):
        sys.exit(f"{c} not found")
    bc = pd.read_csv(os.path.join(c, "cells_x_genes.barcodes.txt"), header=None)[0].astype(str).values
    gid = pd.read_csv(os.path.join(c, "cells_x_genes.genes.txt"), header=None)[0].astype(str).values
    nf = os.path.join(c, "cells_x_genes.genes.names.txt")
    names = pd.read_csv(nf, header=None)[0].astype(str).values if os.path.exists(nf) else gid
    L = {k: _read_mtx(os.path.join(c, f"cells_x_genes.{k}.mtx")).tocsr().astype(np.float32)
         for k in ("mature", "nascent", "ambiguous")}
    A = ad.AnnData(X=(L["mature"] + L["ambiguous"]).tocsr(),
                   obs=pd.DataFrame(index=bc),
                   var=pd.DataFrame({"gene_id": gid, "gene_name": names}, index=names))
    A.layers["spliced"] = L["mature"]
    A.layers["unspliced"] = L["nascent"]
    A.layers["ambiguous"] = L["ambiguous"]
    A.var_names_make_unique()
    out = out or os.path.join(d, "adata_scvelo.h5ad")
    A.write(out)
    tot = sum(float(v.sum()) for v in L.values())
    print(f"{A.n_obs} barcodes x {A.n_vars} genes -> {out}")
    print("  spliced %.0f | unspliced %.0f | ambiguous %.0f | unspliced fraction %.3f"
          % (L["mature"].sum(), L["nascent"].sum(), L["ambiguous"].sum(), L["nascent"].sum() / tot))
    return A


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not args:
        sys.exit(__doc__)
    convert(args[0], args[1] if len(args) > 1 else None, filtered="--filtered" in sys.argv)
