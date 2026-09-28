"""Load kb-python nac output for the three K409 libraries, call cells, save h5ad.

Layers: spliced = mature, unspliced = nascent, ambiguous kept separately.
X = mature + ambiguous (cytoplasmic convention, closest to CellRanger).

Expects, relative to the working directory (override the root with COUNT_ROOT):
    count_<library>/counts_unfiltered/cells_x_genes.{mature,nascent,ambiguous}.mtx
    geo_barcodes/author_cells.json      <- run fetch_author_barcodes.py first
"""
import json
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

COUNT_ROOT = os.environ.get("COUNT_ROOT", ".")
BARCODE_JSON = os.environ.get("AUTHOR_CELLS", "geo_barcodes/author_cells.json")
LIBS = {"tumor_primary": "cutaneous primary tumour",
        "lymph_node": "involved regional lymph node",
        "blood_PBMC": "peripheral blood"}


def load_one(lib):
    d = os.path.join(COUNT_ROOT, f"count_{lib}", "counts_unfiltered")
    bc = pd.read_csv(f"{d}/cells_x_genes.barcodes.txt", header=None)[0].astype(str).values
    genes = pd.read_csv(f"{d}/cells_x_genes.genes.txt", header=None)[0].astype(str).values
    names_f = f"{d}/cells_x_genes.genes.names.txt"
    names = pd.read_csv(names_f, header=None)[0].astype(str).values if os.path.exists(names_f) else genes
    layers = {k: sc.read_mtx(f"{d}/cells_x_genes.{k}.mtx").X.astype(np.float32)
              for k in ("mature", "nascent", "ambiguous")}
    A = ad.AnnData(X=(layers["mature"] + layers["ambiguous"]).tocsr(),
                   obs=pd.DataFrame(index=[f"{lib}:{b}" for b in bc]),
                   var=pd.DataFrame({"gene_id": genes, "gene_name": names}, index=names))
    A.layers["spliced"] = layers["mature"].tocsr()
    A.layers["unspliced"] = layers["nascent"].tocsr()
    A.layers["ambiguous"] = layers["ambiguous"].tocsr()
    A.obs["library"] = lib
    A.obs["tissue"] = LIBS[lib]
    A.obs["barcode"] = bc
    A.var_names_make_unique()
    return A


def knee_threshold(counts, floor=500):
    """Inflection of the log-log barcode-rank curve (kneedle-style), with a UMI floor."""
    v = np.sort(counts)[::-1]
    v = v[v > 0]
    x, y = np.log10(np.arange(1, len(v) + 1)), np.log10(v)
    x0, y0, x1, y1 = x[0], y[0], x[-1], y[-1]
    num = np.abs((y1 - y0) * x - (x1 - x0) * y + (x1 * y0 - y1 * x0))
    idx = int(np.argmax(num / np.hypot(y1 - y0, x1 - x0)))
    return max(float(v[idx]), floor), idx, v


if __name__ == "__main__":
    # Cells are the authors' CellRanger-called barcodes, which keeps this analysis
    # comparable with the published one. The log-log knee is computed as well, but only
    # recorded for the barcode-rank figure -- it is stricter, and in this dataset every
    # knee barcode is a subset of the author call.
    if not os.path.exists(BARCODE_JSON):
        raise SystemExit(f"{BARCODE_JSON} not found -- run fetch_author_barcodes.py first")
    author = {k: set(v) for k, v in json.load(open(BARCODE_JSON)).items()}
    adatas, stats, ranks = {}, {}, {}
    for lib in LIBS:
        A = load_one(lib)                      # loaded once; ranks cached for the figure
        tot = np.asarray(A.X.sum(1)).ravel()
        thr, idx, ranked = knee_threshold(tot)
        ranks[lib] = ranked
        keep = np.isin(A.obs["barcode"].values, list(author[lib]))
        S = A[keep].copy()
        layer_tot = sum(float(S.layers[k].sum()) for k in ("spliced", "unspliced", "ambiguous"))
        stats[lib] = {"barcodes_total": int(A.n_obs), "author_cells": len(author[lib]),
                      "author_cells_found": int(keep.sum()), "knee_threshold": thr,
                      "knee_cells": int((tot >= thr).sum()),
                      "median_umi": float(np.median(np.asarray(S.X.sum(1)).ravel())),
                      "reads_unspliced_frac": float(S.layers["unspliced"].sum() / layer_tot)}
        adatas[lib] = S
        print(lib, stats[lib])
    adata = ad.concat(adatas, label=None, index_unique=None)
    adata.var = adatas["tumor_primary"].var.loc[adata.var_names]
    adata.write("k409_raw_called.h5ad")
    json.dump(stats, open("cell_calling.json", "w"), indent=1)
    np.savez("barcode_ranks.npz", **ranks)
    print("combined:", adata.shape)
