"""Assemble the export bundle written by convert_seurat_to_h5ad.R into h5ad.

    python assemble_h5ad.py <bundle_dir> <out.h5ad>

or from a notebook:

    from assemble_h5ad import assemble
    adata = assemble("data/bundles/innervation", "data/h5ad/innervation.h5ad")

Counts land in both `X` and `layers['counts']`; Seurat reductions become
`obsm['X_<name>']` so scanpy and scFates find them under the usual keys.
"""

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix


def assemble(bundle_dir, out_h5ad=None):
    d = Path(bundle_dir)
    manifest = json.loads((d / "manifest.json").read_text())

    # Matrix Market from Seurat is features x cells; AnnData wants cells x features.
    X = csr_matrix(mmread(d / "counts.mtx").T.tocsr(), dtype=np.float32)
    barcodes = pd.read_csv(d / "barcodes.csv")["barcode"].astype(str)
    features = pd.read_csv(d / "features.csv")["feature"].astype(str)
    assert X.shape == (len(barcodes), len(features)), (
        f"matrix {X.shape} vs {len(barcodes)} barcodes x {len(features)} features"
    )

    obs = pd.read_csv(d / "metadata.csv", index_col=0)
    obs.index = obs.index.astype(str)
    assert obs.index.equals(pd.Index(barcodes)), "metadata rows are not in matrix order"

    adata = ad.AnnData(
        X=X,
        obs=obs,
        var=pd.DataFrame(index=pd.Index(features, name=None)),
    )
    adata.layers["counts"] = adata.X.copy()

    for f in sorted(d.glob("reduction_*.csv")):
        name = f.stem.replace("reduction_", "")
        emb = pd.read_csv(f, index_col=0)
        emb.index = emb.index.astype(str)
        assert emb.index.equals(adata.obs_names), f"{f.name} rows are not in matrix order"
        key = name if name.startswith("X_") else f"X_{name}"
        adata.obsm[key] = emb.to_numpy(dtype=np.float32)

    # Categorical-ise low-cardinality string columns so scanpy plots colour them.
    for c in adata.obs.columns:
        if adata.obs[c].dtype == object and adata.obs[c].nunique() <= 200:
            adata.obs[c] = adata.obs[c].astype("category")

    integral = bool(np.allclose(adata.X.data[:100000], np.round(adata.X.data[:100000])))
    print(f"{manifest['source']}: {adata.n_obs} cells x {adata.n_vars} genes | "
          f"layer={manifest['layer']} | integer counts={integral} | "
          f"obsm={list(adata.obsm)} | obs cols={adata.obs.shape[1]}")

    if out_h5ad is not None:
        Path(out_h5ad).parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(out_h5ad, compression="gzip")
        print(f"wrote {out_h5ad}")
    return adata


if __name__ == "__main__":
    assemble(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
