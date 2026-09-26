"""Build a velocity-ready AnnData directly from the HDCA per-sample .h5 files.

Follows the project pipeline (06_1 -> 06_2 -> 06_4), but reads the
`shoji/Spliced` and `shoji/Unspliced` datasets in Python instead of
going Seurat -> SaveH5Seurat -> Convert(). Same inputs, same result,
one step: `data/h5_files/*.h5` -> AnnData with `spliced`/`unspliced`
layers.

    ad = build_from_shoji(
        "data/h5_files",
        metadata_tsv="data/AnnData/metadata_IN_selection.tsv",
        umap_csv="data/AnnData/UMAP_IN_selection.csv",
    )
"""

from pathlib import Path

import anndata as ad_mod
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, vstack

SHOJI = {"spliced": "shoji/Spliced", "unspliced": "shoji/Unspliced",
         "cellid": "shoji/Cellid", "gene": "shoji/Gene"}


def _load_anndata(src):
    if not isinstance(src, (str, Path)):
        return src
    src = Path(src)
    if src.suffix == ".loom":
        return ad_mod.read_loom(src)
    if src.suffix in (".h5", ".hdf5") and h5py.File(src, "r").get("shoji") is not None:
        return read_shoji_h5(src)
    return ad_mod.read_h5ad(src)


def attach_spliced_unspliced(adata, spliced, unspliced=None, subset=True,
                             min_cell_overlap=0.5, min_gene_overlap=0.5):
    """Add spliced/unspliced layers to an AnnData you already have.

    Use this when your working object carries the annotations and
    embeddings you want to keep, and you only need the splicing counts
    bolted on. Matching is by obs_names and var_names, so the sources
    may be in any cell/gene order and may contain extra cells or genes.

    `spliced` / `unspliced` may each be an AnnData or a path to .h5ad,
    .loom, or a shoji .h5. Pass a single object to `spliced` if it
    already carries both layers (a velocyto loom does).

    With subset=True the returned object is restricted to the shared
    cells and genes -- velocity can only be computed there anyway.
    """
    s = _load_anndata(spliced)
    if unspliced is None:
        if not {"spliced", "unspliced"} <= set(s.layers.keys()):
            raise ValueError("single source must carry spliced+unspliced layers")
        u_mat, s_mat, u_obj = s.layers["unspliced"], s.layers["spliced"], s
    else:
        u_obj = _load_anndata(unspliced)
        s_mat, u_mat = s.X, u_obj.X
        if not s.obs_names.equals(u_obj.obs_names):
            raise ValueError("spliced and unspliced sources have different cells")

    cells = adata.obs_names.intersection(s.obs_names)
    genes = adata.var_names.intersection(s.var_names)
    c_frac = len(cells) / adata.n_obs
    g_frac = len(genes) / adata.n_vars
    print(f"shared cells: {len(cells)}/{adata.n_obs} ({c_frac:.1%})  "
          f"shared genes: {len(genes)}/{adata.n_vars} ({g_frac:.1%})")
    if c_frac < min_cell_overlap:
        raise ValueError(
            f"only {c_frac:.1%} of cells matched -- barcode conventions differ "
            f"(target: {list(adata.obs_names[:2])}, source: {list(s.obs_names[:2])})"
        )
    if g_frac < min_gene_overlap:
        raise ValueError(
            f"only {g_frac:.1%} of genes matched -- gene IDs differ "
            f"(target: {list(adata.var_names[:2])}, source: {list(s.var_names[:2])})"
        )

    out = adata[cells, genes].copy() if subset else adata.copy()
    src_s = s[cells, genes]
    src_u = u_obj[cells, genes] if unspliced is not None else s[cells, genes]
    out.layers["spliced"] = csr_matrix(
        src_s.X if unspliced is not None else src_s.layers["spliced"]
    )
    out.layers["unspliced"] = csr_matrix(
        src_u.X if unspliced is not None else src_u.layers["unspliced"]
    )
    return out


def preprocess_for_dynamical(adata, n_top_genes=2000, n_neighbors=25, n_pcs=30):
    """Preprocessing equivalent to 06_4, rewritten for scVelo 0.3.4.

    06_4 was written against scvelo 0.3.0 (the pin in scVelo.yml). On
    0.3.4 `scv.pp.filter_genes_dispersion` and `scv.pp.log1p` no longer
    exist, and `scv.pp.filter_and_normalize` no longer accepts
    `n_top_genes` (it now only filters and normalises; passing it raises
    `TypeError: normalize_per_cell() got an unexpected keyword argument
    'n_top_genes'`). The scanpy equivalents are used here instead.

    Leaves the object ready for `scv.tl.recover_dynamics`.
    """
    import scanpy as sc
    import scvelo as scv

    scv.pp.filter_genes(adata, min_cells=3)
    scv.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=True)
    sc.pp.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors)
    scv.pp.moments(adata, n_pcs=None, n_neighbors=None)  # reuse the graph above
    return adata


def _decode(arr):
    arr = np.asarray(arr).ravel()
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr])


def _as_cells_by_genes(mat, n_cells, n_genes):
    """Orient a shoji matrix to cells x genes.

    rhdf5 and h5py disagree on axis order for the same file, so decide
    from the shape rather than assuming. Square matrices are ambiguous
    and rejected.
    """
    if mat.shape == (n_cells, n_genes) == (n_genes, n_cells):
        raise ValueError("square matrix -- orientation ambiguous, set it manually")
    if mat.shape == (n_cells, n_genes):
        return mat
    if mat.shape == (n_genes, n_cells):
        return mat.T
    raise ValueError(
        f"matrix {mat.shape} matches neither {(n_cells, n_genes)} nor "
        f"{(n_genes, n_cells)} (Cellid={n_cells}, Gene={n_genes})"
    )


def read_shoji_h5(path):
    """One sample file -> AnnData (X = spliced) with spliced/unspliced layers."""
    with h5py.File(path, "r") as f:
        cells = _decode(f[SHOJI["cellid"]][:])
        genes = _decode(f[SHOJI["gene"]][:])
        spl = _as_cells_by_genes(f[SHOJI["spliced"]][:], len(cells), len(genes))
        uns = _as_cells_by_genes(f[SHOJI["unspliced"]][:], len(cells), len(genes))

    a = ad_mod.AnnData(
        X=csr_matrix(spl.astype("float32")),
        obs=pd.DataFrame(index=pd.Index(cells, dtype=str)),
        var=pd.DataFrame(index=pd.Index(genes, dtype=str)),
    )
    a.layers["spliced"] = a.X.copy()
    a.layers["unspliced"] = csr_matrix(uns.astype("float32"))
    return a


def build_from_shoji(h5_dir, pattern="*.h5", metadata_tsv=None, umap_csv=None,
                     umap_key="X_umap"):
    """Concatenate all sample .h5 files, then optionally subset + annotate.

    `metadata_tsv` is one of the tables written by 06_1 (metadata.tsv,
    metadata_CM/EN/FB/IN.tsv, metadata_IN_selection.tsv); its index
    selects the cells, exactly as 06_4 does. `umap_csv` is the matching
    UMAP export, attached as `obsm[umap_key]`.
    """
    files = sorted(Path(h5_dir).glob(pattern))
    if not files:
        raise FileNotFoundError(f"no files matching {pattern} in {h5_dir}")

    parts = [read_shoji_h5(f) for f in files]
    genes = parts[0].var_names
    for p, f in zip(parts[1:], files[1:]):
        if not p.var_names.equals(genes):
            raise ValueError(f"gene order differs in {f.name} -- align var_names first")

    adata = ad_mod.AnnData(
        X=vstack([p.X for p in parts], format="csr"),
        obs=pd.concat([p.obs for p in parts]),
        var=parts[0].var.copy(),
        layers={
            "spliced": vstack([p.layers["spliced"] for p in parts], format="csr"),
            "unspliced": vstack([p.layers["unspliced"] for p in parts], format="csr"),
        },
    )
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    print(f"{len(files)} files -> {adata.n_obs} cells x {adata.n_vars} genes")

    if metadata_tsv is not None:
        meta = pd.read_table(metadata_tsv, sep="\t", index_col=0)
        keep = adata.obs_names.intersection(meta.index)
        if len(keep) < 0.5 * len(meta):
            raise ValueError(
                f"only {len(keep)}/{len(meta)} metadata cells found in the h5 files "
                f"-- barcode conventions differ (h5: {list(adata.obs_names[:2])}, "
                f"meta: {list(meta.index[:2])})"
            )
        adata = adata[keep].copy()
        adata.obs = adata.obs.join(meta)
        print(f"subset to metadata: {adata.n_obs} cells, "
              f"{meta.shape[1]} annotation columns")

    if umap_csv is not None:
        umap = pd.read_csv(umap_csv, index_col=0)
        missing = adata.obs_names.difference(umap.index)
        if len(missing):
            raise ValueError(f"{len(missing)} cells absent from {umap_csv}")
        adata.obsm[umap_key] = umap.loc[adata.obs_names].to_numpy()

    return adata
