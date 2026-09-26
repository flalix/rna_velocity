"""Fit one compartment end to end, in its own process.

    python src/libs/run_compartment.py cardiomyocytes02 [n_jobs] [n_map]

Run as a SCRIPT, not in the notebook kernel: scFates' GAM step forks one
embedded R per worker and each worker receives a copy of the expression
matrix, so peak memory is roughly n_jobs x (cells x HVGs x 4 bytes) on top
of the object. At 20k cells that OOM-kills a kernel holding other state;
a subprocess death costs nothing.

Writes: results/figures/trends_<name>.png
        results/tables/<name>_gene_modules.csv
        data/h5ad/<name>_fitted.h5ad
"""

from __future__ import annotations

import contextlib
import io
import json
import resource
import sys
import time
import warnings
from pathlib import Path

warnings.simplefilter("ignore")

ROOT_SPEC = {
    "cardiomyocytes02": ("tip", ("Prol_CM", "Immat_CM"), "Cardiomyocytes"),
    "endothelial": ("tip", ("Prol_EC", "Endoc_EC"), "Endothelium"),
    "fibroblasts": ("node", {"EPDC_1", "EPDC_2", "Prol_FB_1", "Prol_FB_2"}, "Fibroblasts"),
    "innervation": ("tip", ("SCP",), "Innervation"),
}


def _rss():
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 1)


def run(name, n_jobs=4, n_map=100, n_top_genes=5000, nodes=20, k=6):
    import matplotlib
    matplotlib.use("Agg")
    import numpy as np
    import scFates as scf

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from libs.assemble_h5ad import assemble
    from libs.paths import BUNDLES, FIGURES, H5AD, TABLES
    from libs.plot_fates import progenitor_node, progenitor_tip
    from libs.trends_figure import trends_figure
    import libs.run_scfates as rs

    mode, spec, title = ROOT_SPEC[name]
    t0 = time.time()

    adata = assemble(BUNDLES / name)
    adata = adata[~adata.obs["cell.types"].astype(str).str.contains("_excl_")].copy()
    adata.obs["cell.types"] = adata.obs["cell.types"].astype(str).astype("category")
    with contextlib.redirect_stdout(io.StringIO()):
        adata = rs.preprocess(adata, n_top_genes=n_top_genes, n_neighbors=25,
                              layer="counts", copy=True)
    # .raw is only needed for gene-level plots; it is a second full-size matrix
    # the GAM workers would otherwise carry. trends_figure reads layers['fitted'].
    adata.raw = None
    assert adata.n_vars <= n_top_genes, f"HVG subset failed: {adata.n_vars} genes"
    print(f"[{name}] preprocessed {adata.shape} | rss {_rss()} GB | {time.time() - t0:.0f}s", flush=True)

    with contextlib.redirect_stdout(io.StringIO()):
        scf.tl.tree(adata, Nodes=nodes, use_rep="X_harmony_subset", ndims_rep=10,
                    method="ppt", ppt_sigma=0.1, ppt_lambda=100, seed=42)
        root = (progenitor_tip(adata, spec) if mode == "tip" else progenitor_node(adata, spec))[0]
        scf.tl.root(adata, root)
        try:
            scf.tl.pseudotime(adata, n_jobs=max(n_jobs * 3, 8), n_map=n_map, seed=42)
        except IndexError:   # empty milestone -> colour step only; numbers are assigned
            adata.uns["seg_colors"] = ["#999999"] * len(adata.uns["graph"]["pp_seg"])
    print(f"[{name}] tree+pseudotime | root {root} | rss {_rss()} GB | {time.time() - t0:.0f}s", flush=True)

    with contextlib.redirect_stdout(io.StringIO()):
        scf.tl.test_association(adata, n_jobs=n_jobs, fdr_cut=1e-4, A_cut=0.3)
    print(f"[{name}] association {int(adata.var['signi'].sum())} genes | rss {_rss()} GB | "
          f"{time.time() - t0:.0f}s", flush=True)

    with contextlib.redirect_stdout(io.StringIO()):
        scf.tl.fit(adata, n_jobs=n_jobs)
    for col in ("seg", "milestones", "edge"):
        adata.obs[col] = adata.obs[col].astype(str).astype("category")

    _, mod = trends_figure(adata, title, str(FIGURES / f"trends_{name}.png"), k=k)
    mod.rename("module").to_frame().join(adata.var[["A", "fdr"]]).to_csv(
        TABLES / f"{name}_gene_modules.csv")
    adata.uns.pop("epg", None)
    adata.write_h5ad(H5AD / f"{name}_fitted.h5ad", compression="gzip")

    g = adata.uns["graph"]
    out = dict(name=name, cells=int(adata.n_obs), assoc=int(mod.size), root=int(g["root"]),
               tips=len(np.ravel(g["tips"])), forks=len(np.ravel(g["forks"])),
               modules=mod.value_counts().sort_index().to_dict(),
               peak_rss_GB=_rss(), minutes=round((time.time() - t0) / 60, 1))
    print(json.dumps(out), flush=True)
    return out


if __name__ == "__main__":
    run(sys.argv[1],
        n_jobs=int(sys.argv[2]) if len(sys.argv) > 2 else 4,
        n_map=int(sys.argv[3]) if len(sys.argv) > 3 else 100)
