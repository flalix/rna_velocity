"""scFates trajectory on an AnnData with no splicing counts.

Follows 06_3_HDCA_heart_scFates_innervation.py. That script loads the
spliced/unspliced objects but never uses them -- scFates works off `X`
alone -- so the same analysis runs on an ordinary counts object.

One deviation worth recording: in 06_3, `X` is the *spliced* matrix
(`adata = spliced_data.copy()`). Running on total counts changes the
embedding slightly, so the tree will not be numerically identical to
the published Fig. 4D.

    from run_scfates import preprocess, fit_curve, plot_trajectory

    adata = preprocess(ad, n_top_genes=5000, n_neighbors=25)
    fit_curve(adata, root_gene="PENK", nodes=30)
    plot_trajectory(adata, cluster_key="clusters_subset",
                    genes=["SOX10", "FOXD3", "MBP", "MPZ", "ASCL1", "PRPH"],
                    out_dir="figures")
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path


def _repair_cached_rpy2(home):
    """Undo a poisoned rpy2 import.

    `rpy2.rinterface_lib.openrlib` evaluates `R_HOME = get_r_home()` at
    import time and stores None WITHOUT raising. The failure only surfaces
    later, in `rpy2.rinterface.initr()`, as

        RuntimeError: Unable to determine R_HOME.

    So any earlier import of rpy2 (directly, or via scFates) in a process
    that lacked R_HOME leaves a cached module holding None, and setting the
    environment variable afterwards is invisible to it -- the same error
    repeats forever. Writing the resolved path onto the cached module fixes
    it in place; a kernel restart is the alternative.
    """
    mod = sys.modules.get("rpy2.rinterface_lib.openrlib")
    if mod is not None and getattr(mod, "R_HOME", None) is None:
        mod.R_HOME = home


def _ensure_r_home():
    """Point rpy2 at an R installation before scFates imports it.

    scFates imports rpy2, which resolves R_HOME and dlopens libR.so at
    import time. With R in a conda env and Python in a separate venv, a
    Jupyter kernel usually has neither R_HOME set nor R on PATH, and the
    import dies with

        Unable to determine R home: [Errno 2] No such file or directory: 'R'

    Setting the variables in a notebook cell only helps if that cell runs
    before the first rpy2 import anywhere in the process -- doing it here
    means importing this module is enough, whatever the cell order.

    Override the search with R_HOME, or point HDCA_R_BIN at an R binary.
    """
    if os.environ.get("R_HOME"):
        home = os.environ["R_HOME"]
        # R_HOME alone is not enough: scFates decides R is installed with
        # shutil.which("R"), i.e. the executable on PATH, and reports
        # "R installation is necessary for ..." when it is absent -- even
        # though rpy2 itself works fine off R_HOME.
        if shutil.which("R") is None and (Path(home) / "bin" / "R").exists():
            os.environ["PATH"] = f"{Path(home) / 'bin'}{os.pathsep}{os.environ.get('PATH', '')}"
        _repair_cached_rpy2(home)
        return home

    roots = ["miniforge3", "miniconda3", "anaconda3", "mambaforge", "micromamba",
             "conda", ".conda", ".claude-science/conda"]
    candidates = [os.environ.get("HDCA_R_BIN"), shutil.which("R")]
    for root in roots:
        candidates += [str(p) for p in sorted((Path.home() / root).glob("envs/*/bin/R"))]
    candidates += [str(p) for p in sorted(Path("/opt").glob("*conda*/envs/*/bin/R"))]
    for r_bin in filter(None, candidates):
        if not Path(r_bin).exists():
            continue
        try:
            home = subprocess.run([r_bin, "RHOME"], capture_output=True, text=True,
                                  timeout=30).stdout.strip().splitlines()[-1]
        except Exception:
            continue
        if home and Path(home).is_dir():
            os.environ["R_HOME"] = home
            # so rpy2's `R RHOME` fallback and any R subprocess agree
            os.environ["PATH"] = f"{Path(r_bin).parent}{os.pathsep}{os.environ.get('PATH', '')}"
            _repair_cached_rpy2(home)
            return home
    return None


_ensure_r_home()

import matplotlib
import numpy as np
import scanpy as sc
import scFates as scf


def preprocess(adata, n_top_genes=5000, n_neighbors=25, n_comps=50,
               layer=None, copy=True):
    """Preprocessing as in 06_3: log10, cell_ranger HVGs, scale, PCA, kNN.

    `layer` selects the count matrix to start from (e.g. "counts");
    None uses `X` as it stands. `.raw` is set before HVG subsetting, so
    gene-level plots still reach the full gene set.
    """
    a = adata.copy() if copy else adata
    if layer is not None:
        a.X = a.layers[layer].copy()

    sc.pp.filter_genes(a, min_cells=3)
    sc.pp.normalize_total(a)
    sc.pp.log1p(a, base=10)
    sc.pp.highly_variable_genes(a, n_top_genes=n_top_genes, flavor="cell_ranger")
    a.raw = a
    if copy:
        a = a[:, a.var.highly_variable].copy()
    else:
        # `a[:, mask].copy()` would rebind a local name and leave the caller's
        # object un-subset -- silently fitting on every gene. Subset in place.
        a._inplace_subset_var(a.var.highly_variable.values)
    sc.pp.scale(a)
    sc.pp.pca(a, n_comps=n_comps)
    sc.pp.neighbors(a, n_neighbors=n_neighbors)
    print(f"preprocessed: {a.n_obs} cells x {a.n_vars} HVGs")
    return a


def fit_curve(adata, root_gene=None, root_node=None, nodes=30,
              use_rep="X_pca", ndims_rep=2, n_map=100, n_jobs=1, seed=42):
    """Fit the principal curve (ElPiGraph) and compute pseudotime.

    Direction is not inferred from the data -- it comes from the root.
    Pass `root_gene` to root at the node with highest mean expression of
    that gene (06_3 uses PENK), or `root_node` to set it by index. The
    choice fixes the direction of every downstream result, so it needs
    to be a biological claim you can defend, not a default.
    """
    scf.tl.curve(adata, Nodes=nodes, use_rep=use_rep, ndims_rep=ndims_rep)

    if root_gene is not None:
        if root_gene not in adata.var_names and (
            adata.raw is None or root_gene not in adata.raw.var_names
        ):
            raise ValueError(
                f"{root_gene!r} is in neither var_names nor raw.var_names -- "
                "it was filtered out; pick another root marker or raise n_top_genes"
            )
        scf.tl.root(adata, root_gene)
    elif root_node is not None:
        scf.tl.root(adata, root_node)
    else:
        raise ValueError("pass root_gene or root_node -- pseudotime needs a root")

    scf.tl.pseudotime(adata, n_jobs=n_jobs, n_map=n_map, seed=seed)
    t = adata.obs["t"]
    print(f"pseudotime: {float(t.min()):.2f} - {float(t.max()):.2f} "
          f"over {adata.obs['seg'].nunique()} segment(s)")
    return adata


def plot_trajectory(adata, cluster_key=None, genes=(), basis="pca",
                    out_dir=None, cmap="RdBu_r"):
    """Pseudotime, clusters and marker panels on the chosen embedding."""
    if out_dir is not None:
        matplotlib.use("Agg")
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        sc.settings.figdir = out_dir

    plot = getattr(sc.pl, basis)
    saved = []

    def _save(suffix):
        return f"_{suffix}.pdf" if out_dir is not None else None

    plot(adata, color="t", show=out_dir is None, save=_save("pseudotime"))
    saved.append(f"{basis}_pseudotime.pdf")
    if cluster_key is not None:
        plot(adata, color=cluster_key, show=out_dir is None, save=_save("clusters"))
        saved.append(f"{basis}_clusters.pdf")

    genes = [g for g in genes
             if g in adata.var_names or (adata.raw is not None
                                         and g in adata.raw.var_names)]
    if genes:
        plot(adata, color=genes, cmap=cmap, show=out_dir is None,
             save=_save("markers"))
        saved.append(f"{basis}_markers.pdf")

    if out_dir is not None:
        print("wrote:", ", ".join(str(out_dir / s) for s in saved))
    return saved


def test_association(adata, n_jobs=1, fdr=1e-4, a_cutoff=0.3):
    """Genes significantly associated with pseudotime, then fitted.

    Not part of 06_3 -- this is the standard scFates follow-on, and the
    output (`adata.var['signi']`, layer 'fitted') is what feeds
    `scf.pl.trends` heatmaps. The keyword is `fdr_cut`, not `fdr`, on
    scFates 1.2.x.

    Requires rpy2 and R's mgcv -- scFates fits the GAMs in R and raises
    "rpy2 installation is necessary for testing feature association to
    the tree" without them. The project's scFates.yml already pins
    r-base=4.3.2 and rpy2=3.5.11; in a fresh env install them with
    `conda install -c conda-forge rpy2 r-mgcv`. Everything above this
    function is pure Python and needs neither.
    """
    scf.tl.test_association(adata, n_jobs=n_jobs, fdr_cut=fdr, A_cut=a_cutoff)
    n = int(adata.var["signi"].sum())
    print(f"{n} genes associated with pseudotime (fdr<{fdr}, A>{a_cutoff})")
    if n:
        scf.tl.fit(adata, n_jobs=n_jobs)
    return adata
