"""Pre-flight check for an AnnData before running the scVelo dynamical model.

`scv.pp.moments()` does NOT raise when spliced/unspliced counts are
absent -- it logs "Skipping moments, because un/spliced counts were not
found." and returns, leaving no Ms/Mu layers. `scv.tl.recover_dynamics()`
then falls back to the raw path (`use_raw=True`, set because Ms/Mu are
missing) and dies with `KeyError: 'unspliced'`.

Run `check_velocity_inputs(ad)` before `scv.pp.moments` to catch that,
and again after, to confirm Ms/Mu were actually written.
"""

import numpy as np

RAW_LAYERS = ("spliced", "unspliced")
MOMENT_LAYERS = ("Ms", "Mu")


def _is_integral(x, n=10000):
    """Do the stored values look like raw counts rather than lognorm data?"""
    from scipy.sparse import issparse

    v = np.asarray(x[:n].data if issparse(x) else np.ravel(np.asarray(x))[:n])
    v = v[np.isfinite(v)]
    return bool(v.size == 0 or np.allclose(v, np.round(v)))


def check_velocity_inputs(adata, verbose=True):
    """Return a dict of findings; print a per-check report when verbose."""
    layers = [k for k in adata.layers.keys() if k is not None]
    report = {
        "n_obs": adata.n_obs,
        "n_vars": adata.n_vars,
        "layers": layers,
        "has_raw_layers": all(k in layers for k in RAW_LAYERS),
        "has_moments": all(k in layers for k in MOMENT_LAYERS),
        "has_neighbors": "neighbors" in adata.uns or "connectivities" in adata.obsp,
        "embeddings": [k for k in adata.obsm.keys() if k.startswith("X_")],
        "is_view": adata.is_view,
    }

    if report["has_raw_layers"]:
        u, s = adata.layers["unspliced"], adata.layers["spliced"]
        tot_u = float(u.sum())
        tot_s = float(s.sum())
        report["unspliced_fraction"] = tot_u / (tot_u + tot_s) if tot_u + tot_s else 0.0
        report["counts_look_integral"] = _is_integral(s) and _is_integral(u)

    if verbose:
        def line(ok, msg):
            print(f"  [{'OK ' if ok else 'FAIL'}] {msg}")

        print(f"AnnData {adata.n_obs} cells x {adata.n_vars} genes"
              f"{' (VIEW -- call .copy() first)' if report['is_view'] else ''}")
        print(f"  layers: {layers or 'none'}")
        line(report["has_raw_layers"],
             "spliced/unspliced layers present (required by scv.pp.moments)")
        if report["has_raw_layers"]:
            frac = report["unspliced_fraction"]
            line(0.05 <= frac <= 0.60,
                 f"unspliced fraction = {frac:.1%} (typical 10-40%; "
                 "far outside that suggests a bad quantification or a mis-merge)")
            if not report["has_moments"]:
                # scv.pp.moments normalises the layers in place, so this
                # check only means anything before it has run.
                line(report["counts_look_integral"],
                     "layer values look like raw counts, not log-normalised data")
        line(report["has_moments"],
             "Ms/Mu present (written by scv.pp.moments; needed by recover_dynamics)")
        line(report["has_neighbors"], "neighbor graph present")
        line(bool(report["embeddings"]),
             f"embedding for plotting: {report['embeddings'] or 'none'}")

    return report


def merge_loom(adata, ldata, strip_loom_barcodes=True, min_overlap=0.5):
    """Merge velocyto/kallisto loom counts into `adata` by barcode.

    velocyto writes obs_names like 'sample:AAACCCAAGAAACCATx'; Cell Ranger
    AnnData usually has 'AAACCCAAGAAACCAT-1'. `strip_loom_barcodes`
    applies that transformation. Check the printed overlap -- a low
    number means the two objects use different barcode conventions
    (per-sample prefixes, -1/-2 suffixes) and need manual harmonisation
    before the merge, not after.
    """
    ldata = ldata.copy()
    ldata.var_names_make_unique()
    if strip_loom_barcodes:
        ldata.obs_names = [
            bc.split(":")[-1].rstrip("x") + "-1" for bc in ldata.obs_names
        ]
    overlap = len(set(adata.obs_names) & set(ldata.obs_names))
    frac = overlap / adata.n_obs if adata.n_obs else 0.0
    print(f"barcode overlap: {overlap}/{adata.n_obs} cells ({frac:.1%})")
    if frac < min_overlap:
        raise ValueError(
            f"only {frac:.1%} of cells matched -- harmonise obs_names first "
            f"(adata: {list(adata.obs_names[:2])}, loom: {list(ldata.obs_names[:2])})"
        )
    import scvelo as scv
    return scv.utils.merge(adata, ldata)
