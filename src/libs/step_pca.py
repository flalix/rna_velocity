"""Normalisation, HVG selection, PCA and PC-covariate diagnostics (pipeline step 3).

Batch correction is a DELIBERATE decision here, not an automatic one -- see
APPLY_HARMONY below before changing it.
"""
import json
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import spearmanr

# ---------------------------------------------------------------------------
# The three libraries are three anatomical sites from ONE patient on identical
# 10x 5' chemistry, so "library" is confounded with tissue of origin -- the
# biology under study. PC1 loads on ribosomal content (rho ~ -0.86) and separates
# blood (rho ~ -0.64): that is the naive-to-activated axis, not a batch effect.
# Harmonising over "library" would delete it, so correction is OFF.
# Set True only when the batch variable is genuinely technical (e.g. several
# patients or several sequencing runs of the same tissue).
APPLY_HARMONY = False
BATCH_KEY = "library"
# ---------------------------------------------------------------------------

adata = sc.read_h5ad("k409_filtered.h5ad")
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3",
                            layer="counts", batch_key=BATCH_KEY)
sc.pp.pca(adata, n_comps=50, mask_var="highly_variable", svd_solver="arpack")

# --- PC vs covariate correlation ---
covs = pd.DataFrame({
    "total_counts": adata.obs["total_counts"].values,
    "n_genes": adata.obs["n_genes_by_counts"].values,
    "pct_mito": adata.obs["pct_counts_mt"].values,
    "pct_ribo": adata.obs["pct_counts_ribo"].values,
    "unspliced_frac": adata.obs["unspliced_frac"].values,
    "doublet_score": adata.obs["doublet_score"].values,
}, index=adata.obs_names)
libs = (adata.obs[BATCH_KEY].cat.categories if hasattr(adata.obs[BATCH_KEY], "cat")
        else adata.obs[BATCH_KEY].unique())
for lib in libs:
    covs[f"is_{lib}"] = (adata.obs[BATCH_KEY] == lib).astype(int).values

P = adata.obsm["X_pca"][:, :30]
rho = np.zeros((P.shape[1], covs.shape[1]))
for j, c in enumerate(covs.columns):
    r, _ = spearmanr(P, covs[c].values)
    rho[:, j] = r[:P.shape[1], -1]
corr = pd.DataFrame(rho, index=[f"PC{i+1}" for i in range(P.shape[1])], columns=covs.columns)
corr.to_csv("pca_covariate_correlation.csv")
var_ratio = adata.uns["pca"]["variance_ratio"]
pd.DataFrame({"pc": np.arange(1, len(var_ratio) + 1), "variance_ratio": var_ratio,
              "cumulative": np.cumsum(var_ratio)}).to_csv("pca_variance.csv", index=False)

tech = corr[["total_counts", "n_genes", "pct_mito", "pct_ribo", "unspliced_frac", "doublet_score"]].abs()
batch = corr[[c for c in corr.columns if c.startswith("is_")]].abs()
diag = {"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars),
        "n_hvg": int(adata.var["highly_variable"].sum()),
        "var_explained_30pc": float(np.sum(var_ratio[:30])),
        "max_abs_rho_technical": float(tech.values.max()),
        "worst_technical_pc": str(tech.max(axis=1).idxmax()),
        "worst_technical_cov": str(tech.max(axis=0).idxmax()),
        "max_abs_rho_library": float(batch.values.max()),
        "worst_library_pc": str(batch.max(axis=1).idxmax()),
        "library_rho_exceeds_0.3_in_top10_pcs": bool(batch.iloc[:10].values.max() > 0.3)}

if APPLY_HARMONY:
    # harmonypy directly: scanpy's harmony_integrate wrapper mis-shapes Z_corr on
    # recent harmonypy versions (ValueError: incorrect shape for X_pca_harmony).
    import harmonypy as hm
    ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, [BATCH_KEY], max_iter_harmony=20)
    Z = np.asarray(ho.Z_corr)
    if Z.shape[0] != adata.n_obs:
        Z = Z.T
    assert Z.shape[0] == adata.n_obs, f"harmony returned {Z.shape}, expected {adata.n_obs} rows"
    adata.obsm["X_pca_harmony"] = Z
    diag["representation"] = "X_pca_harmony"
else:
    diag["representation"] = "X_pca"
diag["harmony_applied"] = bool(APPLY_HARMONY)
diag["batch_correction_rationale"] = (
    "Not applied: the three libraries are three anatomical sites from one patient on "
    "identical chemistry, so the batch variable is confounded with the biology (PC1 = "
    "ribosomal/naive-activated axis). Correcting it would remove the signal under study."
    if not APPLY_HARMONY else f"Harmony applied over {BATCH_KEY}.")

adata.uns["representation"] = diag["representation"]
adata.write("k409_pca.h5ad")
json.dump(diag, open("pca_diagnostics.json", "w"), indent=1)
print(json.dumps(diag, indent=1))
print(corr.iloc[:8].round(2).to_string())
