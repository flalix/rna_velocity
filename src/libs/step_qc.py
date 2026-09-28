"""QC metrics, doublet detection and filtering for the K409 libraries."""
import json
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
import scipy.sparse as sp

adata = sc.read_h5ad("k409_raw_called.h5ad")
adata.var["mt"] = adata.var["gene_name"].str.startswith("MT-")
adata.var["ribo"] = adata.var["gene_name"].str.match(r"^RP[SL]")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)
S = np.asarray(adata.layers["spliced"].sum(1)).ravel()
U = np.asarray(adata.layers["unspliced"].sum(1)).ravel()
Amb = np.asarray(adata.layers["ambiguous"].sum(1)).ravel()
adata.obs["unspliced_frac"] = U / np.maximum(S + U + Amb, 1)

# doublets per library (scrublet via scanpy wrapper)
scores, preds = np.zeros(adata.n_obs), np.zeros(adata.n_obs, bool)
for lib in adata.obs["library"].unique():
    m = (adata.obs["library"] == lib).values
    sub = adata[m].copy()
    sc.pp.scrublet(sub, batch_key=None, random_state=0)
    scores[m] = sub.obs["doublet_score"].values
    preds[m] = sub.obs["predicted_doublet"].values
adata.obs["doublet_score"] = scores
adata.obs["predicted_doublet"] = preds

# thresholds from the observed distributions, applied uniformly across libraries
mt_cut = float(np.percentile(adata.obs["pct_counts_mt"], 98))
mt_cut = min(max(mt_cut, 10.0), 20.0)
gene_lo = 200
umi_lo = 500
keep = ((adata.obs["n_genes_by_counts"] >= gene_lo) & (adata.obs["total_counts"] >= umi_lo)
        & (adata.obs["pct_counts_mt"] <= mt_cut) & (~adata.obs["predicted_doublet"]))
summary = []
for lib in ["tumor_primary", "lymph_node", "blood_PBMC"]:
    m = (adata.obs["library"] == lib).values
    summary.append({"library": lib, "cells_in": int(m.sum()),
                    "fail_genes": int((m & (adata.obs["n_genes_by_counts"] < gene_lo).values).sum()),
                    "fail_umi": int((m & (adata.obs["total_counts"] < umi_lo).values).sum()),
                    "fail_mito": int((m & (adata.obs["pct_counts_mt"] > mt_cut).values).sum()),
                    "doublets": int((m & adata.obs["predicted_doublet"].values).sum()),
                    "cells_out": int((m & keep.values).sum()),
                    "median_umi": float(np.median(adata.obs["total_counts"][m & keep.values])),
                    "median_genes": float(np.median(adata.obs["n_genes_by_counts"][m & keep.values])),
                    "median_pct_mt": float(np.median(adata.obs["pct_counts_mt"][m & keep.values])),
                    "median_unspliced_frac": float(np.median(adata.obs["unspliced_frac"][m & keep.values]))})
qc = pd.DataFrame(summary)
qc.to_csv("qc_summary.csv", index=False)
adata.obs["qc_pass"] = keep.values
adata.write("k409_qc_annotated.h5ad")
filt = adata[keep].copy()
sc.pp.filter_genes(filt, min_cells=3)
filt.write("k409_filtered.h5ad")
json.dump({"mt_cut": mt_cut, "gene_lo": gene_lo, "umi_lo": umi_lo,
           "cells_before": int(adata.n_obs), "cells_after": int(filt.n_obs),
           "genes_after": int(filt.n_vars)}, open("qc_thresholds.json", "w"), indent=1)
print(qc.to_string(index=False))
print("filtered:", filt.shape)
