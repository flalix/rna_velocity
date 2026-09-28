"""scVelo: moments, stochastic + dynamical velocity, confidence, latent time, drivers (step 8).

Runs single-process (n_jobs=1, threading backend, no progress bar) because the sandbox
blocks the AF_UNIX sockets multiprocessing.Manager needs.
"""
import json
import numpy as np, pandas as pd, scanpy as sc, scvelo as scv

adata = sc.read_h5ad("k409_clustered.h5ad")
rep = adata.uns["representation"]
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)  # reuse the 30-PC, 15-NN graph from step_cluster

# baseline: stochastic (second-moment steady-state)
scv.tl.velocity(adata, mode="stochastic")
scv.tl.velocity_graph(adata, n_jobs=1, show_progress_bar=False, backend="threading")
scv.tl.velocity_confidence(adata)
scv.tl.velocity_pseudotime(adata)
adata.obs["velocity_confidence_stochastic"] = adata.obs["velocity_confidence"].values
adata.obs["velocity_length_stochastic"] = adata.obs["velocity_length"].values

# dynamical model on HVGs -> latent time + gene-level likelihoods
hv = adata.var["highly_variable"].values
D = adata[:, hv].copy()
scv.tl.recover_dynamics(D, n_jobs=1, max_iter=20, show_progress_bar=False, backend="threading")
scv.tl.velocity(D, mode="dynamical")
scv.tl.velocity_graph(D, n_jobs=1, show_progress_bar=False, backend="threading")
scv.tl.velocity_confidence(D)
scv.tl.latent_time(D)
scv.tl.paga(D, groups="annotation", use_time_prior="latent_time")

for c in ["latent_time", "velocity_confidence", "velocity_length", "root_cells", "end_points"]:
    if c in D.obs:
        adata.obs[f"dyn_{c}"] = D.obs[c].values
for c in ["fit_likelihood", "fit_alpha", "fit_beta", "fit_gamma", "fit_t_", "fit_scaling"]:
    if c in D.var:
        adata.var.loc[D.var_names, c] = D.var[c].values

drivers = (D.var[["gene_name", "fit_likelihood", "fit_alpha", "fit_beta", "fit_gamma"]]
           .dropna(subset=["fit_likelihood"]).sort_values("fit_likelihood", ascending=False))
drivers.head(50).to_csv("velocity_drivers.csv", index=False)

# per-cluster velocity summary in the dynamical model
summ = (pd.DataFrame({"annotation": D.obs["annotation"].values,
                      "latent_time": D.obs["latent_time"].values,
                      "confidence": D.obs["velocity_confidence"].values,
                      "length": D.obs["velocity_length"].values})
        .groupby("annotation", observed=True).agg(["median", "count"]))
summ.to_csv("velocity_by_cluster.csv")

D.write("k409_velocity_dynamical.h5ad")
adata.write("k409_velocity.h5ad")
out = {"velocity_genes_stochastic": int(adata.var["velocity_genes"].sum()),
       "velocity_genes_dynamical": int(D.var["velocity_genes"].sum()),
       "genes_fitted": int(D.var["fit_likelihood"].notna().sum()),
       "median_confidence_stochastic": float(np.median(adata.obs["velocity_confidence_stochastic"])),
       "median_confidence_dynamical": float(np.median(D.obs["velocity_confidence"])),
       "top_drivers": drivers["gene_name"].head(15).tolist()}
json.dump(out, open("velocity_summary.json", "w"), indent=1)
print(json.dumps(out, indent=1))
print(summ.round(3).to_string())
