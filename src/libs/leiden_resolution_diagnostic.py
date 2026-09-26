"""Leiden resolution diagnostics for a single Visium section.

Builds the two DataFrames used in the resolution figure:
  coh  -- spatial coherence, observed vs label-permutation null
  rdf  -- one row per resolution: k, clusters with marker support,
          subsample ARI, smallest cluster, spatial coherence

Input `A` is the processed AnnData (1067 x 17673) with:
  A.obsm["X_pca"]        50 comps on 2000 scaled HVGs
  A.raw                  log1p-normalised counts (use_raw=True for markers)
  A.obs["domain_short"]  the six res-0.5 domain labels
  A.obsm["spatial"]      Visium array coordinates
"""
from sklearn.metrics import adjusted_rand_score
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

RES = [0.3, 0.5, 0.8, 1.0, 1.4]
rng = np.random.default_rng(0)

def spatial_coherence(lab, ii, jj):
    lab = np.asarray(lab)
    return float((lab[ii] == lab[jj]).mean())

# ---- 3. marker support: clusters with >=3 genes p_adj<1e-3 & logFC>1 ---
def n_supported(res, r):
    tmp = r.copy()
    tmp.obs["L"] = tmp.obs[f"L{res}"]
    sc.tl.rank_genes_groups(tmp, "L", method="wilcoxon", use_raw=True)
    n = 0
    for g in tmp.obs["L"].cat.categories:
        d = sc.get.rank_genes_groups_df(tmp, group=g)
        if ((d.pvals_adj < 1e-3) & (d.logfoldchanges > 1)).sum() >= 3:
            n += 1
    return n

# ---- 4. subsample stability: 5 seeds x 80% of spots -------------------
def subsample_ari(res, r, n_seeds=5, frac=0.8):
    ref = r.obs[f"L{res}"].astype(str).values
    out = []
    for s in range(n_seeds):
        idx = np.random.default_rng(s).choice(r.n_obs, int(frac * r.n_obs), replace=False)
        sub = r[np.sort(idx)].copy()
        sc.pp.neighbors(sub, n_neighbors=15, n_pcs=30)
        sc.tl.leiden(sub, resolution=res, key_added="Ls",
                     flavor="igraph", n_iterations=2, directed=False, random_state=0)
        out.append(adjusted_rand_score(ref[np.sort(idx)], sub.obs["Ls"].astype(str).values))
        del sub
    return float(np.mean(out))


def calc_leiden_diagnostic(ad: sc.AnnData):
    # ---- 1. sweep on a fresh neighbour graph -------------------------------
    r = ad.copy()

    sc.pp.neighbors(r, n_neighbors=15, n_pcs=30)
    sc.tl.umap(r, random_state=0)
    for res in RES:
        sc.tl.leiden(r, resolution=res, key_added=f"L{res}",
                    flavor="igraph", n_iterations=2, directed=False, random_state=0)

    # ---- 2. spatial coherence vs permutation null --------------------------
    sq.gr.spatial_neighbors_grid(r, n_neighs=6)
    Wsp = r.obsp["spatial_connectivities"].tocsr()
    ii, jj = Wsp.nonzero()
    keep = ii < jj                      # each hex edge once
    ii, jj = ii[keep], jj[keep]

    rows = []
    for res in RES:
        lab = r.obs[f"L{res}"].values
        obs = spatial_coherence(lab, ii, jj)
        null = np.mean([spatial_coherence(rng.permutation(lab), ii, jj) for _ in range(20)])
        rows.append(dict(res=res, observed=obs, null=null, ratio=obs / null))
    coh = pd.DataFrame(rows)


    rdf = pd.DataFrame([
        dict(res=res,
            k=r.obs[f"L{res}"].nunique(),
            clusters_with_markers=n_supported(res, r),
            subsample_ARI=round(subsample_ari(res, r), 3),
            min_size=int(r.obs[f"L{res}"].value_counts().min()),
            spatial_coherence=round(coh.loc[coh.res == res, "observed"].item(), 3))
        for res in RES
    ])

    # ---- 5. what the extra clusters at res 0.8 are ------------------------
    ct = pd.crosstab(r.obs["L0.8"], r.obs["domain_short"])
    ct["size"] = ct.sum(axis=1)
    ct["purity"] = (ct.iloc[:, :-1].max(axis=1) / ct["size"]).round(2)

    return r, rdf, ct, coh
