"""Neighbours, Leiden at three resolutions, UMAP and marker-based annotation (step 4).

Labels come from MEAN MARKER EXPRESSION per cluster, not module-score argmax.
Module scores are centred against a background gene set, which in FACS-sorted
lymphocytes mislabels naive CD8 clusters as CD4 (CD4 mRNA is barely captured:
<=0.42 mean in every cluster of this dataset, while IL7R is high in naive CD8).
Every label is written to clustering.json together with the marker table it came
from (cluster_marker_means.csv) so the call can be audited or overridden.
"""
import json
import numpy as np, pandas as pd, scanpy as sc

KEY = "leiden_0.8"
MARKER_PANEL = ["PTPRC", "CD3D", "CD3E", "CD2", "CD8A", "CD8B", "CD4", "IL7R", "SELL", "CCR7",
                "TCF7", "LEF1", "GZMK", "GZMB", "GNLY", "NKG7", "PRF1", "KLRD1", "NCAM1",
                "FCGR3A", "KLRF1", "TYROBP", "FCER1G", "FOXP3", "IL2RA", "CTLA4", "PDCD1",
                "LAG3", "HAVCR2", "LAYN", "TOX", "ITGAE", "CXCR6", "MKI67", "TOP2A", "MS4A1",
                "CD79A", "LYZ", "CD14", "C1QA", "MLANA", "PMEL", "TRDC", "TRGC1"]


def label_cluster(p):
    """p: Series of mean log-normalised expression for one cluster."""
    g = lambda k: float(p.get(k, 0.0))
    if g("MS4A1") > 0.5 or g("CD79A") > 0.5:
        return "B cell"
    if g("LYZ") > 0.8 or g("C1QA") > 0.5 or g("CD14") > 0.5:
        return "myeloid"
    if g("MLANA") > 0.5 or g("PMEL") > 0.5:
        return "malignant"
    cd8 = max(g("CD8A"), g("CD8B"))
    if g("FOXP3") > 0.4 or (g("IL2RA") > 0.6 and cd8 < 0.5):
        return "CD4 Treg"
    # CD4 mRNA is not usable for lineage calling in 10x data, so CD8A/CD8B decide;
    # the middle band is reported as unresolved rather than guessed.
    lineage = "CD8" if cd8 >= 0.9 else ("CD4" if cd8 < 0.4 else "CD4/CD8")
    if g("MKI67") > 0.5 or g("TOP2A") > 0.4:
        state = "cycling"
    elif g("PDCD1") > 1.0 or g("LAG3") > 0.5 or g("HAVCR2") > 0.4 or g("LAYN") > 0.3:
        state = "exhausted"
    elif g("SELL") > 0.8 or g("CCR7") > 0.6 or g("TCF7") > 0.6:
        state = "naive/CM"
    elif g("GNLY") > 1.5 or g("PRF1") > 1.0:
        state = "effector GNLY+"
    elif g("GZMK") > 0.6:
        state = "GZMK+ EM"
    elif g("IL7R") > 1.2:
        state = "memory"
    else:
        state = "effector"
    return f"{lineage} {state}"


adata = sc.read_h5ad("k409_pca.h5ad")
rep = adata.uns["representation"]
sc.pp.neighbors(adata, n_neighbors=15, use_rep=rep, n_pcs=30)
for res in (0.4, 0.8, 1.2):
    sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}", flavor="igraph", n_iterations=2)
sc.tl.umap(adata, min_dist=0.3)

sc.tl.rank_genes_groups(adata, KEY, method="wilcoxon", n_genes=50)
pd.DataFrame({c: adata.uns["rank_genes_groups"]["names"][c][:15]
              for c in adata.obs[KEY].cat.categories}).to_csv("cluster_top_genes.csv", index=False)

panel = [g for g in MARKER_PANEL if g in adata.var_names]
prof = sc.get.obs_df(adata, keys=panel + [KEY]).groupby(KEY, observed=True).mean()
prof.insert(0, "n_cells", adata.obs[KEY].value_counts().sort_index().values)
prof.round(3).to_csv("cluster_marker_means.csv")

labels = {cl: label_cluster(prof.loc[cl]) for cl in prof.index}
# disambiguate clusters that land on the same label, by dominant tissue then size
seen = {}
for cl in prof.index:
    lab = labels[cl]
    if list(labels.values()).count(lab) > 1:
        tis = adata.obs.loc[adata.obs[KEY] == cl, "tissue"].value_counts().idxmax()
        suffix = "blood" if "blood" in str(tis).lower() else "tissue"
        labels[cl] = f"{lab} ({suffix})"
        n = seen.get(labels[cl], 0) + 1
        seen[labels[cl]] = n
        if n > 1:
            labels[cl] = f"{labels[cl]} {n}"
adata.obs["annotation"] = adata.obs[KEY].map(labels).astype("category")
adata.obs["compartment"] = np.where(adata.obs["tissue"].astype(str).str.contains("blood"),
                                    "blood", "tissue (tumour + node)")

comp = pd.crosstab(adata.obs["annotation"], adata.obs["tissue"])
comp.to_csv("cluster_composition.csv")
adata.write("k409_clustered.h5ad")
json.dump({"representation": rep, "n_pcs": 30, "resolution_used": KEY,
           "n_clusters": {f"leiden_{r}": int(adata.obs[f'leiden_{r}'].nunique()) for r in (0.4, 0.8, 1.2)},
           "labels": labels,
           "labelling_note": "mean marker expression per cluster (cluster_marker_means.csv); "
                             "CD4 identity inferred from FOXP3/IL2RA and CD8A/CD8B-low, because "
                             "CD4 mRNA is poorly captured in 10x data"},
          open("clustering.json", "w"), indent=1)
print(comp.to_string())
print(json.dumps(labels, indent=1))
