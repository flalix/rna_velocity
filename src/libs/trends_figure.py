"""Publication-grade trends heatmap for a scFates fit.

`scf.pl.trends` is unusable as a deliverable for two reasons: its gene
labels collide into illegibility past ~15 features, and it draws the
principal graph over whatever `basis` you pass -- so a tree fitted in PCA
space appears as a zigzag over a UMAP. This builds the same information
directly from `layers['fitted']`:

  panel a  genes x pseudotime-bin heatmap, z-scored per gene,
           ordered by the pseudotime at which each gene peaks
  panel b  fitted trends for a named set of genes

    from trends_figure import trends_figure
    trends_figure(adata, "Innervation", "trends_innervation.png",
                  genes=["SOX10", "PLP1", "ASCL1", "CHGB", "STMN2", "TUBB3"])
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans


def gene_modules(adata, n_bins=60, k=6, seed=0):
    """Z-scored fitted profiles for every significant gene, clustered into modules.

    Binning is per segment, not pooled over pseudotime: on a branching tree,
    cells on different branches share t values, so pooling averages distinct
    fates together and produces spurious spikes. Columns are therefore
    (segment, bin-within-segment), segments ordered by mean pseudotime.

    Returns (z, t_mid, modules, seg_of_col): z is genes x column, rows ordered
    by module then within-module peak; modules renumbered so module 1 peaks
    earliest.
    """
    fit = pd.DataFrame(np.asarray(adata.layers["fitted"]),
                       index=adata.obs_names, columns=adata.var_names)
    seg = adata.obs["seg"].astype(str)
    seg_order = adata.obs.groupby(seg, observed=True)["t"].mean().sort_values().index
    per_seg = max(int(round(n_bins / max(len(seg_order), 1))), 4)

    cols, t_parts, seg_parts = [], [], []
    for s in seg_order:
        m = (seg == s).values
        if m.sum() < per_seg:                     # too few cells to bin: one column
            cols.append(fit[m].mean().rename((s, 0)))
            t_parts.append(adata.obs["t"][m].mean()); seg_parts.append(s)
            continue
        b = pd.qcut(adata.obs["t"][m], q=per_seg, labels=False, duplicates="drop")
        g = fit[m].groupby(b, observed=True).mean()
        tm = adata.obs["t"][m].groupby(b, observed=True).mean()
        for bi in g.index:
            cols.append(g.loc[bi].rename((s, bi)))
            t_parts.append(tm.loc[bi]); seg_parts.append(s)

    binned = pd.DataFrame(cols)
    t_mid = pd.Series(t_parts, index=binned.index)
    seg_of_col = pd.Series(seg_parts, index=binned.index)
    z = ((binned - binned.mean()) / binned.std(ddof=0)).T.dropna()   # genes x columns

    km = KMeans(n_clusters=k, n_init=10, random_state=seed).fit(z.values)
    mod = pd.Series(km.labels_, index=z.index)
    peak_t = {m: float(t_mid.iloc[int(np.argmax(z.loc[mod == m].mean(0).values))])
              for m in sorted(mod.unique())}
    remap = {m: i + 1 for i, m in enumerate(sorted(peak_t, key=peak_t.get))}
    mod = mod.map(remap)
    within = z.apply(lambda r: np.argmax(r.values), axis=1)
    order = pd.DataFrame({"mod": mod, "peak": within}).sort_values(["mod", "peak"]).index
    return z.loc[order], t_mid, mod.loc[order], seg_of_col


def trends_figure(adata, title, out_png, n_bins=60, k=6, n_label=3, style_fn=None, figsize=(7.4, 4.2)):
    z, t_mid, mod, seg_of_col = gene_modules(adata, n_bins=n_bins, k=k)
    A = adata.var["A"]
    # column positions where the segment changes
    seg_changes = [i for i in range(1, len(seg_of_col))
                   if seg_of_col.iloc[i] != seg_of_col.iloc[i - 1]]
    seg_spans = {}
    start = 0
    for i in seg_changes + [len(seg_of_col)]:
        seg_spans[seg_of_col.iloc[start]] = (start, i)
        start = i

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.04, 1], width_ratios=[1.35, 1],
                          hspace=0.05, wspace=0.42)

    axt = fig.add_subplot(gs[0, 0])
    axt.imshow(t_mid.values[None, :], aspect="auto", cmap="viridis")
    axt.set_xticks([]); axt.set_yticks([])
    axt.set_title(f"{title} — {z.shape[0]:,} genes associated with pseudotime, "
                  f"{k} modules", fontsize=8, loc="left", pad=4)

    axh = fig.add_subplot(gs[1, 0])
    im = axh.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2, rasterized=True)
    bounds = np.cumsum(mod.value_counts().sort_index().values)
    for b in bounds[:-1]:
        axh.axhline(b - 0.5, color="k", lw=0.6)
    centres = np.concatenate([[0], bounds[:-1]]) + mod.value_counts().sort_index().values / 2
    axh.set_yticks(centres)
    axh.set_yticklabels([f"M{m} ({n})" for m, n in mod.value_counts().sort_index().items()],
                        fontsize=6)
    for c in seg_changes:
        axh.axvline(c - 0.5, color="k", lw=0.5)
        axt.axvline(c - 0.5, color="k", lw=0.5)
    if len(seg_spans) > 1:
        axh.set_xticks([np.mean(v) for v in seg_spans.values()])
        axh.set_xticklabels([f"seg {s}" for s in seg_spans], fontsize=5.5)
        axh.set_xlabel("branch, each binned by pseudotime (colour strip)", fontsize=7)
    else:
        axh.set_xticks([0, z.shape[1] - 1])
        axh.set_xticklabels([f"{t_mid.iloc[0]:.0f}", f"{t_mid.iloc[-1]:.0f}"], fontsize=6)
        axh.set_xlabel("pseudotime", fontsize=7)
    cax = axh.inset_axes([0.0, -0.13, 0.35, 0.035])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("fitted expression (z per gene)", fontsize=6, labelpad=1)
    cb.ax.tick_params(labelsize=6, pad=1)

    # name the strongest genes per module alongside the heatmap
    for m in sorted(mod.unique()):
        idx = np.where(mod.values == m)[0]
        names = A.loc[mod.index[idx]].sort_values(ascending=False).index[:n_label]
        axh.text(z.shape[1] * 1.02, idx.mean(), "\n".join(names), fontsize=5.2,
                 style="italic", va="center", ha="left")

    axl = fig.add_subplot(gs[1, 1])
    cmap = plt.get_cmap("viridis")(np.linspace(0.05, 0.95, k))
    x = np.arange(z.shape[1])
    for i, m in enumerate(sorted(mod.unique())):
        axl.plot(x, z.loc[mod[mod == m].index].mean(0).values,
                 lw=1.5, color=cmap[i], label=f"M{m}")
    axl.axhline(0, color="0.8", lw=0.6, zorder=0)
    for c in seg_changes:
        axl.axvline(c - 0.5, color="0.75", lw=0.5, zorder=0)
    if len(seg_spans) > 1:
        axl.set_xticks([np.mean(v) for v in seg_spans.values()])
        axl.set_xticklabels([f"{s}" for s in seg_spans], fontsize=5.5)
        axl.set_xlabel("branch (segment), binned by pseudotime", fontsize=7)
    else:
        axl.set_xticks([0, z.shape[1] - 1])
        axl.set_xticklabels([f"{t_mid.iloc[0]:.0f}", f"{t_mid.iloc[-1]:.0f}"], fontsize=6)
        axl.set_xlabel("pseudotime", fontsize=7)
    axl.set_ylabel("module mean (z)", fontsize=7)
    axl.set_title("module profiles", fontsize=8, loc="left")
    axl.legend(frameon=False, fontsize=6, labelspacing=0.25, handlelength=1.2, ncol=2)
    if style_fn is not None:
        style_fn(axl)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    return fig, mod
