"""Pseudotime and fate figures for a fitted scFates object.

Why not scf.pl.trajectory / scf.pl.trends directly: both draw the principal
graph over whatever `basis` you pass, so a graph fitted in harmony/PCA space
appears as a zigzag over a UMAP, and their gene labels collide past ~15
features. This module keeps the embedding panels label-light and draws the
graph only when asked for the space it was actually fitted in.

"Fate" in scFates is not a probability: a cell's fate is the terminal branch
of the principal graph it sits on (`obs['seg']`, `obs['milestones']`). For
probabilistic fate assignment use CellRank, which scFates can import via
`scf.tl.cellrank_to_tree`.

    from libs.plot_fates import plot_pseudotime_fates
    fig = plot_pseudotime_fates(IN, "results/figures/pseudotime_fates.png",
                                label_key="cell.types", basis="X_umap_subset")
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def _panel_letter(ax, letter, size=9):
    """Bold panel letter above the axes, as a separate artist.

    Must not use ax.set_title: these panels already carry a title (the gene
    name), and set_title would replace it.
    """
    ax.text(-0.02, 1.06, letter, transform=ax.transAxes, fontsize=size,
            fontweight="bold", va="bottom", ha="right", clip_on=False)


def terminal_fates(adata):
    """Map each cell to the terminal branch it belongs to.

    Returns (labels, info) where labels is a Series over obs_names naming the
    fate ("fate: <dominant label>" for terminal segments, "progenitor//internal"
    otherwise), and info lists per-segment size and whether it is terminal.
    """
    g = adata.uns["graph"]
    tips = {int(t) for t in np.ravel(g["tips"])}
    seg = adata.obs["seg"].astype(str)
    rows, naming = [], {}
    for _, r in g["pp_seg"].iterrows():
        s = str(int(r["n"]))
        is_term = int(r["to"]) in tips
        members = adata.obs.loc[seg == s]
        dom = (members["cell.types"].astype(str).value_counts().idxmax()
               if len(members) and "cell.types" in members else "n/a")
        naming[s] = f"fate: {dom}" if is_term else "transit"
        rows.append(dict(seg=s, n=int((seg == s).sum()), terminal=is_term,
                         dominant=dom, length=float(r["d"])))
    return seg.map(naming).rename("fate"), pd.DataFrame(rows)


def plot_pseudotime_fates(adata, out_png, label_key="cell.types",
                          basis="X_umap_subset", t_key="t", draw_graph=False,
                          title=None, sizes=(8, 7, 6), figsize=(7.2, 4.6)):
    """Four panels: cell types, pseudotime, fate/branch, pseudotime per type.

    Two independent knobs, easy to confuse:

    figsize : (width, height) in INCHES -- the canvas, passed to plt.figure
    sizes   : (base, secondary, tick) in POINTS -- the font ladder
              (titles/axis-labels, legend/annotation, tick labels)

    Enlarging figsize alone makes text relatively smaller, so scale the
    ladder with it: at figsize=(12, 8) use sizes=(13, 12, 10). Marker areas
    are scaled automatically from the width.

    `draw_graph=True` overlays the principal points -- only meaningful when
    `basis` is the representation the graph was fitted in (uns['ppt']/['epg']
    coordinates live in that space, not in a UMAP).
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    k = figsize[0] / 7.2          # scale markers off the reference width
    pt = 3 * k ** 1.6
    fs_base, fs_small, fs_tick = sizes

    try:  # figure-style skill, when loaded
        apply_figure_style(sizes=sizes)  # noqa: F821
        frame = set_frame  # noqa: F821
        letter = panel_letter  # noqa: F821
    except NameError:
        mpl.rcParams.update({"font.size": fs_small, "axes.titlesize": fs_base,
                             "axes.labelsize": fs_base, "xtick.labelsize": fs_tick,
                             "ytick.labelsize": fs_tick, "figure.dpi": 110})
        frame = lambda ax: [ax.spines[s].set_visible(False) for s in ("top", "right")]  # noqa: E731
        letter = lambda ax, l: _panel_letter(ax, l, size=fs_base + 1)  # noqa: E731, E741

    emb = np.asarray(adata.obsm[basis])[:, :2]
    lab = adata.obs[label_key].astype(str)
    t = adata.obs[t_key].values
    fate, seg_info = terminal_fates(adata)

    cats = sorted(lab.unique())
    pal = dict(zip(cats, plt.get_cmap("tab20")(np.linspace(0, 1, max(len(cats), 2)))))
    fates = sorted(fate.unique())
    fpal = dict(zip(fates, plt.get_cmap("Set2")(np.linspace(0, 1, max(len(fates), 2)))))

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1.0], hspace=0.45, wspace=0.12)

    def blank(ax):
        ax.set_xticks([]); ax.set_yticks([]); ax.margins(0.05)
        for s in ax.spines.values():
            s.set_visible(False)

    ax0 = fig.add_subplot(gs[0, 0])
    for c in cats:
        m = (lab == c).values
        ax0.scatter(emb[m, 0], emb[m, 1], s=pt, c=[pal[c]], lw=0, rasterized=True, label=c)
    ax0.set_title("cell type", fontsize=fs_base, loc="left"); blank(ax0); letter(ax0, "a")

    ax1 = fig.add_subplot(gs[0, 1])
    sc1 = ax1.scatter(emb[:, 0], emb[:, 1], s=pt, c=t, cmap="viridis", lw=0, rasterized=True)
    if draw_graph and "graph" in adata.uns:
        F = np.asarray(adata.uns["graph"]["F"])[:2].T
        B = np.asarray(adata.uns["graph"]["B"])
        for i, j in zip(*np.where(np.triu(B) > 0)):
            ax1.plot(F[[i, j], 0], F[[i, j], 1], c="k", lw=0.6, zorder=3)
    cax = ax1.inset_axes([0.62, 0.04, 0.34, 0.035])
    cb = fig.colorbar(sc1, cax=cax, orientation="horizontal")
    cb.set_label("pseudotime", fontsize=fs_tick, labelpad=1); cb.ax.tick_params(labelsize=fs_tick, pad=1)
    ax1.set_title("pseudotime", fontsize=fs_base, loc="left"); blank(ax1); letter(ax1, "b")

    ax2 = fig.add_subplot(gs[0, 2])
    for f in fates:
        m = (fate == f).values
        ax2.scatter(emb[m, 0], emb[m, 1], s=pt, c=[fpal[f]], lw=0, rasterized=True, label=f)
    ax2.set_title(f"branch ({int(seg_info.terminal.sum())} terminal)", fontsize=fs_base, loc="left")
    blank(ax2); letter(ax2, "c")
    ax2.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False,
               fontsize=fs_tick, markerscale=2.2, handletextpad=0.3, labelspacing=0.25)

    ax3 = fig.add_subplot(gs[1, :])
    order = adata.obs.groupby(label_key, observed=True)[t_key].mean().sort_values().index
    rng = np.random.default_rng(0)
    for i, c in enumerate(order):
        v = t[(lab == str(c)).values]
        ax3.scatter(v, i + rng.uniform(-0.22, 0.22, v.size), s=pt,
                    c=[pal[str(c)]], lw=0, rasterized=True)
        ax3.plot([np.median(v)] * 2, [i - 0.34, i + 0.34], c="0.15", lw=1.2, zorder=4)
    ax3.set_yticks(range(len(order))); ax3.set_yticklabels([str(c) for c in order], fontsize=fs_tick)
    ax3.set_xlabel("pseudotime (bar = median)", fontsize=fs_base)
    ax3.set_title("ordering of cell types along the trajectory", fontsize=fs_base, loc="left")
    ax3.margins(x=0.02, y=0.02); frame(ax3); letter(ax3, "d")

    fig.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=3.5 * k, color=pal[c], label=c)
                        for c in cats],
               loc="upper left", bbox_to_anchor=(0.005, 0.055), ncol=min(len(cats), 5),
               frameon=False, fontsize=fs_tick, handletextpad=0.25, columnspacing=0.9)
    if title:
        fig.suptitle(title, fontsize=round(fs_base * 1.15, 1), x=0.005, ha="left")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    return fig, seg_info


def plot_markers(adata, genes, out_png, basis="X_umap_subset", label_key="cell.types",
                 t_key="t", ncol=3, figsize=(7.2, 5.0), sizes=(8, 7, 6), title=None):
    """Embedding by cell type and pseudotime, plus a grid of marker genes.

    This is the layout of innervation_trajectory.png. Expression is read from
    `.raw` when present, so genes dropped by HVG selection still plot.
    """
    import matplotlib.pyplot as plt

    fs_base, fs_small, fs_tick = sizes
    k = figsize[0] / 7.2
    pt = 3 * k ** 1.6
    try:
        apply_figure_style(sizes=sizes)  # noqa: F821
        letter = panel_letter  # noqa: F821
    except NameError:
        letter = lambda ax, l: _panel_letter(ax, l, size=fs_base + 1)  # noqa: E731, E741

    src = adata.raw if adata.raw is not None else adata
    genes = [g for g in genes if g in src.var_names]
    emb = np.asarray(adata.obsm[basis])[:, :2]
    lab = adata.obs[label_key].astype(str)
    cats = sorted(lab.unique())
    pal = dict(zip(cats, plt.get_cmap("tab20")(np.linspace(0, 1, max(len(cats), 2)))))

    nrow = int(np.ceil(len(genes) / ncol))
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1 + nrow, ncol, height_ratios=[1.6] + [1.0] * nrow,
                          hspace=0.30, wspace=0.10)

    def blank(ax):
        ax.set_xticks([]); ax.set_yticks([]); ax.margins(0.05)
        for s in ax.spines.values():
            s.set_visible(False)

    ax0 = fig.add_subplot(gs[0, :max(ncol // 2, 1)])
    for c in cats:
        m = (lab == c).values
        ax0.scatter(emb[m, 0], emb[m, 1], s=pt, c=[pal[c]], lw=0, rasterized=True, label=c)
    ax0.set_title("cell type", fontsize=fs_base, loc="left"); blank(ax0); letter(ax0, "a")

    ax1 = fig.add_subplot(gs[0, max(ncol // 2, 1):])
    s1 = ax1.scatter(emb[:, 0], emb[:, 1], s=pt, c=adata.obs[t_key].values,
                     cmap="viridis", lw=0, rasterized=True)
    cb = fig.colorbar(s1, ax=ax1, fraction=0.04, pad=0.02)
    cb.set_label("pseudotime", fontsize=fs_tick); cb.ax.tick_params(labelsize=fs_tick)
    ax1.set_title("pseudotime", fontsize=fs_base, loc="left"); blank(ax1); letter(ax1, "b")

    for i, g in enumerate(genes):
        ax = fig.add_subplot(gs[1 + i // ncol, i % ncol])
        v = src[:, g].X
        v = np.asarray(v.todense()).ravel() if hasattr(v, "todense") else np.ravel(v)
        o = np.argsort(v)                      # draw expressing cells on top
        ax.scatter(emb[o, 0], emb[o, 1], s=pt * 0.8, c=v[o], cmap="Greys", lw=0, rasterized=True)
        ax.set_title(g, fontsize=fs_small, loc="left", style="italic")
        blank(ax)
        if i == 0:
            letter(ax, "c")
    fig.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=3.5 * k, color=pal[c], label=c)
                        for c in cats],
               loc="upper left", bbox_to_anchor=(0.005, 0.045), ncol=min(len(cats), 5),
               frameon=False, fontsize=fs_tick, handletextpad=0.25, columnspacing=0.9)
    fig.text(0.98, 0.005, "marker expression, log-normalised (dark = high)",
             ha="right", fontsize=fs_tick, color="0.35")
    if title:
        fig.suptitle(title, fontsize=round(fs_base * 1.15, 1), x=0.005, ha="left")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    return fig


def plot_uncertainty(adata, out_png, basis="X_umap_subset", label_key="cell.types",
                     t_key="t", sd_key="t_sd", figsize=(7.4, 2.6), sizes=(8, 7, 6)):
    """Three panels: t_sd on the embedding, t_sd vs pseudotime, t_sd per cell type.

    This is the layout of pseudotime_uncertainty.png. t_sd only exists after
    scf.tl.pseudotime(n_map>1); with n_map=1 every value is 0.
    """
    import matplotlib.pyplot as plt

    fs_base, fs_small, fs_tick = sizes
    k = figsize[0] / 7.4
    pt = 3 * k ** 1.6
    try:
        apply_figure_style(sizes=sizes)  # noqa: F821
        frame = set_frame  # noqa: F821
        letter = panel_letter  # noqa: F821
    except NameError:
        frame = lambda ax: [ax.spines[s].set_visible(False) for s in ("top", "right")]  # noqa: E731
        letter = lambda ax, l: _panel_letter(ax, l, size=fs_base + 1)  # noqa: E731, E741

    emb = np.asarray(adata.obsm[basis])[:, :2]
    tsd = adata.obs[sd_key].values
    t = adata.obs[t_key].values
    rng_t = float(t.max() - t.min())
    lab = adata.obs[label_key].astype(str)

    fig, axes = plt.subplots(1, 3, figsize=figsize,
                             gridspec_kw=dict(wspace=0.45, width_ratios=[1, 1, 1.15]))
    ax = axes[0]
    s = ax.scatter(emb[:, 0], emb[:, 1], s=pt, c=tsd, cmap="magma_r", lw=0, rasterized=True)
    ax.set_title("projection uncertainty\nacross the embedding", fontsize=fs_base, loc="left")
    ax.set_xticks([]); ax.set_yticks([]); ax.margins(0.06)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cax = ax.inset_axes([0.0, -0.10, 0.5, 0.05])
    cb = fig.colorbar(s, cax=cax, orientation="horizontal")
    cb.set_label("t s.d.", fontsize=fs_tick, labelpad=1); cb.ax.tick_params(labelsize=fs_tick, pad=1)
    letter(ax, "a")

    ax = axes[1]
    ax.scatter(t, tsd, s=pt, c="0.35", lw=0, rasterized=True)
    ax.axhline(0.05 * rng_t, color="crimson", lw=0.8, ls="--")
    ax.text(t.min(), 0.05 * rng_t, " 5% of range", va="bottom", fontsize=fs_tick, color="crimson")
    ax.set_xlabel("pseudotime", fontsize=fs_base); ax.set_ylabel("t s.d.", fontsize=fs_base)
    ax.set_title("uncertainty vs position", fontsize=fs_base, loc="left")
    ax.margins(0.04); frame(ax); letter(ax, "b")

    ax = axes[2]
    order = adata.obs.groupby(label_key, observed=True)[sd_key].median().sort_values().index
    rng = np.random.default_rng(0)
    for i, c in enumerate(order):
        v = tsd[(lab == str(c)).values]
        ax.scatter(v, i + rng.uniform(-0.2, 0.2, v.size), s=pt * 0.8, c="0.45", lw=0, rasterized=True)
        ax.plot([np.median(v)] * 2, [i - 0.34, i + 0.34], c="crimson", lw=1.2, zorder=4)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([str(c) for c in order], fontsize=fs_tick)
    ax.set_xlabel("t s.d. (bar = median)", fontsize=fs_base)
    ax.set_title("per cell type", fontsize=fs_base, loc="left")
    ax.margins(x=0.06, y=0.03); frame(ax); letter(ax, "c")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    return fig


def progenitor_tip(adata, prefixes=("SCP",), label_key="cell.types", n_cells=30):
    """Tip node whose nearest cells are richest in progenitor labels -- for CURVES.

    Returns (tip, {tip: fraction}). Use progenitor_node instead when the
    progenitor pool sits in the interior of a tree (e.g. EPDCs in fibroblasts),
    where every tip is a differentiated state and tip-rooting would order the
    lineage backwards.
    """
    R = np.asarray(adata.obsm["X_R"])
    labs = adata.obs[label_key].astype(str).values
    frac = {}
    for tip in map(int, np.ravel(adata.uns["graph"]["tips"])):
        near = labs[np.argsort(R[:, tip])[::-1][:n_cells]]
        frac[tip] = float(np.mean([l.startswith(tuple(prefixes)) for l in near]))
    return max(frac, key=frac.get), frac


def progenitor_node(adata, labels, label_key="cell.types"):
    """Node (tip OR internal) carrying the most assignment weight from `labels`.

    Correct rooting for a tree whose progenitors are interior. Returns
    (node, weight_vector).
    """
    R = np.asarray(adata.obsm["X_R"])
    m = adata.obs[label_key].astype(str).isin(set(labels)).values
    w = R[m].mean(0)
    return int(np.argmax(w)), w
