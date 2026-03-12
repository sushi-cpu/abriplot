"""
abriplot.plot
=============
Core figure-generation logic.  Called by the CLI entry point (abriplot.cli)
but also importable as a library for use in notebooks or custom pipelines.
"""

import os
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import rcParams

warnings.filterwarnings("ignore")

rcParams.update({
    "font.family":       "DejaVu Sans",
    "axes.linewidth":    0.7,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "pdf.fonttype":      42,
    "svg.fonttype":      "none",
})

AMR_PRESENT  = "#2171b5"
VF_PRESENT   = "#cb181d"
AMR_CMAP = LinearSegmentedColormap.from_list("amr", [AMR_PRESENT, AMR_PRESENT])
VF_CMAP  = LinearSegmentedColormap.from_list("vf",  [VF_PRESENT,  VF_PRESENT])

TREE_COLOR   = "#2d3436"
GRID_COL     = "#ffffff"
BG           = "#f8f9fa"
LABEL_COL    = "#2d3436"
ABSENT_STRIP = "#eeeeee"

CLASS_PALETTE = {}


# ── Drug-class helpers ─────────────────────────────────────────────────────────

def resistance_classes(res_string):
    """Return the set of resistance class strings from a semicolon-separated field."""
    if not isinstance(res_string, str) or not res_string.strip():
        return set()
    return {p.strip().lower() for p in res_string.split(";") if p.strip()}


def build_class_palette(amr_meta):
    """Auto-assign a distinct colour to every drug class present in the data."""
    global CLASS_PALETTE
    import matplotlib.cm as _cm
    import matplotlib.colors as _mc

    all_classes = set()
    for res in amr_meta.values():
        all_classes |= resistance_classes(res)
    classes_sorted = sorted(all_classes)

    tab20  = [_mc.to_hex(_cm.get_cmap("tab20" )(i / 20)) for i in range(20)]
    tab20b = [_mc.to_hex(_cm.get_cmap("tab20b")(i / 20)) for i in range(20)]
    tab20c = [_mc.to_hex(_cm.get_cmap("tab20c")(i / 20)) for i in range(20)]
    colours = tab20 + tab20b + tab20c

    CLASS_PALETTE = {cls: colours[i % len(colours)]
                     for i, cls in enumerate(classes_sorted)}

    print(f"  Drug classes found : {len(classes_sorted)}")
    for c in classes_sorted:
        print(f"    {c}")
    return CLASS_PALETTE


# ── Matrix building ────────────────────────────────────────────────────────────

def build_matrix(records, key="gene"):
    """Presence/absence matrix: 1 if gene detected in sample, 0 otherwise."""
    seen = defaultdict(set)
    meta = {}
    for r in records:
        g, s = r[key], r["sample"]
        seen[s].add(g)
        if g not in meta:
            meta[g] = r.get("resistance", r.get("product", ""))
    all_genes   = sorted(meta.keys())
    all_samples = sorted(seen.keys())
    mat = pd.DataFrame(0, index=all_samples, columns=all_genes, dtype=int)
    for s in all_samples:
        for g in seen[s]:
            if g in mat.columns:
                mat.loc[s, g] = 1
    return mat, meta


def order_by_tree(mat, ordered_tips):
    present = [t for t in ordered_tips if t in mat.index]
    missing = [t for t in ordered_tips if t not in mat.index]
    extra   = [s for s in mat.index    if s not in ordered_tips]
    if missing:
        print(f"  WARNING: {len(missing)} tree tip(s) absent from abricate data: "
              f"{missing[:5]}{'...' if len(missing) > 5 else ''}")
    if extra:
        print(f"  WARNING: {len(extra)} sample(s) in abricate data not in tree: "
              f"{extra[:5]}{'...' if len(extra) > 5 else ''}")
    return mat.reindex(present + extra, fill_value=0)


def sort_genes_by_prevalence(mat):
    return mat[(mat > 0).sum(axis=0).sort_values(ascending=False).index]


# ── Class-strip ────────────────────────────────────────────────────────────────

def build_class_strip(amr_mat, amr_meta):
    palette  = build_class_palette(amr_meta)
    classes  = sorted(palette.keys())
    strip    = pd.DataFrame(False, index=classes, columns=amr_mat.columns)
    for g in amr_mat.columns:
        for cls in resistance_classes(amr_meta.get(g, "")):
            if cls in strip.index:
                strip.loc[cls, g] = True
    return strip, classes, palette


def draw_class_strip(ax, strip_df, classes, palette, **_kwargs):
    nr, nc      = len(classes), strip_df.shape[1]
    absent_rgba = np.array(matplotlib.colors.to_rgba(ABSENT_STRIP))

    img = np.zeros((nr, nc, 4))
    for i, cls in enumerate(classes):
        col_rgba = np.array(matplotlib.colors.to_rgba(palette[cls]))
        for j, gene in enumerate(strip_df.columns):
            img[i, j] = col_rgba if strip_df.loc[cls, gene] else absent_rgba

    ax.imshow(img, aspect="auto", interpolation="nearest")
    ax.set_xticks(np.arange(-0.5, nc, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, nr, 1), minor=True)
    ax.grid(which="minor", color=GRID_COL, linewidth=0.4)
    ax.tick_params(which="minor", length=0)
    ax.set_yticks(range(nr))
    ax.set_yticklabels([c.replace("-", "-\n") for c in classes],
                       fontsize=5.2, color=LABEL_COL, va="center")
    ax.tick_params(axis="y", length=0)
    ax.set_xticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.4); sp.set_edgecolor("#cccccc")

    patches = [mpatches.Patch(facecolor=palette[c],
                              label=c.replace("-", " ").title(), linewidth=0)
               for c in classes]
    ax.legend(handles=patches, title="Drug class", title_fontsize=5.5,
              fontsize=5.0, loc="upper left", bbox_to_anchor=(1.01, 1.02),
              frameon=True, framealpha=0.95, edgecolor="#cccccc",
              handlelength=0.85, handleheight=0.80,
              borderpad=0.5, labelspacing=0.25)
    ax.set_title("  Resistance Drug Classes per Gene", fontsize=7.5,
                 fontweight="bold", color="white",
                 backgroundcolor="#444444", pad=3, loc="left", x=0)


# ── Heatmap block ──────────────────────────────────────────────────────────────

def draw_block(ax, mat, cmap, absent_hex, title, title_bg,
               gene_fontsize=5.5, show_yticks=False):
    nr, nc = mat.shape
    if nc == 0:
        ax.axis("off")
        ax.set_title(title, fontsize=8, fontweight="bold", color="white",
                     backgroundcolor=title_bg, pad=4, loc="left", x=0)
        return

    absent_rgba  = np.array(matplotlib.colors.to_rgba(absent_hex))
    present_rgba = np.array(matplotlib.colors.to_rgba(cmap(0.5)))
    img = np.zeros((nr, nc, 4))
    for j in range(nc):
        for i in range(nr):
            img[i, j] = present_rgba if mat.iloc[i, j] > 0 else absent_rgba

    ax.imshow(img, aspect="auto", interpolation="nearest")
    ax.set_xticks(np.arange(-0.5, nc, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, nr, 1), minor=True)
    ax.grid(which="minor", color=GRID_COL, linewidth=0.35, zorder=3)
    ax.tick_params(which="minor", length=0)
    ax.set_xticks(range(nc))
    ax.set_xticklabels(mat.columns, rotation=90, ha="center", va="top",
                       fontsize=gene_fontsize, fontstyle="italic", color=LABEL_COL)
    ax.xaxis.set_ticks_position("bottom")
    if show_yticks:
        ax.set_yticks(range(nr))
        ax.set_yticklabels(mat.index, fontsize=6, color=LABEL_COL)
    else:
        ax.set_yticks([])
    ax.set_title(title, fontsize=8, fontweight="bold", color="white",
                 backgroundcolor=title_bg, pad=4, loc="left", x=0)
    for sp in ax.spines.values():
        sp.set_linewidth(0.5); sp.set_edgecolor("#bbbbbb")


# ── Main figure ────────────────────────────────────────────────────────────────

def make_figure(amr_mat, vf_mat, amr_meta,
                root, leaves, tip_labels_ordered,
                output_path,
                title="Phylogenetic Distribution of AMR and Virulence Genes"):

    n_samples = amr_mat.shape[0]
    n_amr     = amr_mat.shape[1]
    n_vf      = vf_mat.shape[1]

    strip_df, classes, palette = build_class_strip(amr_mat, amr_meta)
    n_classes = len(classes)

    row_h       = 0.22
    strip_row_h = 0.15
    base_h      = max(10, n_samples * row_h)
    strip_h     = max(1.4, n_classes * strip_row_h)
    tree_w      = 2.8
    amr_w       = max(3.0, n_amr * 0.165)
    vf_w        = max(2.0, n_vf  * 0.165)
    gap_w       = 0.35
    right_pad   = 2.2
    total_w     = tree_w + amr_w + gap_w + vf_w + right_pad + 0.5
    total_h     = base_h + strip_h + 3.2

    fig = plt.figure(figsize=(total_w, total_h), dpi=300)
    fig.patch.set_facecolor(BG)

    fig.text(0.5, 0.998, title, ha="center", va="top",
             fontsize=13, fontweight="bold", color=LABEL_COL)
    fig.text(0.5, 0.983,
             "Presence/absence matrix  ·  AMR: CARD  ·  Virulence: VFDB  ·  "
             "Strip = drug classes conferred per gene",
             ha="center", va="top", fontsize=6.5, color="#636e72")

    h_ratios = [strip_h, base_h]
    w_ratios = [tree_w, amr_w, gap_w * 0.04, vf_w]
    gs = gridspec.GridSpec(2, 4, height_ratios=h_ratios, width_ratios=w_ratios,
                           left=0.01, right=0.88, top=0.930, bottom=0.13,
                           hspace=0.04, wspace=0.015)

    ax_strip_space = fig.add_subplot(gs[0, 0]); ax_strip_space.axis("off")
    ax_strip       = fig.add_subplot(gs[0, 1])
    fig.add_subplot(gs[0, 2]).set_visible(False)
    fig.add_subplot(gs[0, 3]).set_visible(False)
    ax_tree = fig.add_subplot(gs[1, 0])
    ax_amr  = fig.add_subplot(gs[1, 1])
    fig.add_subplot(gs[1, 2]).set_visible(False)
    ax_vf   = fig.add_subplot(gs[1, 3])

    draw_class_strip(ax_strip, strip_df, classes, palette)
    draw_tree(ax_tree, root, leaves, tip_labels_ordered)
    ax_tree.set_title("  Phylogeny", fontsize=8, fontweight="bold",
                      color="white", backgroundcolor=TREE_COLOR,
                      pad=4, loc="left", x=0)
    draw_block(ax_amr, amr_mat, AMR_CMAP, "#e8f4fd",
               "  AMR Genes  (CARD)", "#2171b5")
    draw_block(ax_vf, vf_mat, VF_CMAP, "#fff0ee",
               "  Virulence Genes  (VFDB)", "#cb181d")

    mid_x = (ax_amr.get_position().x1 + ax_vf.get_position().x0) / 2
    fig.add_artist(mlines.Line2D(
        [mid_x, mid_x],
        [ax_amr.get_position().y0, ax_strip.get_position().y1],
        transform=fig.transFigure, color="#aaaaaa", lw=0.8, linestyle="--"))

    pa_patches = [
        mpatches.Patch(facecolor=AMR_PRESENT, label="AMR gene present",       linewidth=0),
        mpatches.Patch(facecolor="#e8f4fd",   label="AMR gene absent",        linewidth=0, edgecolor="#aaaaaa"),
        mpatches.Patch(facecolor=VF_PRESENT,  label="Virulence gene present", linewidth=0),
        mpatches.Patch(facecolor="#fff0ee",   label="Virulence gene absent",  linewidth=0, edgecolor="#aaaaaa"),
    ]
    fig.legend(handles=pa_patches, fontsize=6, loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, 0.005), frameon=True, framealpha=0.92,
               edgecolor="#cccccc", handlelength=1.0, handleheight=0.9)

    plt.savefig(output_path, dpi=300, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    print(f"  Saved -> {output_path}")

    base, ext = os.path.splitext(output_path)
    if ext.lower() != ".png":
        png_path = base + ".png"
        plt.savefig(png_path, dpi=300, bbox_inches="tight",
                    facecolor=fig.get_facecolor())
        print(f"  Saved -> {png_path}")

    plt.close(fig)


# ── Tree drawing (also used by CLI directly) ───────────────────────────────────

def draw_tree(ax, root, leaves, tip_labels, color=TREE_COLOR, lw=0.8):
    max_x = max(lf.x for lf in leaves)

    def _draw(node):
        if node.is_leaf():
            return
        ys = [c.y for c in node.children]
        ax.plot([node.x, node.x], [min(ys), max(ys)],
                color=color, lw=lw, solid_capstyle="butt")
        for c in node.children:
            ax.plot([node.x, c.x], [c.y, c.y],
                    color=color, lw=lw, solid_capstyle="butt")
            _draw(c)
    _draw(root)

    for lf, label in zip(leaves, tip_labels):
        ax.text(max_x * 1.02, lf.y, label,
                va="center", ha="left", fontsize=6.5, color=LABEL_COL)

    ax.set_xlim(-max_x * 0.05, max_x * 1.45)
    ax.set_ylim(-0.8, len(leaves) - 0.2)
    ax.invert_yaxis()
    ax.axis("off")
