"""
abriplot.cli
============
Command-line interface entry point.
Registered as the `abriplot` console script in pyproject.toml.
"""

import sys
import argparse

from .io      import load_all
from .newick  import tip_order
from .plot    import (build_matrix, order_by_tree,
                      sort_genes_by_prevalence, make_figure,
                      resistance_classes)


def build_parser():
    ap = argparse.ArgumentParser(
        prog="abriplot",
        description=(
            "abriplot — publication-quality AMR + Virulence heatmap\n"
            "with phylogenetic tree from abricate TSV output.\n\n"
            "Combines CARD (*__card.tsv) and VFDB (*__vfdb.tsv) results\n"
            "into a single figure ordered by a Newick phylogeny."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples\n"
            "--------\n"
            "  # Basic run\n"
            "  abriplot --input_dir results/ --tree tree.nwk --output fig.pdf\n\n"
            "  # Filter to top genes, set custom title\n"
            "  abriplot --input_dir results/ --tree tree.nwk --output fig.pdf \\\n"
            "           --top_amr 40 --top_vf 30 \\\n"
            "           --title 'AMR genes in V. cholerae 2002-2019'\n"
        ),
    )

    # ── Required ───────────────────────────────────────────────────────────
    req = ap.add_argument_group("required arguments")
    req.add_argument(
        "--input_dir", required=True, metavar="DIR",
        help="Directory containing *__card.tsv and *__vfdb.tsv abricate output files",
    )
    req.add_argument(
        "--tree", required=True, metavar="FILE",
        help="Newick tree file (.nwk / .tree)",
    )

    # ── Output ─────────────────────────────────────────────────────────────
    out = ap.add_argument_group("output")
    out.add_argument(
        "--output", default="abriplot_figure.pdf", metavar="FILE",
        help="Output path — pdf, png, or svg [default: abriplot_figure.pdf]. "
             "A PNG copy is always saved alongside a PDF/SVG.",
    )
    out.add_argument(
        "--title",
        default="Phylogenetic Distribution of AMR and Virulence Genes",
        metavar="TEXT",
        help="Figure title [default: 'Phylogenetic Distribution of AMR and Virulence Genes']",
    )

    # ── Filtering ──────────────────────────────────────────────────────────
    flt = ap.add_argument_group("filtering")
    flt.add_argument(
        "--min_identity", type=float, default=80.0, metavar="FLOAT",
        help="Minimum %% identity threshold [default: 80]",
    )
    flt.add_argument(
        "--min_coverage", type=float, default=80.0, metavar="FLOAT",
        help="Minimum %% coverage threshold [default: 80]",
    )
    flt.add_argument(
        "--top_amr", type=int, default=None, metavar="INT",
        help="Keep only the N most prevalent AMR genes [default: all]",
    )
    flt.add_argument(
        "--top_vf", type=int, default=None, metavar="INT",
        help="Keep only the N most prevalent virulence genes [default: all]",
    )

    return ap


def main():
    ap = build_parser()
    args = ap.parse_args()

    print(f"\n{'='*60}")
    print("  abriplot")
    print(f"{'='*60}")
    for k, v in vars(args).items():
        print(f"  {k:<16}: {v}")
    print(f"{'='*60}\n")

    # ── Parse tree ─────────────────────────────────────────────────────────
    with open(args.tree) as fh:
        nwk = fh.read().strip()
    ordered_tips, root, leaves = tip_order(nwk)
    print(f"  Tree tips  : {len(ordered_tips)}")

    # ── Load data ──────────────────────────────────────────────────────────
    amr_records, vf_records = load_all(
        args.input_dir, args.min_identity, args.min_coverage
    )
    if not amr_records and not vf_records:
        sys.exit("ERROR: No records passed the identity/coverage filters. "
                 "Check --min_identity and --min_coverage.")

    # ── Build matrices ─────────────────────────────────────────────────────
    amr_mat, amr_meta = build_matrix(amr_records)
    vf_mat,  _        = build_matrix(vf_records)

    amr_mat = order_by_tree(amr_mat, ordered_tips)
    vf_mat  = order_by_tree(vf_mat,  ordered_tips)
    vf_mat  = vf_mat.reindex(amr_mat.index, fill_value=0)

    amr_mat = sort_genes_by_prevalence(amr_mat)
    vf_mat  = sort_genes_by_prevalence(vf_mat)

    if args.top_amr:
        amr_mat = amr_mat.iloc[:, :args.top_amr]
    if args.top_vf:
        vf_mat  = vf_mat.iloc[:,  :args.top_vf]

    print(f"\n  Samples  : {amr_mat.shape[0]}")
    print(f"  AMR genes: {amr_mat.shape[1]}")
    print(f"  VF genes : {vf_mat.shape[1]}\n")

    # ── Plot ───────────────────────────────────────────────────────────────
    make_figure(
        amr_mat, vf_mat, amr_meta,
        root, leaves, amr_mat.index.tolist(),
        args.output,
        title=args.title,
    )

    # ── Summary ────────────────────────────────────────────────────────────
    n = amr_mat.shape[0]
    print("\n── Top 15 AMR genes by prevalence ──")
    for g, c in (amr_mat > 0).sum().sort_values(ascending=False).head(15).items():
        cls = ", ".join(sorted(resistance_classes(amr_meta.get(g, ""))))
        print(f"  {g:<30} {c:>3}/{n} ({100*c/n:.0f}%)  [{cls or 'unclassified'}]")
    print()


if __name__ == "__main__":
    main()
