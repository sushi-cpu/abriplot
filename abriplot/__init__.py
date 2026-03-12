"""
abriplot
========
Publication-quality AMR + Virulence gene heatmap with phylogenetic tree
from abricate TSV output.

Quick start
-----------
From the command line::

    abriplot --input_dir results/ --tree tree.nwk --output fig.pdf

As a library::

    from abriplot.io     import load_all
    from abriplot.newick import tip_order
    from abriplot.plot   import build_matrix, order_by_tree, sort_genes_by_prevalence, make_figure

"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("abriplot")
except PackageNotFoundError:
    __version__ = "unknown"

__all__ = ["__version__"]
