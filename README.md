# abriplot

**Publication-quality AMR and Virulence gene heatmap with phylogenetic tree from [abricate](https://github.com/tseemann/abricate) output.**

[![PyPI version](https://badge.fury.io/py/abriplot.svg)](https://pypi.org/project/abriplot/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)

---

## What it does

`abriplot` takes the TSV files produced by abricate (CARD and VFDB databases) and a Newick phylogenetic tree, and produces a single publication-ready figure with:

- **Phylogenetic tree** — branch lengths proportional, samples ordered by topology
- **AMR gene heatmap** — presence/absence matrix from CARD, genes sorted by prevalence
- **Drug-class annotation strip** — one row per resistance class above the AMR block; multi-class genes appear in all relevant rows; classes auto-detected from your data
- **Virulence gene heatmap** — presence/absence matrix from VFDB
- PDF + PNG exported automatically

![Example figure](docs/example_figure.png)

---

## Installation

```bash
pip install abriplot
```

No Biopython required — pure Python Newick parser included.

---

## Quick start

```bash
abriplot \
    --input_dir  /path/to/abricate_results/ \
    --tree       /path/to/tree.nwk \
    --output     figure.pdf
```

This produces `figure.pdf` and `figure.png` in the current directory.

---

## All options

```
usage: abriplot [-h] --input_dir DIR --tree FILE [--output FILE] [--title TEXT]
                [--min_identity FLOAT] [--min_coverage FLOAT]
                [--top_amr INT] [--top_vf INT]

required arguments:
  --input_dir DIR       Directory containing *__card.tsv and *__vfdb.tsv files
  --tree FILE           Newick tree file (.nwk / .tree)

output:
  --output FILE         Output path — pdf, png, or svg [default: abriplot_figure.pdf]
  --title TEXT          Figure title [default: 'Phylogenetic Distribution of AMR
                        and Virulence Genes']

filtering:
  --min_identity FLOAT  Minimum % identity threshold [default: 80]
  --min_coverage FLOAT  Minimum % coverage threshold [default: 80]
  --top_amr INT         Keep only the N most prevalent AMR genes [default: all]
  --top_vf INT          Keep only the N most prevalent virulence genes [default: all]
```

### Examples

```bash
# Custom title for a publication
abriplot --input_dir results/ --tree tree.nwk --output fig1.pdf \
         --title "AMR gene distribution in V. cholerae (2002–2019)"

# Tighter figure — top 40 AMR genes, top 25 VF genes
abriplot --input_dir results/ --tree tree.nwk --output fig1.pdf \
         --top_amr 40 --top_vf 25

# Stricter filtering
abriplot --input_dir results/ --tree tree.nwk --output fig1.pdf \
         --min_identity 95 --min_coverage 90
```

---

## Input file format

abricate output files must follow the naming convention:

```
<sample>__card.tsv
<sample>__vfdb.tsv
```

These are produced automatically when you run abricate with the `--db card` and `--db vfdb` flags. Example abricate command:

```bash
abricate --db card  assembly.fasta > sample__card.tsv
abricate --db vfdb  assembly.fasta > sample__vfdb.tsv
```

---

## Use as a library

```python
from abriplot.io     import load_all
from abriplot.newick import tip_order
from abriplot.plot   import (build_matrix, order_by_tree,
                              sort_genes_by_prevalence, make_figure)

with open("tree.nwk") as f:
    ordered_tips, root, leaves = tip_order(f.read())

amr_records, vf_records = load_all("results/", min_id=80, min_cov=80)

amr_mat, amr_meta = build_matrix(amr_records)
vf_mat,  _        = build_matrix(vf_records)

amr_mat = order_by_tree(amr_mat, ordered_tips)
vf_mat  = order_by_tree(vf_mat,  ordered_tips).reindex(amr_mat.index, fill_value=0)

amr_mat = sort_genes_by_prevalence(amr_mat)
vf_mat  = sort_genes_by_prevalence(vf_mat)

make_figure(amr_mat, vf_mat, amr_meta, root, leaves,
            amr_mat.index.tolist(), "figure.pdf",
            title="My custom title")
```

---

## Requirements

- Python ≥ 3.8
- numpy ≥ 1.21
- pandas ≥ 1.3
- matplotlib ≥ 3.5


---

## License

MIT © [Your Name]
