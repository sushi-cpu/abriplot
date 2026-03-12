"""
Basic tests for abriplot.
Run with:  pytest tests/ -v
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pandas as pd
import pytest

from abriplot.newick import parse_newick, assign_x, assign_y, tip_order
from abriplot.plot   import (resistance_classes, build_matrix,
                              order_by_tree, sort_genes_by_prevalence,
                              build_class_palette)
from abriplot.io     import extract_sample


# ── Newick parser ──────────────────────────────────────────────────────────────

NWK_SIMPLE = "((A:0.1,B:0.2):0.3,C:0.4);"
NWK_REAL   = ("(IDH_11763_2018_:0.0003,(IDH_7451_2015_:0.018,"
              "L_98411_2006_:0.019)100:0.004);")

def test_parse_newick_simple():
    root = parse_newick(NWK_SIMPLE)
    assign_x(root)
    leaves = assign_y(root)
    names = [lf.name for lf in leaves]
    assert set(names) == {"A", "B", "C"}

def test_tip_order_count():
    names, root, leaves = tip_order(NWK_SIMPLE, strip_trailing_underscore=False)
    assert len(names) == 3

def test_tip_order_strip_underscore():
    names, _, _ = tip_order(NWK_REAL, strip_trailing_underscore=True)
    assert all(not n.endswith("_") for n in names)

def test_tip_order_no_strip():
    names, _, _ = tip_order(NWK_REAL, strip_trailing_underscore=False)
    assert any(n.endswith("_") for n in names)


# ── resistance_classes ─────────────────────────────────────────────────────────

def test_resistance_classes_card_style():
    result = resistance_classes("fluoroquinolone;macrolide;penam")
    assert result == {"fluoroquinolone", "macrolide", "penam"}

def test_resistance_classes_empty():
    assert resistance_classes("") == set()
    assert resistance_classes(None) == set()

def test_resistance_classes_single():
    assert resistance_classes("aminoglycoside") == {"aminoglycoside"}

def test_resistance_classes_whitespace():
    result = resistance_classes(" fluoroquinolone ; macrolide ")
    assert "fluoroquinolone" in result
    assert "macrolide" in result


# ── build_matrix ───────────────────────────────────────────────────────────────

RECORDS = [
    {"sample": "S1", "gene": "blaTEM", "resistance": "penam",           "identity": 99, "coverage": 100},
    {"sample": "S1", "gene": "gyrA",   "resistance": "fluoroquinolone", "identity": 95, "coverage": 100},
    {"sample": "S2", "gene": "blaTEM", "resistance": "penam",           "identity": 98, "coverage": 100},
    {"sample": "S3", "gene": "tetA",   "resistance": "tetracycline",    "identity": 99, "coverage": 100},
]

def test_build_matrix_shape():
    mat, meta = build_matrix(RECORDS)
    assert mat.shape == (3, 3)   # 3 samples × 3 genes

def test_build_matrix_binary():
    mat, _ = build_matrix(RECORDS)
    assert set(mat.values.flatten().tolist()).issubset({0, 1})

def test_build_matrix_presence():
    mat, _ = build_matrix(RECORDS)
    assert mat.loc["S1", "blaTEM"] == 1
    assert mat.loc["S3", "blaTEM"] == 0

def test_build_matrix_meta():
    _, meta = build_matrix(RECORDS)
    assert "blaTEM" in meta
    assert "penam" in meta["blaTEM"]


# ── order_by_tree ──────────────────────────────────────────────────────────────

def test_order_by_tree_ordering():
    mat, _ = build_matrix(RECORDS)
    ordered = order_by_tree(mat, ["S3", "S1", "S2"])
    assert ordered.index.tolist() == ["S3", "S1", "S2"]

def test_order_by_tree_missing_tip():
    mat, _ = build_matrix(RECORDS)
    # S4 is not in the data — should be silently skipped
    ordered = order_by_tree(mat, ["S1", "S2", "S3", "S4"])
    assert "S4" not in ordered.index


# ── sort_genes_by_prevalence ───────────────────────────────────────────────────

def test_sort_genes_by_prevalence():
    mat, _ = build_matrix(RECORDS)
    sorted_mat = sort_genes_by_prevalence(mat)
    # blaTEM appears in 2 samples, should come first
    assert sorted_mat.columns[0] == "blaTEM"


# ── build_class_palette ────────────────────────────────────────────────────────

def test_build_class_palette_completeness():
    _, meta = build_matrix(RECORDS)
    palette = build_class_palette(meta)
    assert "penam"           in palette
    assert "fluoroquinolone" in palette
    assert "tetracycline"    in palette

def test_build_class_palette_distinct_colours():
    _, meta = build_matrix(RECORDS)
    palette = build_class_palette(meta)
    colours = list(palette.values())
    assert len(colours) == len(set(colours))   # all colours unique


# ── extract_sample ─────────────────────────────────────────────────────────────

def test_extract_sample_card():
    assert extract_sample("IDH_4533_2013__card.tsv") == "IDH_4533_2013"

def test_extract_sample_vfdb():
    assert extract_sample("SC_103_2003__vfdb.tsv") == "SC_103_2003"

def test_extract_sample_path():
    assert extract_sample("/data/results/IDH_4533_2013__card.tsv") == "IDH_4533_2013"
