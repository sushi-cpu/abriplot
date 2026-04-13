"""
abriplot.io
===========
Parsers for abricate TSV output files.
"""

import os
import re
import glob

import pandas as pd


def extract_sample(filepath):
    """Derive a clean sample name from an abricate TSV filename."""
    return re.sub(r"_(card|resfinder|vfdb)\.tsv$", "",
                  os.path.basename(filepath))


def parse_tsv(fp, min_id=80.0, min_cov=80.0):
    """
    Read one abricate TSV file, apply identity/coverage filters,
    and return a DataFrame or None if no rows pass.
    """
    try:
        df = pd.read_csv(fp, sep="\t", low_memory=False)
        df.columns = [c.lstrip("#") for c in df.columns]
        if df.empty or "GENE" not in df.columns:
            return None
        df["%IDENTITY"] = pd.to_numeric(df["%IDENTITY"], errors="coerce")
        df["%COVERAGE"] = pd.to_numeric(df["%COVERAGE"], errors="coerce")
        df = df[(df["%IDENTITY"] >= min_id) & (df["%COVERAGE"] >= min_cov)]
        return df if not df.empty else None
    except Exception as e:
        print(f"  WARNING: could not parse {fp}: {e}")
        return None


def load_all(input_dir, min_id=80.0, min_cov=80.0):
    """
    Load all *__card.tsv and *__vfdb.tsv files from input_dir.

    Returns:
        amr_records : list of dicts  {sample, gene, resistance, identity, coverage}
        vf_records  : list of dicts  {sample, gene, product,    identity, coverage}
    """
    amr, vf = [], []

    for fp in sorted(glob.glob(os.path.join(input_dir, "*__card.tsv"))):
        df = parse_tsv(fp, min_id, min_cov)
        if df is None:
            continue
        s = extract_sample(fp)
        for _, r in df.iterrows():
            amr.append({
                "sample":     s,
                "gene":       str(r["GENE"]).strip(),
                "resistance": (str(r.get("RESISTANCE", ""))
                               if pd.notna(r.get("RESISTANCE")) else ""),
                "identity":   r["%IDENTITY"],
                "coverage":   r["%COVERAGE"],
            })

    for fp in sorted(glob.glob(os.path.join(input_dir, "*__vfdb.tsv"))):
        df = parse_tsv(fp, min_id, min_cov)
        if df is None:
            continue
        s = extract_sample(fp)
        for _, r in df.iterrows():
            vf.append({
                "sample":   s,
                "gene":     str(r["GENE"]).strip(),
                "product":  str(r.get("PRODUCT", "")).strip(),
                "identity": r["%IDENTITY"],
                "coverage": r["%COVERAGE"],
            })

    print(f"  CARD records : {len(amr)}")
    print(f"  VFDB records : {len(vf)}")
    return amr, vf
