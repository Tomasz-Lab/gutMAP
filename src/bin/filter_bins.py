#!/usr/bin/env python3
"""Copy the CheckM2 bins with high completeness and low contamination to a new directory."""
import argparse
import os
import shutil

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkm_report")
    parser.add_argument("bins_dir")
    parser.add_argument("out_dir")
    parser.add_argument("min_completeness", type=float)
    parser.add_argument("max_contamination", type=float)
    args = parser.parse_args()

    df = pd.read_csv(args.checkm_report, sep="\t", header=0, index_col=None, on_bad_lines="skip")
    df = df.dropna(subset=["Completeness", "Contamination"])
    hq = df[(df["Completeness"] >= args.min_completeness) & (df["Contamination"] <= args.max_contamination)]

    os.makedirs(args.out_dir, exist_ok=True)
    for bin_id in hq["Name"]:
        source = os.path.join(args.bins_dir, f"{bin_id}.fa")
        if os.path.exists(source):
            shutil.copy(source, os.path.join(args.out_dir, f"{bin_id}.fa"))


if __name__ == "__main__":
    main()
