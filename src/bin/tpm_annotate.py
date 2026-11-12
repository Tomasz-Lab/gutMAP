#!/usr/bin/env python3
"""Add read counts and TPM from samtools idxstats to a Bakta annotation TSV."""
import argparse
from io import StringIO

import pandas as pd


def read_idxstats(path):
    genes = []
    with open(path) as f:
        for line in f:
            gene_id, length, mapped_reads, _ = line.strip().split("\t")
            if gene_id == "*" or int(length) == 0:
                continue
            genes.append((gene_id, int(length), int(mapped_reads)))
    return genes


def tpm_table(genes):
    rpk = [reads / (length / 1000.0) for _, length, reads in genes]
    # Not sum(): since Python 3.12 it uses compensated summation, and the TPM values then
    # differ from the production run in the last digit.
    rpk_sum = 0.0
    for value in rpk:
        rpk_sum += value
    scaling_factor = rpk_sum / 1_000_000.0 if rpk_sum > 0 else 0
    return pd.DataFrame([
        {"gene_id": gene_id, "length": length, "read_count": reads,
         "tpm": gene_rpk / scaling_factor if scaling_factor > 0 else 0}
        for (gene_id, length, reads), gene_rpk in zip(genes, rpk)
    ], columns=["gene_id", "length", "read_count", "tpm"])


def read_bakta_tsv(path):
    """Bakta writes comment lines first. The last '#' line is the header."""
    with open(path) as f:
        lines = f.readlines()
    header_idx = max(i for i, line in enumerate(lines) if line.startswith("#"))
    header = lines[header_idx].lstrip("#").strip()
    table = header + "\n" + "".join(lines[header_idx + 1:])
    return pd.read_csv(StringIO(table), sep="\t", skipinitialspace=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("idxstats")
    parser.add_argument("bakta_tsv")
    parser.add_argument("output")
    args = parser.parse_args()

    abundances = tpm_table(read_idxstats(args.idxstats))
    annotation = read_bakta_tsv(args.bakta_tsv)
    merged = annotation.merge(abundances, left_on="Locus Tag", right_on="gene_id", how="inner")
    merged.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
