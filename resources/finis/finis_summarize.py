#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Summarize simplified results of all samples in the pipeline
"""

from argparse import ArgumentParser
from os import listdir
import pandas as pd

def parse_args():
    """
    Parse commandline arguments.
    Returns object command with command line values as attributes
    """
    parser = ArgumentParser(description = "Summarize results for all samples")
    parser.add_argument("-i", "--indir",
                        type     = str,
                        required = True,
                        default  = "results/variants/simplified_data",
                        help     = "Input directory to simplified results")
    parser.add_argument("-o", "--outdir",
                        type     = str,
                        default  = "results/",
                        help     = "Output directory")
    parser.add_argument("-g", "--genes", 
                        nargs    = "*",
                        default  = None,
                        help     = "Genes specified in SYMBOL column to keep")
    return parser.parse_args()

def main():
    # Parse command-line arguments
    args  = parse_args()
    
    # Parse remaining command-line arguments
    indir  = args.indir
    outdir = args.outdir
    genes  = args.genes

    # Concatenate all simplified data dataframes
    df_list = []
    for csv in listdir(indir):
        prefix = csv.split("_")[0]
        df     = pd.read_csv(f"{indir}/{csv}")
        df["sample_name"] = prefix
        df.drop_duplicates(subset = ["HGVSg"], inplace = True, ignore_index = True)
        df_list.append(df)
        df     = None

    merged_df = pd.concat(df_list, axis = 0, ignore_index = True)

    # Iterate through rows and decide which rows to keep
    keep_indices = []
    for index, row in merged_df.iterrows():
        clinvar  = merged_df.loc[index, "ClinVar_CLNSIG"]
        spliceai = merged_df.loc[index, "SpliceAI_SpliceAI_highest_score"]
        if not pd.isna(clinvar):
            if "pathogenic" in clinvar.lower():
                keep_indices.append(index)
        if spliceai >= 0.2:
            keep_indices.append(index)

    keep_indices = list(set(keep_indices))
    merged_df = merged_df.iloc[keep_indices]
    merged_df = merged_df[merged_df.SYMBOL.isin(genes)]

    # Export
    merged_df.to_csv(f"{outdir}/summarized_results.csv", index = False)
    
    # Return 0 for successful run
    return 0

# Run main() if this script is called from command-line    
if __name__ == "__main__":
    main()