"""
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) CHER_WEI_YUAN, 21JAN2021 
License     : MIT_LICENSE 
Maintainer  : E0031403@U.NUS.EDU 
Portability : POSIX

This programme annotates VCF converted to TSV format by GATK's VariantsToTable 
with gNOMAD v3.1.2 gene-specific data, filters the TSV based on user-supplied
tags and cleans & plot the final output
"""

from argparse import ArgumentParser, BooleanOptionalAction
import dask.dataframe as dd
from numpy import nan
import sys

PROGRAM_NAME = "finis"

def exit_with_error(message, exit_status):
    """Print an error message to stderr, prefixed by the program name and "ERROR".
    Then exit program with supplied exit status.
    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    """
    error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)

def parse_args():
    """
    Parse commandline arguments.
    Returns object command with command line values as attributes
    """
    parser = ArgumentParser(description = "Processes VEP VCF with CSQ header")
    parser.add_argument("--vcf_tsv",
                        type = str,
                        required = True,
                        help = "Sample VCF coverted to tsv by GATK " +\
                               "VariantsToTable")
    parser.add_argument("--out_full",
                        type = str,
                        required = True,
                        help = "Output path and file prefix, e.g. /mnt/e/out " +\
                               "Set to '0' if full output is not required")
    parser.add_argument("--out_simple",
                        type = str,
                        required = True,
                        help = "Output path and file prefix, e.g. /mnt/e/out " +\
                               "Set to '0' if simple output is not required")
    parser.add_argument("--skip_spliceai",
                        action = BooleanOptionalAction,
                        help = "Add flag to not process SpliceAI column")
    parser.add_argument("--column_name",
                        type = str,
                        default = "SpliceAI_SpliceAI",
                        help = "Column name to process")
    parser.add_argument("--chunksize",
                        type = int,
                        default = 1000000,
                        help = "Rows to write at a time")

    return parser.parse_args()

def removeComments(inputFileName, outputFileName):
    input = open(inputFileName, "r")
    output = open(outputFileName, "w")
    removed = open(outputFileName + ".linesRemoved", "w")
    output.write(input.readline())
    for line in input:
        if not line.lstrip().startswith("##"):
            output.write(line)
        else:
            removed.write(line)
    input.close()
    output.close()
    removed.close()

def str_to_numeric(string_number):
    """
    Convert float from string to numeric
    """
    # Some SpliceAI entries have "." instead of probabilities
    if string_number == ".":
        return 0
    try:
        return float(string_number)
    except:
        return int(string_number)

def parse_splice_predictions(vcfdf, column_name):
    """
    Parse SpliceAI scores
    Store highest score obtained per variant across all score categories
    """
    # Replace cells with "-" as nan
    vcfdf.replace({"-": nan}, inplace = True)
    vcfdf[column_name + "_highest_score"] = nan
    for index, row in vcfdf.iterrows():
        ## Processing SpliceAI column
        # Each row in SpliceAI column contains one SpliceAI record set per 
        # variant, which means there can be many record sets. 
        # Each score set is separated by comma
        # e.g. GCTCTCT|SLC25A13|0.00|0.00|0.00|0.00|6|-45|0|-45,
        #      GCTCTCTCT|SLC25A13|0.00|0.00|0.00|0.00|30|-45|0|-45,
        #      GCTCTCTCTCT|SLC25A13|0.00|0.00|0.00|0.00|-12|-45|48|-45
        row_high_scores = []
        gene_name = vcfdf.loc[index, "SYMBOL"]
        try:
            if dd.isnull(vcfdf.loc[index, column_name]):
                continue
            variant_records = vcfdf.loc[index, column_name].split(",")
            for record_set in variant_records:
                record = record_set.split("|")
                # Ensure highest score for each row is for the correct gene
                if record[1] != gene_name:
                    continue
                row_high_scores += [max([str_to_numeric(record[2]),   
                                        # DS_AG score (acceptor gain)
                                        str_to_numeric(record[3]),   
                                        # DS_AL score (acceptor loss)
                                        str_to_numeric(record[4]),   
                                        # DS_DG score (donor gain)
                                        str_to_numeric(record[5])])] 
                                        # DS_DL score (donor loss)
        # Catch error when SpliceAI entry is nan
        except AttributeError:
            pass

        if row_high_scores:
            vcfdf.loc[index, column_name + "_highest_score"] = max(row_high_scores)
        else:
            pass
    return vcfdf

def main():
    # Parse command-line arguments
    args = parse_args()
    
    # Parse remaining command-line arguments
    vcf_tsv = args.vcf_tsv
    out_full = args.out_full
    out_simple = args.out_simple
    skip_spliceai = args.skip_spliceai
    column_name = args.column_name
    chunksize = args.chunksize
    
    # Load dataframes and remove lines starting with double hex
    removeComments(vcf_tsv, vcf_tsv + ".noHex.tsv")
    try:
        vcf_df = dd.read_csv(vcf_tsv + ".noHex.tsv", sep = "\t", header = 1)
        if skip_spliceai:
            pass
        else:
            vcf_df = parse_splice_predictions(vcf_df, column_name)
    except KeyError:
        vcf_df = dd.read_csv(vcf_tsv + ".noHex.tsv", sep = "\t")
        if skip_spliceai:
            pass
        else:
            vcf_df = parse_splice_predictions(vcf_df, column_name)     
           
    # Output full csv
    if out_full != "0":
        vcf_df.to_csv(out_full, index = False)

    # Output simplified csv
    if out_simple != "0":
        if skip_spliceai:
            vcf_df = vcf_df[[
                "HGVSg",
                "HGVSc",
                "HGVSp",
                "SYMBOL",
                "ZYG",
                "Existing_variation",
                "Consequence",
                "gnomADg_AF",
                "gnomADg_EAS_AF",
                "ClinVar_CLNSIG"]]
        else:
            vcf_df = vcf_df[[
                "HGVSg",
                "HGVSc",
                "HGVSp",
                "SYMBOL",
                "ZYG",
                "Existing_variation",
                "Consequence",
                "gnomADg_AF",
                "gnomADg_EAS_AF",
                "ClinVar_CLNSIG",
                "SpliceAI_SpliceAI_highest_score"]]
    vcf_df.to_csv(out_simple, index = False)
    
# Run main() if this script is called from command-line    
if __name__ == "__main__":
    main()