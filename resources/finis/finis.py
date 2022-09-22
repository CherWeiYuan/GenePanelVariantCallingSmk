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

from argparse import ArgumentParser
from logging import basicConfig
from logging import info
from logging import error
from logging import DEBUG
from numpy import nan
import pandas as pd
import plotly.express as px
import sys

PROGRAM_NAME = "finis"
FILTER_TAG_KEY_ERROR = 1
FILTER_TAG_VALUE_ERROR = 2

def init_logging(log_filename):
    """If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv
    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    """
    if log_filename is not None:
        basicConfig(filename=log_filename,
                            level=DEBUG,
                            filemode="w",
                            format="%(asctime)s %(levelname)s - %(message)s",
                            datefmt="%Y-%m-%dT%H:%M:%S%z")
        info("program started")
        info("command line: %s", " ".join(sys.argv))

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
    parser = ArgumentParser(description = "Annotates a VCF using gNOMAD gene-" +\
                                          "specific data")
    parser.add_argument("--gnomadCSV",
                        type = str,
                        required = True,
                        help = "gNOMAD csv can be downloaded using region " +\
                               "view, e.g. https://gnomad.broadinstitute.org/" +\
                               "region/1-93992792-94121192?dataset=gnomad_r3")
    parser.add_argument("--vcf_tsv",
                        type = str,
                        required = True,
                        help = "Sample VCF coverted to tsv by GATK " +\
                               "VariantsToTable")
    parser.add_argument("--csq_headers",
                        type = str,
                        default = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|SOURCE|HGVS_OFFSET|HGVSg|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|ClinVar|ClinVar_ClinVar|ClinVar_vcf|ClinVar_exact|ClinVar_0|ClinVar_AF_ESP|ClinVar_AF_EXAC|ClinVar_AF_TGP|ClinVar_ALLELEID|ClinVar_CLNDN|ClinVar_CLNDNINCL|ClinVar_CLNDISDB|ClinVar_CLNDISDBINCL|ClinVar_CLNHGVS|ClinVar_CLNREVSTAT|ClinVar_CLNSIG|ClinVar_CLNSIGCONF|ClinVar_CLNSIGINCL|ClinVar_CLNVC|ClinVar_CLNVCSO|ClinVar_CLNVI|ClinVar_DBVARID|ClinVar_GENEINFO|ClinVar_MC|ClinVar_ORIGIN|ClinVar_RS|ClinVar_SSR|NARD|NARD_AF|NARD_AF_MNG|NARD_AF_KOR|NARD_AF_JPN|NARD_AF_CHN|NARD_AF_HKG",
                        help = "INFO field for VEP's CSQ, e.g. " +\
                               "'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|SOURCE|HGVS_OFFSET|HGVSg|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|ClinVar|ClinVar_ClinVar|ClinVar_vcf|ClinVar_exact|ClinVar_0|ClinVar_AF_ESP|ClinVar_AF_EXAC|ClinVar_AF_TGP|ClinVar_ALLELEID|ClinVar_CLNDN|ClinVar_CLNDNINCL|ClinVar_CLNDISDB|ClinVar_CLNDISDBINCL|ClinVar_CLNHGVS|ClinVar_CLNREVSTAT|ClinVar_CLNSIG|ClinVar_CLNSIGCONF|ClinVar_CLNSIGINCL|ClinVar_CLNVC|ClinVar_CLNVCSO|ClinVar_CLNVI|ClinVar_DBVARID|ClinVar_GENEINFO|ClinVar_MC|ClinVar_ORIGIN|ClinVar_RS|ClinVar_SSR|NARD|NARD_AF|NARD_AF_MNG|NARD_AF_KOR|NARD_AF_JPN|NARD_AF_CHN|NARD_AF_HKG'")
    parser.add_argument("--filter_tags",
                        nargs = "+",
                        required = True,
                        help = """
                        Filtering criteria in string format.
                        
                        Rules:
                        -If filtering is not required, do not use --filter_tags
                        -Each tag follows a Pandas.Dataframe.query string format 
                        -Each filtering tag needs to be in quotation marks; if the condition is a string, it needs it's own quotation marks; integers or float does not need quotation marks; e.g. "tag1 == 'string'" "tag2 < 0.05" \
                        -Backquotes must be used for long VCF fields with white space, such as "\`Allele Frequency East Asian\` < 0.01"
                        -To use more than one tag, do --filter_tags "X < 0.01" "Y < 0.01" "Z < 0.01"
                        -If you want to filter to obtain entries absent in population database, you can do "X == @nan"
                        
                        Examples:
                            Absence rule: --filter_tags "\`Allele Frequency East Asian\` == @nan"
                            1% rule: --filter_tags "\`Allele Frequency East Asian\` < 0.01"
                            ClinVar: "ClinVar_CLNSIG == 'Likely pathogenic' or ClinVar_CLNSIG == 'Pathogenic'"
                            Combined: --filter_tags "\`Allele Frequency East Asian\` < 0.01" "`Allele Frequency East Asian` == @nan" "CLNSIG == 'Likely pathogenic' or CLNSIG == 'Pathogenic'"
                        """)
    parser.add_argument("--out",
                        type = str,
                        default = "out",
                        help = "Output path and file prefix, e.g. /mnt/e/out")
    return parser.parse_args()

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

def parse_csq(vcfdf, csq):
    """
    Parse CSQ column
    """
    vcfdf.rename(columns={"CHROM": "Chromosome", "POS": "Position", 
                          "REF": "Reference", "ALT": "Alternate",
                          "CSQ": csq}, inplace = True)
    
    # Get CSQ column names
    csq_names = vcfdf[csq].name.split("|")
    
    # Get length of CSQ columns
    csq_len = len(csq_names)
    
    # Get CSQ values as a list per row
    csq_df = vcfdf[csq].str.split("|").to_frame()
    
    # VEP CSQ has excess columns, and they do not have names
    # We have to discard the columns without names,
    # which occur immediately after the columns with names 
    # Note that each row is still a list, where each element is equivalent
    # to a CSV column
    for index, row in csq_df.iterrows():
        csq_df.loc[index, csq] = csq_df.loc[index, csq][:csq_len]
    
    # Create dataframe with only CSQ data
    csq_df = pd.DataFrame(csq_df[csq].to_list(), columns = csq_names)

    return csq_df

def parse_splice_predictions(vcfdf, gene_names):
    """
    Parse SpliceAI scores
    Store highest score obtained per variant across all score categories
    """
    vcfdf["SpliceAI_highest_score"] = 0

    for index, row in vcfdf.iterrows():
        ## Processing SpliceAI column
        # Each row in SpliceAI column contains one SpliceAI record set per 
        # variant, which means there can be many record sets. 
        # Each score set is separated by comma
        # e.g. GCTCTCT|SLC25A13|0.00|0.00|0.00|0.00|6|-45|0|-45,
        #      GCTCTCTCT|SLC25A13|0.00|0.00|0.00|0.00|30|-45|0|-45,
        #      GCTCTCTCTCT|SLC25A13|0.00|0.00|0.00|0.00|-12|-45|48|-45
        row_high_scores = []
        try:
            if pd.isnull(vcfdf.loc[index, "SpliceAI"]):
                continue
            variant_records = vcfdf.loc[index, "SpliceAI"].split(",")
            for record_set in variant_records:
                record = record_set.split("|")
                gene = record[1]
                if gene in gene_names:
                    row_high_scores += [max([str_to_numeric(record[2]),   
                                            # DS_AG score (acceptor gain)
                                            str_to_numeric(record[3]),   
                                            # DS_AL score (acceptor loss)
                                            str_to_numeric(record[4]),   
                                            # DS_DG score (donor gain)
                                            str_to_numeric(record[5])])] 
                                            # DS_DL score (donor loss)
                else:
                    row_high_scores += [0]
                    info(f"SpliceAI produced entry for gene {gene} which is " +\
                         f"found in the VCF's gene list: {gene_names}")
        # Catch error when SpliceAI entry is nan
        except AttributeError:
            pass
        vcfdf.loc[index, "SpliceAI_highest_score"] = max(row_high_scores)
        
    return vcfdf

def parse_vcf(vcf_tsv, csq):
    """
    Parameters
    ----------
    vcf_tsv : TYPE str
        DESCRIPTION. Path to vcf.tsv file to be annotated
    csq : TYPE str
        DESCRIPTION. CSQ INFO field in VEP-annotated VCF file     

    Returns
    -------
    vcfdf : TYPE Pandas dataframe
        DESCRIPTION. Contains VEP-annotated vcf entries in dataframe format
    """
    vcfdf = pd.read_csv(vcf_tsv, sep = "\t")

    # Parse csq column
    csq_df = parse_csq(vcfdf, csq)    

    # Add csq_df columns to vcfdf since row order of csq_df is 
    # not changed from manipulating vcfdf
    vcfdf = pd.concat([vcfdf, csq_df], axis = 1)

    # Parse SpliceAI and Pangolin columns
    gene_names = vcfdf["GENE"].unique()
    vcfdf = parse_splice_predictions(vcfdf, gene_names)

    # Change all blank cells into np.nan
    vcfdf = vcfdf.replace(r"^\s*$", nan, regex = True)
    
    # Set column type for comparison with gNOMAD genomes csv
    vcfdf = vcfdf.astype({"Chromosome": str, 
                          "Position" : int,
                          "Reference": str,
                          "Alternate": str})
    
    return vcfdf

def parse_gNOMAD(gnomadcsv):
    """
    Parameters
    ----------
    gnomadcsv : TYPE str
        DESCRIPTION. Path to gNOMAD csv file
        
    Returns
    -------
    gnomadOut : TYPE Pandas dataframe
        DESCRIPTION. Entries found in both vcf and gNOMAD
    """
    gnomadf = pd.read_csv(gnomadcsv, low_memory=False)
    gnomadf = gnomadf.astype({"Chromosome": str, "Position" : int,
                              "Reference": str, "Alternate": str,
                              "Allele Count": int, "Allele Number": int})
    return gnomadf

def update_allele_freq(gnomadf):
    """
    Calculates allele frequency for every population and add them as new columns
    """
    population_list = ["Other", "Latino/Admixed American", "European (Finnish)",
                       "Amish", "East Asian", "Middle Eastern", "South Asian",
                       "African/African American", "Ashkenazi Jewish",
                       "European (non-Finnish)"]
    for pop in population_list:
        gnomadf[f"Allele Frequency {pop}"] = \
            gnomadf[f"Allele Count {pop}"]/ gnomadf[f"Allele Number {pop}"]
    return gnomadf

def merge(vcfdf, gnomadf):
    """
    Parameters
    ----------
    vcfdf : TYPE Pandas dataframe
        DESCRIPTION. Contains vcf entries
    gnomadf : TYPE Pandas dataframe
        DESCRIPTION. Contains gNOMAD entries

    Returns
    -------
    gnomadOut : TYPE Pandas dataframe
        DESCRIPTION. Entries found in both vcf and gNOMAD
    """
    gnomadgenome_col = ["Chromosome", "Position","Reference","Alternate",
                        "rsIDs","Source","Filters - exomes","Filters - genomes",
                        "Transcript","HGVS Consequence","Protein Consequence",
                        "Transcript Consequence","VEP Annotation",
                        "ClinVar Clinical Significance","ClinVar Variation ID",
                        "Flags","Allele Count","Allele Number","Allele Frequency",
                        "Homozygote Count","Hemizygote Count","Allele Count Other",
                        "Allele Number Other","Homozygote Count Other",
                        "Hemizygote Count Other",
                        "Allele Count Latino/Admixed American",
                        "Allele Number Latino/Admixed American",
                        "Homozygote Count Latino/Admixed American",
                        "Hemizygote Count Latino/Admixed American",
                        "Allele Count European (Finnish)",
                        "Allele Number European (Finnish)",
                        "Homozygote Count European (Finnish)",
                        "Hemizygote Count European (Finnish)",
                        "Allele Count Amish","Allele Number Amish",
                        "Homozygote Count Amish","Hemizygote Count Amish",
                        "Allele Count East Asian","Allele Number East Asian",
                        "Homozygote Count East Asian","Hemizygote Count East Asian",
                        "Allele Count Middle Eastern","Allele Number Middle Eastern",
                        "Homozygote Count Middle Eastern","Hemizygote Count Middle Eastern",
                        "Allele Count African/African American",
                        "Allele Number African/African American",
                        "Homozygote Count African/African American",
                        "Hemizygote Count African/African American",
                        "Allele Count South Asian","Allele Number South Asian",
                        "Homozygote Count South Asian","Hemizygote Count South Asian",
                        "Allele Count Ashkenazi Jewish","Allele Number Ashkenazi Jewish",
                        "Homozygote Count Ashkenazi Jewish","Hemizygote Count Ashkenazi Jewish",
                        "Allele Count European (non-Finnish)","Allele Number European (non-Finnish)",
                        "Homozygote Count European (non-Finnish)","Hemizygote Count European (non-Finnish)",
                        "Allele Frequency Other", "Allele Frequency Latino/Admixed American",
                        "Allele Frequency European (Finnish)", "Allele Frequency Amish", "Allele Frequency East Asian",
                        "Allele Frequency Middle Eastern","Allele Frequency African/African American", 
                        "Allele Frequency South Asian","Allele Frequency Ashkenazi Jewish", 
                        "Allele Frequency European (non-Finnish)"]
    merged_df = pd.merge(vcfdf,
                         gnomadf[gnomadgenome_col],
                         on = ["Chromosome", "Position", 
                               "Reference", "Alternate"], 
                         how = "left")
    return merged_df

def AF_to_numeric(merged_df):
    """
    Changes all allele frequency fields to numeric type
    """
    merged_df["NARD_AF_HKG"] = merged_df["NARD_AF_HKG"].str.split(",")\
                               .map(lambda x: float(x[0]))
    AF = ["AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF", "NARD_AF",
          "NARD_AF_MNG", "NARD_AF_KOR", "NARD_AF_JPN" "NARD_AF_CHN", 
          "NARD_AF_HKG", "gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", 
          "gnomAD_ASJ_AF", "gnomAD_EAS_AF", "gnomAD_FIN_AF", "gnomAD_NFE_AF", 
          "gnomAD_OTH_AF", "gnomAD_SAS_AF", "Allele Frequency European (Finnish)", 
          "Allele Frequency Amish",  "Allele Frequency East Asian",
          "Allele Frequency Middle Eastern", "Allele Frequency South Asian"
          "Allele Frequency African/African American", 
          "Allele Frequency Ashkenazi Jewish", 
          "Allele Frequency European (non-Finnish)"]
    try:
        merged_df[AF] = merged_df[AF].apply(pd.to_numeric)
    except:
        info("WARNING: Allele frequency dtype changing has failed.")
    return merged_df
        
def filter_tsv(naf_merged_df, filter_tag):
    """
    Parameters
    ----------
    naf_merged_df : TYPE Pandas dataframe
        DESCRIPTION. Dataframe containing VEP-gNOMAD-exomes and gNOMAD genome
                     annotation from merge(vcfdf, gnomadf). Allele frequencies
                     are changed to numeric dtype by AF_to_numeric(merged_df).
    filter_tag : TYPE str
        DESCRIPTION. Filter tags in the format usable by Pandas DataFrame.query

    Returns
    -------
    naf_merged_df : TYPE TYPE Pandas dataframe
    """
    try:
        naf_merged_df.query(filter_tag, inplace = True) 
    except KeyError:
        exit_with_error(f"Your filter tag {filter_tag} is faulty. " +\
                        "Check if it is correct (make sure there is no " +\
                        "white space, correct spelling and be CAPs sensitive", 
                        FILTER_TAG_KEY_ERROR)
    except ValueError:
        exit_with_error("Filter tag {filter_tag} is faulty. "+\
                        "Check if you enclose one entire tag" +\
                        'with double quotation (i.e. {rule: "tag ' +\
                        '< int"} or {rule: "tag name" == "string"})',
                        FILTER_TAG_VALUE_ERROR)
    return naf_merged_df

def csv_out(merged_df, out):
    """
    Output gNOMAD genome annotated tsv
    
    Also output gNOMAD genome csv if allele frequency is updated
    """
    merged_df.to_csv(f"{out}_gnomadAnnotated.csv", index = False)

def plot(df, filter_name, out): 
    df["Consequence"] = df["Consequence"].replace(nan, "None")
    df["ClinVar_CLNSIG"] = df["ClinVar_CLNSIG"].replace(nan, "None")
    fig = px.histogram(df, x = "Consequence", color = "ClinVar_CLNSIG", 
                       title = "VCF.TSV Summary", template = "seaborn",
                       labels = {"ClinVar_CLNSIG": "ClinVar significance", 
                                 "Consequence": "Variant consequence"},)
    fig.write_html(f"{out}_{filter_name}.html")

def main():
    # Parse command-line arguments
    args = parse_args()
    
    # Parse remaining command-line arguments
    gnomadcsv = args.gnomadCSV
    vcf_tsv = args.vcf_tsv
    filter_tags = args.filter_tags
    out = args.out
    csq = args.csq_headers
    
    # Initiate logging
    init_logging(f"{out}_logging.txt")
    
    # Load dataframes and merge
    vcfdf = parse_vcf(vcf_tsv, csq)
    gnomadf = update_allele_freq(parse_gNOMAD(gnomadcsv))
    merged_df = AF_to_numeric(merge(vcfdf, gnomadf))

    # Write out raw csv
    csv_out(merged_df, f"{out}_raw")
    
    # Plot raw
    plot(merged_df, "raw", out)
    
    # Filter if filter_tags are provided
    if filter_tags: 
        # Create new empty df
        total_filtered_df = merged_df.iloc[:0,:].copy()
        total_filtered_df["Filter_tag"] = "" 
        
        # Append filtered df to new df
        for tag in filter_tags:
            filtered_df = filter_tsv(merged_df.copy(), tag)
            
            # If the filtering criteria works, filtered dataframe is not empty
            # Otherwise, log to show dataframe filtering failed
            if filtered_df.empty:
                info(f"Filter tag {tag} is accepted by the programme but " +\
                     "filtering produced no output because filtering " +\
                     "criteria is not met.")
            else:
                info(f"Filter tag {tag} is successfully used")
                filtered_df.loc[:, "Filter_tag"] = tag
                total_filtered_df = total_filtered_df.append(filtered_df)
        
        # Remove identical rows but keep all filter tags
        remove = []
        total_filtered_df = total_filtered_df.sort_values(
            by = ["Chromosome", "Position", "Reference", "Alternate"],
            axis = 0).reset_index(drop = True)
        
        df_cols = len(total_filtered_df.columns)
        for i in total_filtered_df.index[:-1]:
            current_row = total_filtered_df.iloc[i, 0:df_cols-1]
            next_row = total_filtered_df.iloc[i + 1, 0:df_cols-1]
            if current_row.equals(next_row):
                remove += [i]
                total_filtered_df.loc[i + 1, "Filter_tag"] += \
                    f",{total_filtered_df.loc[i, 'Filter_tag']}"      
        total_filtered_df = total_filtered_df.drop(list(set(remove)))
                
        # Write out filtered csv
        plot(total_filtered_df, "filtered", out)
        csv_out(total_filtered_df, f"{out}_gnomadAnnotated.csv")
    
    info("Annotation complete")
    
# Run main() if this script is called from command-line    
if __name__ == "__main__":
    main()
