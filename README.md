# Snakemake Gene Panel Variant Calling Pipeline

## Set up Snakemake
```mamba create --name snakemake -c bioconda -y snakemake```

## Set up Snakemake directory
```
# Clone Git repository into vc directory
git clone https://github.com/CherWeiYuan/GenePanelVariantCallingSmk.git

# Download GATK resources
cd GenePanelVariantCallingSmk/resources/genome
gsutil -m cp \
  "gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz" \
  "gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
  "gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi" \
  "gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz" \
  "gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi" \
  .

# ClinVar
cd GenePanelVariantCallingSmk/resources/vep_cache/ClinVar
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.md5

# VEP
curl -o GenePanelVariantCallingSmk/resources/vep_cache/homo_sapiens_vep_108_GRCh38.tar.gz \
ftp://ftp.ensembl.org/pub/release-108/variation/indexed_vep_cache/homo_sapiens_vep_108_GRCh38.tar.gz
tar xzf GenePanelVariantCallingSmk/resources/vep_cache/homo_sapiens_vep_108_GRCh38.tar.gz
```

## Run pipeline
Important notes:
* Input are paired-end fastq files in the file name format: 
  - Read 1 as {sample_name}_R1_001.fastq.gz
  - Read 2 as {sample_name}_R2_001.fastq.gz
* Reads will be mapped to GRCh38 reference genome

Using 4 cores and max memory of 10 GB
```
snakemake --cores 4 --use-conda --resources mem_mb=10000 --retries 3
```

## Output
Output will be sent to "results" folder in the working directory.
There are two output folders: "all_data" and "simplified_data"

Generally, the CSV reports in "simplified_data" will suffice:\

| **Header**        | **Description**     |
|:------------:|:---------------------------------------|
| ```HGVSg```       | HGVS genomic      |
|```HGVSc```        | HGVS coding       |
|```HGVSp```        | HGVS protein      |
|```SYMBOL```       | Gene name         |
|```Gene```         | ENSEMBL gene ID   |
|```Existing_variation```| Variant ID, most frequently rsID |
|```Feature``` | ENSEMBL feature ID     |
|```Feature_type``` | Feature type, such as Transcript or Regulatory Feature|
|```Consequence```	| VEP calculated consequence, such as frameshift_variant, intron_variant, start_loss, missense_variant, intron_variant or stop_gained |
|```AF``` | Global allele frequency|
|```EAS_AF``` | East Asian allele frequency |
|```gnomADe_AF```  | gNOMAD v2.1 exome (WES) population database allele frequency   |
|```gnomADe_EAS_AF```   | gNOMAD v2.1 exome (WES) population database allele frequency for East Asians  |
|```gnomADg_AF``` | gNOMAD v3.2.1 genome (WGS) population database allele frequency |
|```gnomADg_EAS_AF```   | gNOMAD v3.2.1 genome (WGS) population database allele frequency for East Asians    |
|```ClinVar_CLNSIG```	  | ClinVar clinical significance, such as Pathogenic, Likely_pathogenic or Conflicting_interpretations_of_pathogenicity   |
|```SpliceAI_SpliceAI_highest_score``` | SpliceAI score, define as the highest delta score among the four categories: acceptor gain/ loss, donor gain/loss|

Should more information be required, "all_data" folder will provide a comprehensive report for each sample
