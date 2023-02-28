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

## To run pipeline
Using 10 cores and max memory of 20 GB
```
snakemake --cores 10 --use-conda --resources mem_mb=20000 --retries 3
```
