dir=~/Desktop/vc/resources

# GATK resources
mkdir -p $dir/genome
cd $dir/genome
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
mkdir -p $dir/vep_cache/ClinVar
cd $dir/vep_cache/ClinVar
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20221211.vcf.gz.md5
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20221211.vcf.gz.tbi

# VEP
mkdir -p $dir/vep_cache/
curl -o $dir/vep_cache/homo_sapiens_vep_108_GRCh38.tar.gz \
ftp://ftp.ensembl.org/pub/release-108/variation/indexed_vep_cache/homo_sapiens_vep_108_GRCh38.tar.gz
tar xzf $dir/vep_cache/homo_sapiens_vep_108_GRCh38.tar.gz

# Snakemake 
mamba create --name snakemake -y
mamba activate snakemake
mamba install -c bioconda snakemake -y
mamba deactivate

# Pangolin (optional)
mamba create --name pangolin pip python=3.8 -y
mamba activate pangolin
mamba install -c bioconda pyvcf -y
mamba install pytorch torchvision torchaudio pytorch-cuda -c pytorch -c nvidia -y
mamba install -c bioconda gffutils -y
mamba install -c conda-forge biopython -y
mamba install -c bioconda pyfastx -y
mamba install -c anaconda pandas -y
mamba install -c anaconda numpy -y
mamba deactivate

cd $dir
git clone https://github.com/tkzeng/Pangolin.git
cd $dir/resources/Pangolin
mamba activate pangolin
pip install .
mamba deactivate

curl -o $dir/Pangolin/gencode.v38.annotation.db \
https://www.dropbox.com/sh/6zo0aegoalvgd9f/AADOhGYJo8tbUhpscp3wSFj6a/gencode.v38.annotation.db

# GATK4 (optional)
mamba create --name gatk4 -y
mamba activate gatk4
mamba install -c bioconda gatk4 -y
mamba install -c anaconda pandas -y
mamba install -c anaconda numpy -y
mamba deactivate

# CI-SpliceAI (optional)
mamba create -n cis_use python=3.8 tensorflow-gpu cudnn keras -c conda-forge -y

mamba activate cis_use
mkdir -p $dir/resources/
git clone https://github.com/YStrauch/CI-SpliceAI__Annotation.git
cd $dir/resources/CI-SpliceAI__Annotation
pip install .
mamba deactivate
