# To run: snakemake --use-conda --retries 3 --cores 48 
# Caveats: 
#   only for grch38 (fasta and spliceai setting)
#   Nonmenclature: {sample}_R1_001.fastq.gz, {sample}_R2_001.fastq.gz

configfile: "config/config.yaml"

# Get prefix of fastq files
from os import listdir
from sys import exit
sample_names = []
for name in listdir("fastq"):
    if "_R1_001.fastq.gz" in name:
        sample_names.append(name[:-16])
    elif "_R2_001.fastq.gz" in name:
        sample_names.append(name[:-16])
    else:
        exit("Sample name is not properly formatted as " +\
        "{sample}_R1_001.fastq.gz, {sample}_R2_001.fastq.gz" +\
        "for read 1 and read 2, respectively")

rule all:
    input:
        "results/quality_control/multiqc_report.html"

rule fastp:
    input:
        read1="fastq/{sample}_R1_001.fastq.gz",
        read2="fastq/{sample}_R2_001.fastq.gz"
    output:
        read1=temp("fastq/clean_{sample}_R1_001.fastq.gz"),
        read2=temp("fastq/clean_{sample}_R2_001.fastq.gz"),
        report="results/quality_control/{sample}.fastp.html"
    conda:
        "envs/fastp.yaml"
    threads:
        config["threads"]
    shell:
        "fastp -i {input.read1} -I {input.read2} "
        "-o {output.read1} -O {output.read2} "
        "-h {output.report} --detect_adapter_for_pe "
        "--qualified_quality_phred 15 --unqualified_percent_limit 40 "
        "--length_required 15 --threads {threads} "

rule dragmap:
    input:
        read1="fastq/clean_{sample}_R1_001.fastq.gz",
        read2="fastq/clean_{sample}_R2_001.fastq.gz"
    output:
        temp("sam/{sample}.sam")
    conda:
        "envs/dragmap.yaml"
    shell:
        "dragen-os -r resources/reference "
        "-1 {input.read1} -2 {input.read2} "
        "--output-directory sam/ "
        "--output-file-prefix {wildcards.sample}"

rule sortsam:
    input:
        "sam/{sample}.sam"
    output:
        temp("bam/{sample}_raw.sorted.bam")
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk SortSam --I {input} --O {output} --SORT_ORDER coordinate"

rule markduplicates:
    input:
        "bam/{sample}_raw.sorted.bam"
    output:
        "bam/{sample}.sorted.removed_duplicates.bam"
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk MarkDuplicatesSpark --remove-all-duplicates I {input} -O {output}"

rule calibratedrag:
    input:
        genome="resources/genome/Homo_sapiens_assembly38_masked.fasta",
        bam="bam/{sample}.sorted.removed_duplicates.bam",
        str_table="resources/genome/Homo_sapiens_assembly38_masked.tsv"
    output:
        temp("{sample}.dragstr_model.txt")
    conda:
        "envs/gatk4.yaml"
    threads:
        config["threads"]
    shell:
        "gatk CalibrateDragstrModel -R {input.genome} -I {input.bam} "
        "-str {input.str_table} --threads {threads} -O {output}"

rule haplotypecaller:
    input:
        genome="resources/genome/Homo_sapiens_assembly38_masked.fasta",
        bam="bam/{sample}.sorted.removed_duplicates.bam",
        interval_list=config["intervals_list"]
    output:
        temp("vcf/{sample}_raw.vcf")
    conda:
        "envs/gatk4.yaml"
    threads:
        config["threads"]
    shell:
        "gatk HaplotypeCaller -R {input.genome} -I {input.bam} -O {output} "
        "--dragen-mode true --dragstr-params-path "
        "{wildcards.sample}.dragstr_model.txt "
        "--native-pair-hmm-threads {threads}"

rule hardfilter:
    input:
        "vcf/{sample}_raw.vcf"
    output:
        temp("vcf/{sample}.hardfiltered.vcf")
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk VariantFiltration -V {input} --filter-expression 'QUAL < 10.4139'"
        "--filter-name DRAGENHardQUAL -O {output}"

rule vep:
    input:
        genome="resources/genome/Homo_sapiens_assembly38_masked.fasta",
        vcf="vcf/{sample}.hardfiltered.vcf"
    output:
        temp("vcf/{sample}.vep.vcf")
    conda:
        "envs/vep.yaml"
    shell:
        "vep -i {input.vcf} o- {output} --cache --vcf "
        "--dir_cache /opt/vep/.vep/resources/VEP_population_database" 
        "--species 'homo_sapiens' --force_overwrite --af_1kg --af_gnomad "
        "--fasta {input.genome} --hgvs --offline "
        "--custom /opt/vep/.vep/resources/VEP_population_database/ClinVar/"
        "clinvar.vcf.gz,ClinVar,vcf,exact,0,ClinVar,vcf,exact,0,CLNDN,"
        "CLNDNINCL,CLNDISDB,CLNDISDBINCL,CLNHGVS,CLNREVSTAT,CLNSIG,CLNSIGCONF,"
        "CLNSIGINCL,CLNVC,CLNVCSO,CLNVI,DBVARID,GENEINFO,MC,ORIGIN,SSR"                     

rule spliceai:
    input:
        genome="resources/genome/Homo_sapiens_assembly38_masked.fasta",
        vcf="vcf/{sample}.vep.vcf"
    output:
        temp("vcf/{sample}.spliceai.vcf")
    conda:
        "envs/spliceai.yaml"
    shell:
        "spliceai -I {input.vcf} -O {output} -R {input.genome} -A grch38"

rule pangolin:
    input:
        genome="resources/genome/Homo_sapiens_assembly38_masked.fasta",
        vcf="vcf/{sample}.spliceai.vcf",
        annotation_db="resources/Pangolin/gencode.v38.annotation.db"
    output:
        "vcf/pangolin/{sample}.spliceai.pangolin.vcf"
    conda:
        "envs/pangolin.yaml"
    shell:
        "pangolin {input.vcf} {input.genome} {input.annotation_db} "
        "{wildcards.sample}"

rule variantstotable:
    input:
        "vcf/pangolin/{sample}.spliceai.pangolin.vcf"
    output:
        "vcf/tsv/{sample}.tsv"
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk VariantsToTable -V {input} -F CHROM -F POS -F ID -F REF -F ALT "
        "-F TYPE -F QUAL -F DP -F MQ -F MQRankSum -F QD -F ReadPosRankSum "
        "-F SOR -F GT -F GQ -F DP -F AD -F VAF -F PL -F FILTER -F CSQ -F CLNDN "
        "-F CLNDNINCL -F CLNDISDB -F CLNDISDBINCL -F CLNHGVS -F CLNREVSTAT "
        "-F CLNSIG -F CLNSIGCONF -F CLNSIGINCL -F CLNVC -F CLNVCSO -F CLNVI "
        "-F DBVARID -F GENEINFO -F MC -F ORIGIN -F SSR -F SpliceAI -F Pangolin "
        "-GF AD -O {output}"
        "touch vcf/tsv/all_variants.tsv"

rule multiqc:
    input:
        expand("vcf/tsv/{sample}.tsv", sample = sample_names)
    output:
        "results/quality_control/multiqc_report.html"
    conda:
        "envs/multiqc.yaml"
    shell:
        "multiqc -o {output} --force --ignore resources -o {output} ."

onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred")
    shell("mail -s 'an error occurred' {config['email']} < {log}")
