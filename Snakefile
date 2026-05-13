from pathlib import Path

configfile: "config/config.yaml"

REFERENCE = config["reference"]
GENOME_DIR = config["input_dir"]
BWA_INDEX_SUFFIXES = ["amb", "ann", "bwt", "pac", "sa"]
PEGAS_SCRIPT = "scripts/pegas_haplotype_analysis.R"

if config.get("samples"):
    SAMPLES = config["samples"]
else:
    SAMPLES = glob_wildcards(str(Path(GENOME_DIR) / "{sample}.fasta")).sample

if not SAMPLES:
    raise ValueError(
        "No samples found. Add sample IDs to config['samples'] or provide FASTA files in input_dir."
    )


rule all:
    input:
        "results/variants/filtered_variants.vcf.gz",
        "results/variants/filtered_variants.vcf.gz.tbi",
        "results/plink/ld_decay.ld",
        "results/pegas/haplotype_summary.tsv",


rule index_reference:
    input:
        REFERENCE
    output:
        expand("{ref}.{suffix}", ref=REFERENCE, suffix=BWA_INDEX_SUFFIXES)
    shell:
        "bwa index {input}"


rule align_sample:
    input:
        ref=REFERENCE,
        ref_index=rules.index_reference.output,
        genome=lambda wc: str(Path(GENOME_DIR) / f"{wc.sample}.fasta"),
    output:
        bam=temp("results/alignment/{sample}.sorted.bam")
    threads:
        config.get("threads", 4)
    log:
        "results/logs/alignment/{sample}.log"
    shell:
        (
            "mkdir -p results/alignment results/logs/alignment && "
            "bwa mem -t {threads} {input.ref} {input.genome} 2> {log} | "
            "samtools sort -@ {threads} -o {output.bam}"
        )


rule index_bam:
    input:
        bam="results/alignment/{sample}.sorted.bam"
    output:
        bai="results/alignment/{sample}.sorted.bam.bai"
    shell:
        "samtools index {input.bam}"


rule call_variants:
    input:
        ref=REFERENCE,
        bams=expand("results/alignment/{sample}.sorted.bam", sample=SAMPLES),
        bais=expand("results/alignment/{sample}.sorted.bam.bai", sample=SAMPLES),
    output:
        vcf="results/variants/raw_variants.vcf.gz",
        tbi="results/variants/raw_variants.vcf.gz.tbi",
    params:
        min_mapping_quality=config.get("variant_calling", {}).get("min_mapping_quality", 20),
    threads:
        config.get("threads", 4)
    shell:
        (
            "mkdir -p results/variants && "
            "bcftools mpileup -Ou -q {params.min_mapping_quality} -f {input.ref} {input.bams} | "
            "bcftools call -mv -Oz -o {output.vcf} && "
            "tabix -p vcf {output.vcf}"
        )


rule filter_variants:
    input:
        vcf="results/variants/raw_variants.vcf.gz",
        tbi="results/variants/raw_variants.vcf.gz.tbi",
    output:
        vcf="results/variants/filtered_variants.vcf.gz",
        tbi="results/variants/filtered_variants.vcf.gz.tbi",
    params:
        min_quality=config.get("variant_calling", {}).get("min_quality", 30),
    shell:
        (
            "bcftools filter -i 'QUAL>={params.min_quality}' -Oz -o {output.vcf} {input.vcf} && "
            "tabix -p vcf {output.vcf}"
        )


rule plink_ld_decay:
    input:
        vcf="results/variants/filtered_variants.vcf.gz",
        tbi="results/variants/filtered_variants.vcf.gz.tbi",
    output:
        "results/plink/ld_decay.ld"
    params:
        out_prefix="results/plink/ld_decay",
        ld_window=config.get("plink", {}).get("ld_window", 99999),
        ld_window_kb=config.get("plink", {}).get("ld_window_kb", 1000),
    shell:
        (
            "mkdir -p results/plink && "
            "plink --vcf {input.vcf} --double-id --allow-extra-chr "
            "--r2 --ld-window {params.ld_window} --ld-window-kb {params.ld_window_kb} "
            "--ld-window-r2 0 --out {params.out_prefix}"
        )


rule pegas_haplotype_analysis:
    input:
        vcf="results/variants/filtered_variants.vcf.gz"
    output:
        "results/pegas/haplotype_summary.tsv"
    params:
        script=PEGAS_SCRIPT
    shell:
        "mkdir -p results/pegas && Rscript {params.script} {input.vcf} {output}"
