from pathlib import Path

configfile: "config/config.yaml"

REFERENCE = config["reference"]
GENOME_DIR = config["input_dir"]
FASTA_EXTENSIONS = config.get("fasta_extensions", [".fasta", ".fa", ".fna"])
MINIMAP2_INDEX = REFERENCE + ".mmi"
PEGAS_SCRIPT = "scripts/pegas_haplotype_analysis.R"


def resolve_sample_genome(sample):
    for extension in FASTA_EXTENSIONS:
        candidate = Path(GENOME_DIR) / f"{sample}{extension}"
        if candidate.exists():
            return str(candidate)
    raise ValueError(
        f"Could not find genome FASTA for sample '{sample}' in {GENOME_DIR} "
        f"with extensions: {', '.join(FASTA_EXTENSIONS)}. "
        f"Expected format: {{sample}}<extension> (e.g., {sample}.fasta)."
    )


def discover_samples():
    discovered = []
    for extension in FASTA_EXTENSIONS:
        pattern = str(Path(GENOME_DIR) / "{sample}") + extension
        discovered.extend(glob_wildcards(pattern).sample)
    return sorted(set(discovered))


if config.get("samples"):
    SAMPLES = config["samples"]
else:
    SAMPLES = discover_samples()

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
        MINIMAP2_INDEX
    conda:
        "envs/alignment.yaml"
    shell:
        "minimap2 -d {output} {input}"


rule align_sample:
    input:
        ref_index=MINIMAP2_INDEX,
        genome=lambda wc: resolve_sample_genome(wc.sample),
    output:
        bam=temp("results/alignment/{sample}.sorted.bam")
    threads: 40
    params:
        preset=config.get("minimap2_preset", "asm5"),
    log:
        "results/logs/alignment/{sample}.log"
    conda:
        "envs/alignment.yaml"
    shell:
        (
            "mkdir -p results/alignment results/logs/alignment && "
            "minimap2 -a -x {params.preset} -t {threads} {input.ref_index} {input.genome} 2> {log} | "
            "samtools sort -@ {threads} -o {output.bam}"
        )


rule index_bam:
    input:
        bam="results/alignment/{sample}.sorted.bam"
    output:
        bai=temp("results/alignment/{sample}.sorted.bam.bai")
    conda:
        "envs/alignment.yaml"
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
    threads: 40
    conda:
        "envs/variants.yaml"
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
    conda:
        "envs/variants.yaml"
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
    conda:
        "envs/plink.yaml"
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
    conda:
        "envs/pegas.yaml"
    shell:
        "mkdir -p results/pegas && Rscript {params.script} {input.vcf} {output}"
