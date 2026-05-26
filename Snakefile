from pathlib import Path

configfile: "config/config.yaml"

REFERENCE = config["reference"]
GENOME_DIR = config["input_dir"]
FASTA_EXTENSIONS = config.get("fasta_extensions", [".fasta", ".fa", ".fna"])
MINIMAP2_INDEX = REFERENCE + ".mmi"
REFERENCE_FAI  = REFERENCE + ".fai"
PEGAS_SCRIPT = "scripts/pegas_haplotype_analysis.R"
RECOMBINATION_THRESHOLD = config.get("pegas", {}).get("recombination_threshold", 0.9)
PLINK_PLOT_SCRIPT = "scripts/plink_ld_plots.R"
SFS_NUC_SCRIPT = "scripts/plot_sfs_nuc.R"
OUTPUT_DIR = config["output_dir"]


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
    ref_resolved = str(Path(REFERENCE).resolve())
    discovered = []
    for extension in FASTA_EXTENSIONS:
        pattern = str(Path(GENOME_DIR) / "{sample}") + extension
        for s in glob_wildcards(pattern).sample:
            candidate = str((Path(GENOME_DIR) / f"{s}{extension}").resolve())
            if candidate != ref_resolved:
                discovered.append(s)
    return sorted(set(discovered))


SAMPLES = discover_samples()

if not SAMPLES:
    raise ValueError(
        "No samples found. Add sample IDs to config['samples'] or provide FASTA files in input_dir."
    )


rule all:
    input:
        f"{OUTPUT_DIR}/variants/filtered_variants.vcf.gz",
        f"{OUTPUT_DIR}/variants/filtered_variants.vcf.gz.tbi",
        f"{OUTPUT_DIR}/plink/ld_decay.ld",
        f"{OUTPUT_DIR}/plink/ld_heatmap.pdf",
        f"{OUTPUT_DIR}/plink/ld_decay_plot.pdf",
        f"{OUTPUT_DIR}/pegas/haplotype_summary.tsv",
        f"{OUTPUT_DIR}/pegas/haplotype_network.pdf",
        f"{OUTPUT_DIR}/sfs/sfs_nuc_density_minor_alleles.pdf",
        f"{OUTPUT_DIR}/sfs/sfs_nuc_density_all_alleles.pdf",
        f"{OUTPUT_DIR}/sfs/sfs_nuc_sfs.csv",


rule index_reference:
    input:
        REFERENCE
    output:
        MINIMAP2_INDEX
    conda:
        "envs/alignment.yaml"
    shell:
        "minimap2 -d {output} {input}"


rule faidx_reference:
    input:
        REFERENCE
    output:
        REFERENCE_FAI
    conda:
        "envs/alignment.yaml"
    shell:
        "samtools faidx {input}"


rule align_sample:
    input:
        ref_index=MINIMAP2_INDEX,
        ref=REFERENCE,
        genome=lambda wc: resolve_sample_genome(wc.sample),
    output:
        paf=f"{OUTPUT_DIR}/alignment/{{sample}}.paf.gz",
    threads: 40
    params:
        preset=config.get("minimap2_preset", "asm5"),
    log:
        f"{OUTPUT_DIR}/logs/alignment/{{sample}}.log"
    conda:
        "envs/alignment.yaml"
    shell:
        (
            f"mkdir -p {OUTPUT_DIR}/alignment {OUTPUT_DIR}/logs/alignment && "
            "minimap2 -c --cs -x {params.preset} -t {threads} {input.ref_index} {input.genome} 2> {log} | "
            "sort -k6,6 -k8,8n | bgzip > {output.paf}"
        )


rule call_variants_paftools:
    input:
        paf=f"{OUTPUT_DIR}/alignment/{{sample}}.paf.gz",
        ref=REFERENCE,
        fai=REFERENCE_FAI,
    output:
        vcf=f"{OUTPUT_DIR}/variants/per_sample/{{sample}}.vcf.gz",
        tbi=f"{OUTPUT_DIR}/variants/per_sample/{{sample}}.vcf.gz.tbi",
    conda:
        "envs/alignment.yaml"
    shell:
        (
            f"mkdir -p {OUTPUT_DIR}/variants/per_sample && "
            "zcat {input.paf} | paftools.js call -l 50 -L 50 -s {wildcards.sample} -f {input.ref} - | "
            "bgzip > {output.vcf} && "
            "tabix -p vcf {output.vcf}"
        )


rule merge_vcfs:
    input:
        vcfs=expand(f"{OUTPUT_DIR}/variants/per_sample/{{sample}}.vcf.gz", sample=SAMPLES),
        tbis=expand(f"{OUTPUT_DIR}/variants/per_sample/{{sample}}.vcf.gz.tbi", sample=SAMPLES),
    output:
        vcf=f"{OUTPUT_DIR}/variants/raw_variants.vcf.gz",
        tbi=f"{OUTPUT_DIR}/variants/raw_variants.vcf.gz.tbi",
    conda:
        "envs/variants.yaml"
    shell:
        (
            f"mkdir -p {OUTPUT_DIR}/variants && "
            "bcftools merge --force-samples -Oz -o {output.vcf} {input.vcfs} && "
            "tabix -p vcf {output.vcf}"
        )


rule filter_variants:
    input:
        vcf=f"{OUTPUT_DIR}/variants/raw_variants.vcf.gz",
        tbi=f"{OUTPUT_DIR}/variants/raw_variants.vcf.gz.tbi",
    output:
        vcf=f"{OUTPUT_DIR}/variants/filtered_variants.vcf.gz",
        tbi=f"{OUTPUT_DIR}/variants/filtered_variants.vcf.gz.tbi",
    params:
        min_quality=config.get("variant_calling", {}).get("min_quality", 30),
    conda:
        "envs/variants.yaml"
    shell:
        (
            f"bcftools filter -i 'QUAL>={params.min_quality}' -Oz -o {output.vcf} {input.vcf} && "
            f"tabix -p vcf {output.vcf}"
        )


rule fill_ref_genotypes:
    input:
        vcf=f"{OUTPUT_DIR}/variants/filtered_variants.vcf.gz",
        tbi=f"{OUTPUT_DIR}/variants/filtered_variants.vcf.gz.tbi",
    output:
        vcf=f"{OUTPUT_DIR}/variants/filtered_variants_refgt.vcf.gz",
        tbi=f"{OUTPUT_DIR}/variants/filtered_variants_refgt.vcf.gz.tbi",
    conda:
        "envs/variants.yaml"
    shell:
        (
            "bcftools view -e 'COUNT(GT=\"mis\")=N_SAMPLES' {input.vcf} | "
            "bcftools +setGT -- -t . -n 0 | "
            "bgzip > {output.vcf} && "
            "tabix -p vcf {output.vcf}"
        )


rule plink_ld_decay:
    input:
        vcf=f"{OUTPUT_DIR}/variants/filtered_variants_refgt.vcf.gz",
        tbi=f"{OUTPUT_DIR}/variants/filtered_variants_refgt.vcf.gz.tbi",
    output:
        f"{OUTPUT_DIR}/plink/ld_decay.ld"
    threads: 40
    params:
        out_prefix=f"{OUTPUT_DIR}/plink/ld_decay",
        ld_window=config.get("plink", {}).get("ld_window", 99999),
        ld_window_kb=config.get("plink", {}).get("ld_window_kb", 1000),
    conda:
        "envs/plink.yaml"
    shell:
        (
            f"mkdir -p {OUTPUT_DIR}/plink && "
            "plink --threads {threads} --vcf {input.vcf} --double-id --allow-extra-chr --memory 8000 "
            "--mac 1 --r2 --ld-window {params.ld_window} --ld-window-kb {params.ld_window_kb} "
            "--ld-window-r2 0 --out {params.out_prefix}"
        )


rule plink_ld_plots:
    input:
        ld=f"{OUTPUT_DIR}/plink/ld_decay.ld",
        fai=REFERENCE_FAI,
    output:
        heatmap=f"{OUTPUT_DIR}/plink/ld_heatmap.pdf",
        decay=f"{OUTPUT_DIR}/plink/ld_decay_plot.pdf",
    params:
        script=PLINK_PLOT_SCRIPT,
        out_dir=f"{OUTPUT_DIR}/plink",
    conda:
        "envs/plink.yaml"
    shell:
        "Rscript {params.script} {input.ld} {params.out_dir} {input.fai}"


rule pegas_haplotype_analysis:
    input:
        vcf=f"{OUTPUT_DIR}/variants/filtered_variants_refgt.vcf.gz"
    output:
        tsv=f"{OUTPUT_DIR}/pegas/haplotype_summary.tsv",
        pdf=f"{OUTPUT_DIR}/pegas/haplotype_network.pdf",
    params:
        script=PEGAS_SCRIPT,
        recombination_threshold=RECOMBINATION_THRESHOLD
    conda:
        "envs/pegas.yaml"
    shell:
        f"mkdir -p {OUTPUT_DIR}/pegas && Rscript {{params.script}} {{input.vcf}} {{output.tsv}} {{output.pdf}} {params.recombination_threshold}"

rule sfs_nuc_plot:
    input:
        vcf=f"{OUTPUT_DIR}/variants/filtered_variants_refgt.vcf.gz",
        tbi=f"{OUTPUT_DIR}/variants/filtered_variants_refgt.vcf.gz.tbi",
    output:
        minor_pdf=f"{OUTPUT_DIR}/sfs/sfs_nuc_density_minor_alleles.pdf",
        both_pdf=f"{OUTPUT_DIR}/sfs/sfs_nuc_density_all_alleles.pdf",
        csv=f"{OUTPUT_DIR}/sfs/sfs_nuc_sfs.csv",
    params:
        script=SFS_NUC_SCRIPT,
        out_prefix=f"{OUTPUT_DIR}/sfs/sfs_nuc",
        gff=config.get("reference_gff", ""),
    conda:
        "envs/sfs_nuc.yaml"
    shell:
        f"mkdir -p {OUTPUT_DIR}/sfs && "
        "if [ -n \"{params.gff}\" ]; then "
        "Rscript {params.script} {input.vcf} {params.gff} {params.out_prefix}; "
        "else "
        "Rscript {params.script} {input.vcf} {params.out_prefix}; "
        "fi"
