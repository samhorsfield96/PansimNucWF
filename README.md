# PansimNucWF

Workflows for analysis of PansimNuc simulated genomes.

## Snakemake workflow

This repository contains a Snakemake pipeline that:

1. Aligns many individual genome FASTA files to a reference (`minimap2`)
2. Calls joint variants across samples (`samtools`)
3. Produces filtered VCF output for downstream analysis
4. Runs LD-decay outputs with `PLINK`
5. Runs haplotype summaries with `pegas` (R)

### Expected inputs

- Reference genome FASTA: `resources/reference/genome.fasta` (default; configurable)
- Individual genome FASTA files named `{sample}.fasta`, `{sample}.fa`, or `{sample}.fna` under `resources/genomes/` (default; configurable)

Configure paths and parameters in `config/config.yaml`.
You can customize accepted sample FASTA suffixes with `fasta_extensions`.

### Run

```bash
snakemake --cores 8 --use-conda
```

Or dry-run:

```bash
snakemake -n --use-conda
```

### Main outputs

- `results/variants/filtered_variants.vcf.gz`
- `results/plink/ld_decay.ld`
- `results/pegas/haplotype_summary.tsv`
