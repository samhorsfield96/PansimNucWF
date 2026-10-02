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

### Installation

Install the required packages using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)/[mamba](https://github.com/mamba-org/mamba):

```
git clone https://github.com/samhorsfield96/PansimNucWF.git
cd PansimNucWF
mamba env create -n PansimNucWF "snakemake>=9.19.0"
mamba activate PansimNucWF
```

### Run

```bash
snakemake --cores 8 --use-conda 
```

If using mamba, use

```bash
snakemake --cores 8 --use-conda --conda-frontend mamba
```

Or dry-run:

```bash
snakemake -n --use-conda
```

### Main outputs

The paths below are under the configured `output_dir` (shown as `results/`).

- Variants: `results/variants/filtered_variants.vcf.gz` and its `.tbi` index
- LD analysis: `results/plink/ld_decay.ld`, `ld_heatmap.pdf`, `ld_decay_plot_mean.pdf`, and `ld_decay_plot_per_dist.pdf`
- Site-frequency spectrum: `results/sfs/sfs_nuc_density_minor_alleles.pdf`, `sfs_nuc_density_all_alleles.pdf`, and `sfs_nuc_sfs.csv`

Haplotype outputs depend on `simulated`:

- With `simulated: false`: `results/pegas/haplotype_summary.tsv` and `haplotype_network.pdf`
- With `simulated: true`: `results/haplotypes/haplotypes_haplotype_summary.csv`, `haplotypes_haplotype_freq.pdf`, `haplotypes_haplotype_composition.pdf`, `haplotypes_per_haplotype_composition.pdf`, and `haplotypes_sel_coeff_composition.pdf`; plus `results/gene_freq/gene_frequencies_frequency_dist.pdf`, `gene_frequencies_distribution.csv`, and `gene_frequencies_per_gene.csv`

Additional outputs are included when their corresponding conditions are met: `results/sv/sv_plot.pdf` (`simulated: true`, `plot_SVs: true`); `results/synteny/synteny_plot.pdf` (`simulated: false`, `plot_SVs: true`); gene-distance plots and CSVs (`simulated: true`, `plot_gene_dists: true`); TE copy-number CSVs (`simulated: true`, `plot_TEs: true`); and DFE plots when `simulated: true` and `resources/genomes/selection_samples.csv` exists.
