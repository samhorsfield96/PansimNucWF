import subprocess
import os
import multiprocessing
import random

def normalize_fasta_chromosomes(input_fasta: str, output_dir: str) -> str:
    """
    Rename all chromosomes in a FASTA file to chr1, chr2, etc. to ensure
    consistent chromosome IDs across all genomes for syri/plotsr.
    Returns path to the normalized FASTA.
    """
    basename = os.path.basename(input_fasta)
    output_fasta = os.path.join(output_dir, basename)
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        chrom_count = 0
        for line in infile:
            if line.startswith('>'):
                chrom_count += 1
                outfile.write(f'>chr{chrom_count}\n')
            else:
                outfile.write(line)
    return output_fasta

def run_minimap2_alignment(genome1 : str, genome2 : str, output_dir : str, minimap2_params : str, log_file : str) -> str:
    """
    Run minimap2 to align genome2 to genome1 and save the output in output_dir.
    """

    file_prefix = os.path.splitext(os.path.basename(genome1))[0] + "_vs_" + os.path.splitext(os.path.basename(genome2))[0]
    output_file_prefix = os.path.join(output_dir, file_prefix)
    
    command = f"minimap2 {minimap2_params} -t 1 --eqx {genome1} {genome2} 2>>{log_file} | samtools sort -O BAM - > {output_file_prefix}.bam 2>>{log_file}"
    subprocess.run(command, shell=True, check=True)

    # index the bam file
    command = f"samtools index {output_file_prefix}.bam >> {log_file} 2>&1"
    subprocess.run(command, shell=True, check=True)

    # call structural variant calling script here
    command = f"syri -c {output_file_prefix}.bam -r {genome1} -q {genome2} -F B --dir {output_dir} --prefix {file_prefix}_ >> {log_file} 2>&1"
    subprocess.run(command, shell=True, check=True)

    return output_file_prefix

def run_progressive_minimap2(reference : str, fasta_dir : str, output_dir : str, max_alignments : int, log_file : str, minimap2_params : str = "-ax asm5", threads : int = 1, FASTA_EXTENSIONS = (".fasta", ".fa", ".fna")):
    """
    Run minimap2 to align the sequences in the fasta_dir to the reference genome and save the output in output_dir.
    """

    # generate ordered list of fasta files in the fasta_dir
    file_list = []
    for file in os.listdir(fasta_dir):
        if file.endswith(tuple(FASTA_EXTENSIONS)):
            if file != os.path.basename(reference):  # Exclude the reference genome from the file list
                file_list.append(os.path.join(fasta_dir, file))
                
    # randomly sample to get max_alignments number of files if there are more than max_alignments
    if len(file_list) > max_alignments:
        file_list = random.sample(file_list, max_alignments)

    file_list.sort()  # Sort the file list to ensure consistent order

    # add reference genome to start of the file list
    file_list.insert(0, reference)

    # Normalize chromosome IDs across all FASTAs so syri/plotsr use consistent names
    norm_dir = os.path.join(output_dir, "normalized_fasta")
    os.makedirs(norm_dir, exist_ok=True)
    file_list = [normalize_fasta_chromosomes(f, norm_dir) for f in file_list]

    # generate all iterative pairs of genomes for alignment
    pairs = []
    for i in range(len(file_list) - 1):
        pairs.append((file_list[i], file_list[i+1]))

    # parrallel processing of minimap2 alignments
    output_file_prefixes = []
    with multiprocessing.Pool(processes=threads) as pool:
        for output_file_prefix in pool.starmap(run_minimap2_alignment, [(pair[0], pair[1], output_dir, minimap2_params, log_file) for pair in pairs]):
            output_file_prefixes.append(output_file_prefix)
        
    return output_file_prefixes, file_list

def run_plot_synteny(reference : str, fasta_dir : str, output_dir : str, max_alignments : int, log_file : str, minimap2_params : str = "-ax asm5", threads : int = 1, FASTA_EXTENSIONS = (".fasta", ".fa", ".fna")):
    
    # create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # generate alignments
    output_file_prefixes, file_list = run_progressive_minimap2(reference, fasta_dir, output_dir, max_alignments, log_file, minimap2_params, threads, FASTA_EXTENSIONS)

    # generate list of genome names for plotsr
    genome_names = []
    for file in file_list:
        genome_names.append((file, os.path.splitext(os.path.basename(file))[0]))
    
    with open(os.path.join(output_dir, "genome_names.txt"), "w") as f:
        f.write("#file\tname\n")
        for name in genome_names:
            f.write(f"{name[0]}\t{name[1]}\n")

    # run plotsr
    plotsr_file_commands = "--sr " + " --sr ".join([f"{output_file_prefixes[i]}_syri.out" for i in range(len(output_file_prefixes))])
    command = f"plotsr {plotsr_file_commands} --genomes {os.path.join(output_dir, 'genome_names.txt')} -o {os.path.join(output_dir, 'synteny_plot.pdf')} >> {log_file} 2>&1"
    subprocess.run(command, shell=True, check=True)

run_plot_synteny(snakemake.input.reference, snakemake.params.fasta_dir, snakemake.params.output_dir, snakemake.params.max_alignments, snakemake.log[0], snakemake.params.minimap2_params, snakemake.threads, snakemake.params.FASTA_EXTENSIONS)
