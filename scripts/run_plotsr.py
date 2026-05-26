import subprocess
import os
import multiprocessing

def run_minimap2_alignment(genome1 : str, genome2 : str, output_dir : str, minimap2_params : str):
    """
    Run minimap2 to align genome2 to genome1 and save the output in output_dir.
    """
    output_file = os.path.join(output_dir, f"alignment_{os.path.basename(genome1)}_vs_{os.path.basename(genome2)}")
    command = f"minimap2 -ax asm5 -t 1 --eqx {genome1} {genome2} > {output_file}.sam | samtools sort -O BAM - > {output_file}.bam"
    subprocess.run(command, shell=True, check=True)

    # index the bam file
    command = f"samtools index {output_file}.bam"
    subprocess.run(command, shell=True, check=True)

    # call structural variant calling script here
    command = f"syri -c {output_file}.bam -r {genome1} -q {genome2} -F B --prefix {output_file}"
    subprocess.run(command, shell=True, check=True)

    return output_file_prefix

def run_progressive_minimap2(reference : str, fasta_dir : str, output_dir : str, minimap2_params : str = "-ax asm5", threads : int = 1, FASTA_EXTENSIONS = (".fasta", ".fa", ".fna")):
    """
    Run minimap2 to align the sequences in the fasta_dir to the reference genome and save the output in output_dir.
    """

    # generate ordered list of fasta files in the fasta_dir
    file_list = []
    for file in os.listdir(fasta_dir):
        if file.endswith(tuple(FASTA_EXTENSIONS)):
            if file != reference:  # Exclude the reference genome from the file list
                file_list.append(os.path.join(fasta_dir, file))
    
    file_list.sort()  # Sort the file list to ensure consistent order

    # add reference genome to start of the file list
    file_list.insert(0, reference)

    # generate all iterative pairs of genomes for alignment
    pairs = []
    for i in range(len(file_list) - 1):
        pairs.append((file_list[i], file_list[i+1]))

    # parrallel processing of minimap2 alignments
    with multiprocessing.Pool(processes=threads) as pool:
        output_file_prefix = pool.starmap(run_minimap2_alignment, [(pair[0], pair[1], output_dir, minimap2_params) for pair in pairs])
        
    return output_file_prefix, file_list

def run_plot_synteny(reference : str, fasta_dir : str, output_dir : str, minimap2_params : str = "-ax asm5", threads : int = 1, FASTA_EXTENSIONS = (".fasta", ".fa", ".fna")):
    
    # create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # generate alignments
    output_file_prefix, file_list = run_progressive_minimap2(reference, fasta_dir, output_dir, minimap2_params, threads, FASTA_EXTENSIONS)

    # generate list of genome names for plotsr
    genome_names = []
    for file in file_list:
        genome_names.append((file, os.path.basename(file)))
    
    with open(os.path.join(output_dir, "genome_names.txt"), "w") as f:
        for name in genome_names:
            f.write(f"{name[0]}\t{name[1]}\n")

    # run plotsr
    plotsr_file_commands = "--sr " + " --sr ".join([f"{output_file_prefix[i]}_syri.out" for i in range(len(output_file_prefix))])
    command = f"plotsr {plotsr_file_commands} --genome-names {os.path.join(output_dir, 'genome_names.txt')} -o {os.path.join(output_dir, 'synteny_plot.pdf')}"
    subprocess.run(command, shell=True, check=True)

run_plot_synteny(snakemake.input.reference, snakemake.input.fasta_dir, snakemake.params.output_dir, snakemake.params.minimap2_params, snakemake.threads, snakemake.params.FASTA_EXTENSIONS)
