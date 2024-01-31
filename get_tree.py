#!/usr/bin/env python
import argparse
import os
import shutil
import subprocess
import sys
import tempfile

from Bio import SeqIO
from Bio import Phylo


# Reading the FASTA file and filtering sequences
def filter_sequences(input_fasta, names_to_keep):
    for record in SeqIO.parse(input_fasta, "fasta"):
        if record.id in names_to_keep:
            yield record


# Function to run MAFFT
def run_mafft(input_fasta: str):
    # MAFFT command
    mafft_cmd = [f"{os.path.dirname(sys.executable)}/mafft", "--auto", input_fasta]
    return subprocess.run(mafft_cmd, capture_output=True, text=True)


def run_trimal(input_aln: str, output_phy: str):
    # trimal command
    trimal_cmd = [
        f"{os.path.dirname(sys.executable)}/trimal",
        "-in",
        input_aln,
        "-phylip",
        "-out",
        output_phy,
    ]
    return subprocess.run(trimal_cmd, capture_output=True, text=True)


# Function to run RAxML
def run_raxml(alignment_file: str, output_dir: str):
    # RAxML command
    raxml_cmd = [
        f"{os.path.dirname(sys.executable)}/raxmlHPC-PTHREADS-SSE3",
        "-s",
        alignment_file,
        "-m",
        "GTRGAMMA",
        "-f",
        "a",
        "-p",
        "42",
        "-x",
        "60",
        "-N",
        "100",
        "-n",
        "whatever",
        "-T",
        str(16),
        "-w",
        output_dir,
    ]
    return subprocess.run(raxml_cmd, capture_output=True, text=True, cwd=".")


def process_sequences(
    zotu_fasta_path: str,
    alignment_out: str,
    tree_out: str,
    seq_names: list[str] = None,
) -> int:
    if seq_names is None:
        input_fasta = zotu_fasta_path
    else:
        # with tempfile.NamedTemporaryFile(
        #     mode="w", delete=False, suffix=".fasta"
        # ) as temp_fasta:
        temp_fasta = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta")
        # Filter and write sequences to a temp file
        SeqIO.write(
            filter_sequences(zotu_fasta_path, seq_names),
            temp_fasta,
            "fasta",
        )
        input_fasta = temp_fasta.name

    # Run MAFFT on the temp file and check if it ran successfully
    mafft_result = run_mafft(input_fasta)
    if mafft_result.returncode != 0:
        print("Error in running MAFFT:", mafft_result.stderr)
        os.remove(temp_fasta.name)
        return 1
    temp_aln = tempfile.NamedTemporaryFile(mode="w", suffix=".aln")
    temp_aln.write(mafft_result.stdout)

    # Run trimal on the alignment file and check if it ran successfully
    trimal_result = run_trimal(temp_aln.name, alignment_out)
    temp_aln.close()
    if trimal_result.returncode != 0:
        print("Error in running trimal:", trimal_result.stderr)
        return 1
    else:
        print(
            "MAFFT and trimal completed successfully. Alignment saved as:",
            alignment_out,
        )

    # Run RAxML, output to a temporary directory, copy the tree file to the output path
    with tempfile.TemporaryDirectory() as temp_dir:
        raxml_result = run_raxml(alignment_out, temp_dir)
        tree_file = os.path.join(temp_dir, "RAxML_bestTree.whatever")
        shutil.copy(tree_file, tree_out)

    if seq_names is not None:
        temp_fasta.close()

    # Check if RAxML ran successfully
    if raxml_result.returncode != 0:
        print("Error in running RAxML:", raxml_result.stderr)
        return 1
    else:
        print("RAxML completed successfully. Tree saved in Newick format as:", tree_out)

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a phylogenetic tree from a FASTA file of sequences and a list of sequence names"
    )
    parser.add_argument(
        "-i",
        "--input_fasta",
        help="Input FASTA file",
        required=True,
        type=str,
        metavar="FASTA",
    )
    args = parser.parse_args()

    # Read the sequence names file

    # Process the sequences
    # prefix is dirname + anything before first dot
    prefix = os.path.join(
        os.path.dirname(args.input_fasta),
        os.path.basename(args.input_fasta).split(".")[0],
    )
    aln_file = prefix + ".phy"
    tree_file = prefix + ".newick"
    process_sequences(args.input_fasta, aln_file, tree_file)

    # Read the tree file and print the tree
    # tree = Phylo.read(args.output, "newick")
    # Phylo.draw_ascii(tree)
