#!/usr/bin/env python
import argparse
import os
import shutil
import subprocess
import tempfile
from typing import Iterable

from Bio import SeqIO
from Bio import SeqRecord


# Reading the FASTA file and filtering sequences
def filter_sequences(
    input_fasta: str, names_to_keep: list[str]
) -> Iterable[SeqRecord.SeqRecord]:
    for record in SeqIO.parse(input_fasta, "fasta"):
        if record.id in names_to_keep:
            yield record


# Function to run MAFFT
def run_mafft(input_fasta: str):
    # MAFFT command
    mafft_cmd = ["mafft", "--auto", input_fasta]
    return subprocess.run(mafft_cmd, capture_output=True, text=True)


def run_trimal(input_aln: str, output_phy: str, tree_method: str = "raxml"):
    # trimal command
    if tree_method == "raxml":
        trimal_cmd = [
            "trimal",
            "-in",
            input_aln,
            "-phylip",
            "-out",
            output_phy,
        ]
    else:
        # output fasta
        trimal_cmd = [
            "trimal",
            "-in",
            input_aln,
            "-out",
            output_phy,
        ]
    return subprocess.run(trimal_cmd, capture_output=True, text=True)


# Function to run RAxML
def run_raxml(alignment_file: str, output_dir: str):
    # RAxML command
    raxml_cmd = [
        "raxmlHPC-PTHREADS-SSE3",
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


def run_fasttree(alignment_file: str, output_tree: str):
    """Run FastTre on nucleotide alignment file and save the tree by taking input from
    stdin and redirecting stdout to a file.
    """
    # FastTree command
    fasttree_cmd = [
        "FastTree",
        "-nt",
        "-gtr",
    ]
    print(f"Running FastTree with command: {' '.join(fasttree_cmd)}")
    # return subprocess.run(fasttree_cmd, capture_output=True, text=True, cwd=".")
    with open(alignment_file, "r") as aln, open(output_tree, "w") as tree:
        return subprocess.run(fasttree_cmd, stdin=aln, stdout=tree, text=True, cwd=".")


def process_sequences(
    zotu_fasta_path: str,
    alignment_out: str,
    tree_out: str,
    seq_names: list[str] | None = None,
    tree_method: str = "raxml",
) -> int:
    if seq_names is None:
        input_fasta = zotu_fasta_path
    else:
        # with tempfile.NamedTemporaryFile(
        #     mode="w", delete=False, suffix=".fasta"
        # ) as temp_fasta:
        temp_fasta = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta")
        # Filter and write sequences to a temp file
        _ = SeqIO.write(
            list(filter_sequences(zotu_fasta_path, seq_names)),
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
    temp_aln.flush()

    # Run trimal on the alignment file and check if it ran successfully
    trimal_result = run_trimal(temp_aln.name, alignment_out, tree_method)
    temp_aln.close()
    if trimal_result.returncode != 0:
        print("Error in running trimal:", trimal_result.stderr)
        return 1
    else:
        print(
            "MAFFT and trimal completed successfully. Alignment saved as:",
            alignment_out,
        )

    if tree_method == "raxml":
        # Run RAxML, output to a temporary directory, copy the tree file to the output path
        with tempfile.TemporaryDirectory() as temp_dir:
            raxml_result = run_raxml(alignment_out, temp_dir)
            tree_file = os.path.join(temp_dir, "RAxML_bestTree.whatever")
        # Check if RAxML ran successfully
        if raxml_result.returncode != 0:
            print("Error in running RAxML:", raxml_result.stderr)
            return 1
        else:
            print(
                "RAxML completed successfully. Tree saved in Newick format as:",
                tree_out,
            )
            shutil.copy(tree_file, tree_out)
    elif tree_method == "fasttree":
        fasttree_result = run_fasttree(alignment_out, tree_out)
        if fasttree_result.returncode != 0:
            print("Error in running FastTree:", fasttree_result.stderr)
            return 1
        else:
            print(
                "FastTree completed successfully. Tree saved in Newick format as:",
                tree_out,
            )
    else:
        raise ValueError("Invalid tree method:", tree_method)

    if seq_names is not None:
        temp_fasta.close()

    return 0


def main() -> None:
    # if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a phylogenetic tree from a FASTA file of sequences and a list of sequence names"
    )
    _ = parser.add_argument(
        "-i",
        "--input_fasta",
        help="Input FASTA file",
        required=True,
        type=str,
        metavar="FASTA",
    )
    _ = parser.add_argument(
        "-m",
        "--tree_method",
        help="Method to generate the tree",
        required=False,
        type=str,
        default="fasttree",
        choices=["raxml", "fasttree"],
    )
    args = parser.parse_args()

    # Read the sequence names file

    # Process the sequences
    # prefix is dirname + anything before first dot
    prefix = os.path.join(
        os.path.dirname(args.input_fasta),
        os.path.basename(args.input_fasta).split(".")[0],
    )
    aln_file = prefix + (".phy" if args.tree_method == "raxml" else ".alnfna")
    tree_file = prefix + ".newick"
    process_sequences(args.input_fasta, aln_file, tree_file, tree_method="fasttree")

    # Read the tree file and print the tree
    # tree = Phylo.read(args.output, "newick")
    # Phylo.draw_ascii(tree)
