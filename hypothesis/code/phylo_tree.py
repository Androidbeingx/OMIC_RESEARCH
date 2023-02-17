import os
import subprocess
from Bio import Phylo


def join_fasta_files(fasta_files: list[str] | tuple[str, ...], output_dir: str, file_name: str) -> str:
    """
    Combine the contents of multiple input FASTA files into a single output file.
    """

    # Check if any of the input files don't exist
    for file_path in fasta_files:
        if not os.path.exists(file_path):
            raise IOError(f"File not found: {file_path}")

    # Create the output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)

    # Construct the output file path
    output_file_path: str = os.path.join(output_dir, file_name)

    # Create the output file and copy the contents of the input files to it
    with open(output_file_path, "w") as out_file:
        for file_path in fasta_files:
            # Open each input FASTA file and copy its contents to the output file
            with open(file_path, "r") as in_file:
                out_file.write(in_file.read())

    print(f"Joined {len(fasta_files)} FASTA files into {output_dir}")

    return output_file_path


def run_clustal(input_file_path: str, output_msa_path: str, output_tree_path: str) -> None:
    """
    Run the Clustal Omega multiple sequence alignment tool on the specified input file.
    """

    try:
        subprocess.run(['clustalo', '-i', input_file_path, '-o', output_msa_path, f'--guidetree-out={output_tree_path}'])
    except FileNotFoundError as err:
        print(err)
        print('Could not execute Clustalo on the system. Please, check if it is intalled.')


def display_phylo_tree(tree_path: str, format: str) -> None:
    """
    Display a phylogenetic tree in ASCII or graphical format.
    """

    # Read the phylogenetic tree from file using Bio.Phylo
    tree = Phylo.read(tree_path, "newick")
    Phylo.draw_ascii(tree)
    Phylo.draw(tree)


def main() -> None:
    "Example of usage"

    fasta_files = [
    "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_guanche.fasta",
    "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_sechellia.fasta",
    "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_simulans.fasta",
    "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_ananassae.fasta",
    "../data/genbanks-and-fastas/alcohol_dehydrogenase-homo_sapiens.fasta",
    "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_persimilis.fasta",
    "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_teissieri.fasta",
    "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_melanogaster.fasta",
    "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_subobscura.fasta",
    "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_yakuba.fasta",
    "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_suzukii.fasta",
    "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_erecta.fasta"
    ]

    multi_fasta_path: str = join_fasta_files(fasta_files, '../data/multi-fasta', 'multi-fasta.fasta')
    output_msa_path: str = '../data/multi-fasta/msa_result.clustal'
    output_tree_path: str = '../data/multi-fasta/tree.dnd'

    run_clustal(multi_fasta_path, output_msa_path, output_tree_path)
    display_phylo_tree(output_tree_path, 'newick')

if __name__ == '__main__':
    main()
