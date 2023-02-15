import os
from io import StringIO
from Bio import AlignIO, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

def join_fasta_files(fasta_files: list[str] | tuple[str, ...], output_dir: str, file_name: str) -> str:
    """
    Joins multiple FASTA files into a single multi-FASTA file.
    Returns the file path of the output file.

    Example Usage:
        fasta_files = ['./seq1.fasta', './seq2.fasta', './seq3.fasta']
        output_file = './multi_seq.fasta'
        join_fasta_files(fasta_files, output_file)
    """

    # Check if any of the input files don't exist
    for file_path in fasta_files:
        if not os.path.exists(file_path):
            raise IOError(f"File not found: {file_path}")

    os.makedirs(output_dir, exist_ok=True)

    file_path: str = os.path.join(output_dir, file_name)

    # Create the output file
    with open(file_path, 'w') as out_file:
        for file_path in fasta_files:
            # Open each input FASTA file and copy its contents to the output file
            with open(file_path, 'r') as in_file:
                out_file.write(in_file.read())

    print(f"Joined {len(fasta_files)} FASTA files into {output_dir}")

    return file_path

from Bio import SeqIO

def adjust_sequence_lengths(filepath, length=None):
    """
    Given a filepath to a multi-FASTA file containing protein sequences, this
    function checks if all the sequences have the same length. If not, it
    truncates or pads the sequences to have the same length. If the argument
    'length' is not provided, the function will use the length of the longest
    sequence.

    Returns a string with the adjusted sequences in FASTA format.
    """
    # Parse the multi-FASTA file
    records = list(SeqIO.parse(filepath, "fasta"))

    # Find the length of the longest sequence
    max_length = max(len(record.seq) for record in records)

    # Set the desired length (if not provided)
    if length is None:
        length = max_length

    # Truncate or pad the sequences to have the same length
    for record in records:
        seq = record.seq
        if len(seq) < length:
            seq += "-" * (length - len(seq))
        elif len(seq) > length:
            seq = seq[:length]
        record.seq = seq

    # Return the adjusted sequences as a string in FASTA format
    return "".join(record.format("fasta") for record in records)

def create_phylogenetic_tree_from_file(filepath):
    """
    Given a filepath to a multi-FASTA file containing protein sequences, this function creates a phylogenetic tree
    using the Neighbor Joining algorithm.
    """

    # Create an alignment object from the protein sequences
    alignment = AlignIO.read(filepath, "fasta")

    # Create a distance matrix using the alignment
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)

    # Construct the phylogenetic tree using the Neighbor Joining algorithm
    constructor = DistanceTreeConstructor()

    # Build tree
    tree = constructor.nj(distance_matrix)

    # Draw the tree using the Phylo library
    Phylo.draw_ascii(tree)

if __name__ == '__main__':
    multi_fasta = adjust_sequence_lengths("../data/multi-fasta/multi-fasta-result.fasta")

    with open('testy.fasta', 'w') as f:
        f.write(multi_fasta)

    create_phylogenetic_tree_from_file('testy.fasta')
    # create_phylogenetic_tree_from_file("../data/multi-fasta/multi-fasta-result.fasta")
