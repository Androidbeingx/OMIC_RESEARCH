"""
Main script of the first thesis.
"""

import re
import itertools
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from download_ncbi_query import download_gene_files
from aligments_fasta import use_blosum62_matrix, AlignmentResults


# Global types
# ==============================
FileName = str
NCBIQuery = str

ProteinName = str
Alignments = dict[tuple[ProteinName, ProteinName], AlignmentResults]
# ==============================

def download_fasta_and_genbank(queries: dict[FileName, NCBIQuery]) -> list[str]:
    """Download FASTA and GenBank files from NCBI based on given queries."""

    # Download gene files from NCBI database for each query and save the paths for successful downloads
    result_paths: list[str] = []
    for file_name, ncbi_query in queries.items():

        # Download the gene files and get their paths
        genbank_path, fasta_path = download_gene_files(
            database="protein",
            query=ncbi_query,
            email="majerdaniel93@gmail.com",
            output_dir="../data/genbanks-and-fastas",
            file_name=file_name,
        )

        # Check if the download was successful and add the paths to the result_paths list
        if genbank_path is not None and fasta_path is not None:
            result_paths.append(genbank_path)
            result_paths.append(fasta_path)

    # Return the list of paths for the successfully downloaded files
    return result_paths


def filter_fasta_paths(paths: list[str]) -> list[str]:
    """
    Filter a list of file paths to only include those with a '.fasta' file extension.
    """

    # Filter the given list of file paths to only keep the paths with .fasta extension
    return [path for path in paths if path.lower().endswith(".fasta")]


def get_alignments(file_paths: tuple[str, ...] | list[str]) -> Alignments:
    """Get pairwise alignments for protein sequences in the given files."""

    # Generate pairwise alignments between all the given files
    # Combinations of all possible file path pairs
    file_path_combinations: tuple[tuple[str, str]] = tuple(
        combination for combination in itertools.combinations(file_paths, 2)
    )
    # Create an empty dictionary to store the alignment results
    alignments: Alignments = {}

    # Loop through all the file path pairs and perform pairwise alignment for each pair
    for paths in file_path_combinations:

        # Read the protein sequences from the fasta files
        protein_1: SeqRecord = SeqIO.read(paths[0], "fasta")
        protein_2: SeqRecord = SeqIO.read(paths[1], "fasta")

        # Extract the protein names from the sequence descriptions
        protein_name_1: str = " ".join(protein_1.description.split(" ")[1:])
        protein_name_2: str = " ".join(protein_2.description.split(" ")[1:])

        # Calculate the alignment score and the aligned sequences for the given protein sequences
        alignment_result: AlignmentResults = use_blosum62_matrix(
            protein_1.seq, protein_2.seq
        )
        # Store the alignment result in the dictionary using the protein names as keys
        alignments[protein_name_1, protein_name_2] = alignment_result

    # Return the dictionary containing the alignment results
    return alignments


def save_alignments_as_files(alignments: Alignments) -> None:

    # TODO: This is to build header annotations for the alignments.
    # NAME: <name of gene / protein>
    # SPECIES: <taxonomy name>
    # TYPE: <If it's protein, nc o what..>

    taxonomy_name_regex = re.compile(r'(\[([^\[\]]*)\])$')
    gene_name_regex = re.compile(r'^[a-z\-\s]+[^\[]')

    for pair in list(alignments.keys()):
        # print(str(pair))
        print(list(taxonomy_name_regex.finditer(str(pair[0])))[0].group(2))
        # print(list(gene_name_regex.finditer(str(pair[0]))))

    # x = alignments[)[0]]['alignment'][0]
    # print(x)

    # with open('testy.alignment', 'w') as f:
    #     f.write(str(x))


def main() -> None:

    # Define the NCBI queries for each gene of interest
    queries: dict[FileName, NCBIQuery] = {
        "alcohol_dehydrogenase-drosophila_melanogaster": 'alcohol dehydrogenase[Protein Name] AND "Drosophila melanogaster"[Organism] AND refseq[filter]',
        "alcohol_dehydrogenase-drosophila_simulans": 'alcohol dehydrogenase[Protein Name] AND "Drosophila simulans"[Organism] AND refseq[filter]',
        "alcohol_dehydrogenase-drosophila_sechellia": 'alcohol dehydrogenase[Protein Name] AND "Drosophila sechellia"[Organism] AND refseq[filter]',
        "alcohol_dehydrogenase-drosophila_yakuba": 'alcohol dehydrogenase[Protein Name] AND "Drosophila yakuba"[Organism] AND refseq[filter]',
        "alcohol_dehydrogenase-drosophila_erecta": 'alcohol dehydrogenase[Protein Name] AND "Drosophila erecta"[Organism] AND refseq[filter]',
        "alcohol_dehydrogenase-drosophila_ananassae": 'alcohol dehydrogenase[Protein Name] AND "Drosophila ananassae"[Organism] AND refseq[filter]',
        "alcohol_dehydrogenase-drosophila_suzukii": 'alcohol dehydrogenase[Protein Name] AND "Drosophila suzukii"[Organism] AND refseq[filter]',
        "alcohol_dehydrogenase-drosophila_teissieri": 'alcohol dehydrogenase[Protein Name] AND "Drosophila teissieri"[Organism] AND refseq[filter]',
        "alcohol_dehydrogenase-homo_sapiens": 'ADH1B[Gene Name] AND "Homo sapiens"[Organism] AND refseq[filter]',
    }

    # These paths are for testing. Simply avoiding the unnecessary repeated downloads.
    # ==============================================================================
    fasta_paths: tuple[str, ...] = (
        "../data/genbanks-and-fastas/alcohol_dehydrogenase-homo_sapiens.fasta",
        "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_teissieri.fasta",
        "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_suzukii.fasta",
        "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_ananassae.fasta",
        "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_erecta.fasta",
        "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_yakuba.fasta",
        "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_sechellia.fasta",
        "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_simulans.fasta",
        "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_melanogaster.fasta"
    )
    # ==============================================================================

    # # Download the genbank and fasta files of each first result / query.
    # paths: list[str] = download_fasta_and_genbank(queries)

    # # Extract the paths for only the fasta files.
    # fasta_paths: list[str] = filter_fasta_paths(paths)

    # Get alignment results.
    alignments: Alignments = get_alignments(fasta_paths)
    save_alignments_as_files(alignments)

if __name__ == "__main__":
    main()
