"""
Main script of the first thesis.
"""

import re
import os
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


def get_alignments(target_fasta: str, file_paths: tuple[str, ...] | list[str]) -> Alignments:
    """Get pairwise alignments for protein sequences in the given files."""

    # Generate pairwise alignments between all the given files
    # Combinations of all possible file path pairs
    # file_path_combinations: tuple[tuple[str, str]] = tuple(
    #     combination for combination in itertools.combinations(file_paths, 2)
    # )

    # Filter out our target fasta from the rest.
    filtered_paths: list[str] = [path for path in file_paths if path != target_fasta]

    # Create an empty dictionary to store the alignment results
    alignments: Alignments = {}

    # Loop through all the file path pairs and perform pairwise alignment for each pair
    for path in filtered_paths:

        # Read the protein sequences from the fasta files
        protein_1: SeqRecord = SeqIO.read(target_fasta, "fasta")
        protein_2: SeqRecord = SeqIO.read(path, "fasta")

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


def save_alignments_as_files(alignments: Alignments, output_dir: str) -> list[str]:
    '''
    Save the protein sequence alignments as text files.
    '''

    # Regular expressions for extracting taxonomy name and gene name from the input
    taxonomy_name_regex = re.compile(r'(\[([^\[\]]*)\])$')
    gene_name_regex = re.compile(r'^[a-z\-\s]+[^\[,]')

    # A list to store the output file paths
    output_file_paths: list[str] = []

    # Iterate over the alignments and save each alignment as a separate file
    for gene_pair, result in alignments.items():
        # Extract the taxonomy names and gene name from the input using regex
        taxonomy_name_1: str = list(taxonomy_name_regex.finditer(str(gene_pair[0])))[0].group(2)
        taxonomy_name_2: str = list(taxonomy_name_regex.finditer(str(gene_pair[1])))[0].group(2)
        gene_name: str = list(gene_name_regex.finditer(str(gene_pair[1])))[0].group(0).strip('_')

        # Get scores
        score: float = result['score']
        percent_identity: float = result['percent_identity']

        # Extract the alignment content from the result dictionary
        alignment_content: str = str(result['alignment'][0])

        # Format the taxonomy names and gene name for creating the output file name
        formated_taxonomy_name_1: str =  taxonomy_name_1.replace(' ', '_').lower()
        formated_taxonomy_name_2: str =  taxonomy_name_2.replace(' ', '_').lower()
        formated_gene_name: str = gene_name.replace(' ', '_').lower()

        # Create the output file name and path
        file_name: str = f"{formated_gene_name}-{formated_taxonomy_name_1}-{formated_taxonomy_name_2}.alignment"
        file_path: str = os.path.join(output_dir, file_name)

        # Create the output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        # Write the alignment content to the output file
        with open(file_path, 'w') as file:
            file.write(
                f"GENE NAME: {gene_name}\n" +\
                f"TAXONOMY NAME 1: {taxonomy_name_1}\n" +\
                f"TAXONOMY NAME 2: {taxonomy_name_2}\n" +\
                f"TYPE: PROTEIN\n" +\
                f"SCORE: {score}\n" +\
                f"PERCENT IDENTITY: {percent_identity}\n\n" +\
                f"{alignment_content}"
            )

        # Add the file path to the output list
        output_file_paths.append(file_path)

    return output_file_paths


def main() -> None:

    # Define the NCBI queries for each gene of interest
    queries: dict[FileName, NCBIQuery] = {
        "alcohol_dehydrogenase-drosophila_persimilis": 'alcohol dehydrogenase[Protein Name] AND "Drosophila persimilis"[Organism] AND refseq[filter]',
        "alcohol_dehydrogenase-drosophila_guanche": 'alcohol dehydrogenase[Protein Name] AND "Drosophila guanche"[Organism] AND refseq[filter]',
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
        '../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_guanche.fasta',
        '../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_sechellia.fasta',
        '../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_simulans.fasta',
        '../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_ananassae.fasta',
        '../data/genbanks-and-fastas/alcohol_dehydrogenase-homo_sapiens.fasta',
        '../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_persimilis.fasta',
        '../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_teissieri.fasta',
        '../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_melanogaster.fasta',
        '../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_subobscura.fasta',
        '../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_yakuba.fasta',
        '../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_suzukii.fasta',
        '../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_erecta.fasta',
    )
    # ==============================================================================

    # # Download the genbank and fasta files of each first result / query.
    # paths: list[str] = download_fasta_and_genbank(queries)

    # # Extract the paths for only the fasta files.
    # fasta_paths: list[str] = filter_fasta_paths(paths)

    target_fasta: str = '../data/genbanks-and-fastas/alcohol_dehydrogenase-homo_sapiens.fasta'

    # Get alignment results.
    alignments: Alignments = get_alignments(target_fasta, fasta_paths)
    alignment_paths: list[str] = save_alignments_as_files(alignments, '../data/alignments')

    # This is for the end!!!
    # Leave it commented for now.
    # ===========================
    # multi_fasta_path: str = join_fasta_files(fasta_paths, '../data/multi-fasta', 'multi-fasta-result.fasta')
    # ===========================

if __name__ == "__main__":
    main()
