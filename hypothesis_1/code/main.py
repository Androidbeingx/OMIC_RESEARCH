"""
Main script of the first thesis.
"""
import itertools
import os
import re
from typing import TypedDict

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from aligments_fasta import AlignmentResults
from aligments_fasta import use_blosum62_matrix
from annotation_regexs import ANNOTATION_REGEXS
from annotation_regexs import get_info_with_regex
from download_ncbi_query import download_gene_files
from phylo_tree import join_fasta_files, run_clustal, display_phylo_tree

# Global types
# ==============================
FileName = str
NCBIQuery = str

ProteinName = str
Alignments = dict[tuple[ProteinName, ProteinName], AlignmentResults]


class AlignmentResultTable(TypedDict):
    accession: list[str]  # NCBI accession number
    organism: list[str]  # Organism name
    chromosome: list[str]  # Chromosome number
    product: list[str]  # Gene product description
    gene: list[str]  # Gene name
    cds: list[str]  # CDS sequence
    definition: list[str]  # Gene definition
    score: list[float]  # Alignment score
    identity_percent: list[float]  # Alignment percent identity

# ==============================


TAXA_NAME_FROM_DEFINITION_REGX = re.compile(r"(\[([^\[\]]*)\])$")
GENE_NAME_FROM_DEFINITION_REGX = re.compile(r"^[a-z\-\s]+[^\[,]")


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


def get_files_with_extension(directory: str, extension: str) -> list[str]:
    """
    Gets files from directory with the given extension.
    """
    files: list[str] = []
    for file in os.listdir(directory):
        if file.endswith(f'.{extension.strip(".")}'):
            files.append(os.path.join(directory, file))

    return files


def get_alignments(
    target_fasta: str, file_paths: tuple[str, ...] | list[str]
) -> Alignments:
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
    """
    Save the protein sequence alignments as text files.
    """

    # A list to store the output file paths
    output_file_paths: list[str] = []

    # Iterate over the alignments and save each alignment as a separate file
    for gene_pair, result in alignments.items():
        # Extract the taxonomy names and gene name from the input using regex
        taxonomy_name_1: str = list(
            TAXA_NAME_FROM_DEFINITION_REGX.finditer(str(gene_pair[0]))
        )[0].group(2)
        taxonomy_name_2: str = list(
            TAXA_NAME_FROM_DEFINITION_REGX.finditer(str(gene_pair[1]))
        )[0].group(2)
        gene_name: str = (
            list(GENE_NAME_FROM_DEFINITION_REGX.finditer(str(gene_pair[1])))[0]
            .group(0)
            .strip("_")
        )

        # Get scores
        score: float = result["score"]
        percent_identity: float = result["percent_identity"]

        # Extract the alignment content from the result dictionary
        alignment_content: str = str(result["alignment"][0])

        # Format the taxonomy names and gene name for creating the output file name
        formated_taxonomy_name_1: str = taxonomy_name_1.replace(" ", "_").lower()
        formated_taxonomy_name_2: str = taxonomy_name_2.replace(" ", "_").lower()
        formated_gene_name: str = gene_name.replace(" ", "_").lower()

        # Create the output file name and path
        file_name: str = f"{formated_gene_name}-{formated_taxonomy_name_1}-{formated_taxonomy_name_2}.alignment"
        file_path: str = os.path.join(output_dir, file_name)

        # Create the output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        # Write the alignment content to the output file
        with open(file_path, "w") as file:
            file.write(
                f"GENE NAME: {gene_name}\n"
                + f"TAXONOMY NAME 1: {taxonomy_name_1}\n"
                + f"TAXONOMY NAME 2: {taxonomy_name_2}\n"
                + f"TYPE: PROTEIN\n"
                + f"SCORE: {score}\n"
                + f"PERCENT IDENTITY: {percent_identity}\n\n"
                + f"{alignment_content}"
            )

        # Add the file path to the output list
        output_file_paths.append(file_path)

    return output_file_paths


def save_dataframe(output_dir: str, file_name: str, dataframe: pd.DataFrame) -> str:
    """
    Saves the given dataframe to output path as CSV, separated by ";".
    """
    os.makedirs(output_dir, exist_ok=True)
    file_path: str = os.path.join(output_dir, file_name)
    dataframe.to_csv(file_path, sep=";")

    return file_path


def create_dataframe_data(
    file_paths: list[str] | tuple[str, ...],
    alignments: Alignments,
    target_organism: str,
) -> AlignmentResultTable:
    """
    Create dictionary to identify what each value is.
    (* The function will ignore the "target_organism".)

    Return: Information organized with your key.
    """

    # Initialize an empty dictionary to store the output data
    dict_info: AlignmentResultTable = {
        "accession": [],  # NCBI accession number
        "organism": [],  # Organism name
        "chromosome": [],  # Chromosome number
        "product": [],  # Gene product description
        "gene": [],  # Gene name
        "cds": [],  # CDS sequence
        "definition": [],  # Gene definition
        "score": [],  # Alignment score
        "identity_percent": [],  # Alignment percent identity
    }

    # Extract information from annotation files
    for path in file_paths:
        annotations: list[str] = get_info_with_regex(path, ANNOTATION_REGEXS)
        # Ignore records that belong to the target organism
        if annotations[1].strip().lower() != target_organism.strip().lower():
            dict_info["accession"].append(annotations[0])
            dict_info["organism"].append(annotations[1])
            dict_info["chromosome"].append(annotations[2])
            dict_info["product"].append(annotations[3])
            dict_info["gene"].append(annotations[4])
            dict_info["cds"].append(annotations[5])
            dict_info["definition"].append(annotations[6])

    # Extract information from the alignment data
    for gene_pair, info in alignments.items():
        # Extract the taxonomy name of the query sequence
        query_taxa_name: str = list(
            TAXA_NAME_FROM_DEFINITION_REGX.finditer(str(gene_pair[1]))
        )[0].group(2)

        # Find the corresponding alignment score and percent identity for each organism in the output data
        for organism in dict_info["organism"]:
            if query_taxa_name.strip().lower() == organism.strip().lower():
                dict_info["score"].append(info["score"])
                dict_info["identity_percent"].append(info["percent_identity"])

    # Return the output dictionary
    return dict_info


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

    # ==============================================================================

    # # Download the genbank and fasta files of each first result / query.
    paths: list[str] = download_fasta_and_genbank(queries)

    # Get file paths with the given extensions.
    fasta_paths: list[str] = get_files_with_extension(
        "../data/genbanks-and-fastas", "fasta"
    )
    genbank_paths: list[str] = get_files_with_extension(
        "../data/genbanks-and-fastas", "gb"
    )

    # The file of the target sequence.
    target_fasta: str = (
        "../data/genbanks-and-fastas/alcohol_dehydrogenase-homo_sapiens.fasta"
    )

    # Get alignment results.
    alignments: Alignments = get_alignments(target_fasta, fasta_paths)
    alignment_paths: list[str] = save_alignments_as_files(
        alignments, "../data/alignments"
    )

    # Create dictionary of all the gathered informations to be further used in a dataframe.
    table_data: AlignmentResultTable = create_dataframe_data(
        genbank_paths, alignments, "Homo sapiens"
    )

    result_df: pd.DataFrame = (
        pd.DataFrame(table_data)
        .sort_values(by="score", ascending=False)
        .reset_index(drop=True)
    )
    result_df_path: str = save_dataframe(
        "../data/csv", "alignment-results.csv", result_df
    )

    # Phylogenetic tree from the result!
    # ===========================
    multi_fasta_path: str = join_fasta_files(fasta_paths, '../data/multi-fasta', 'multi-fasta.fasta')

    output_msa_path: str = '../data/multi-fasta/msa_result.clustal'
    output_tree_path: str = '../data/multi-fasta/tree.dnd'

    run_clustal(multi_fasta_path, output_msa_path, output_tree_path)
    display_phylo_tree(output_tree_path, 'newick')
    # ===========================


if __name__ == "__main__":
    main()
