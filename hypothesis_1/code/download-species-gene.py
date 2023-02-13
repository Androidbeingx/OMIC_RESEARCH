"""
Module for downloading FASTA and Genbank files of a given Gene and Species from NCBI.

The module also includes two helper functions:
- __generate_query_from_args builds a query from the given gene symbol, gene full name, species, and filter, to be
sent against the NCBI database.

- __assign_gene checks if a gene symbol or full name has been given as an argument and assigns the gene accordingly.

Usage example:
    download_gene_files(
        email="majerdaniel93@gmail.com",
        species="Drosophila melanogaster",
        gene_full_name="alcohol dehydrogenase",
        filter="mRNA",
    )

"""

from Bio import Entrez
from typing import TypedDict, TextIO
import os


class SearchResult(TypedDict):
    Count: str
    RetMax: str
    RetStart: str
    IdList: list[str]
    TranslationSet: list[dict[str, str]]
    TranslationStack: list[dict[str, str]]
    QueryTranslation: str


def __generate_query_from_args(
    gene_symbol: str, gene_full_name: str, species: str, filter: str
) -> str:
    """
    This function is a sub-functions of "download_gene_files(...)".
    Builds a query to be sent against NCBI DB.
    Argument for "args" parameter MUST contain the following keys:
        - gene_symbol: str
        - gene_full_name: str
        - species: str
        - filter
    """

    search_query: str = f""

    if gene_symbol != "" and gene_full_name != "":
        raise ValueError("""Can't use both "gene_symbol" and "gene_full_name"!""")
    elif gene_symbol == "" and gene_full_name == "":
        raise ValueError(
            """Both "gene_symbol" and "gene_full_name" have empty values! Must use one of them!"""
        )
    elif gene_symbol != "":
        search_query += f"{gene_symbol}[Gene]"
    elif gene_full_name != "":
        search_query += f"{gene_full_name} gene"

    search_query += f" AND {species}[Organism]"

    if filter != "":
        search_query += search_query + f" AND {filter} [Filter]"

    return search_query


def __assign_gene(gene_symbol: str, gene_full_name: str) -> str:
    """
    This function is a sub-functions of "download_gene_files(...)".
    Checks if gene symbol of full name has been given as an argument.
    Argument for "args" parameter MUST contain the following keys:
        - gene_symbol: str
        - gene_full_name: str
    """

    gene: str = ""

    if gene_symbol != "" and gene_full_name != "":
        raise ValueError("""Can't use both "gene_symbol" and "gene_full_name"!""")
    elif gene_symbol == "" and gene_full_name == "":
        raise ValueError(
            """Both "gene_symbol" and "gene_full_name" have empty values! Must use one of them!"""
        )
    elif gene_symbol != "":
        gene = gene_symbol
    elif gene_full_name != "":
        gene = gene_full_name

    return gene


def download_gene_files(
    email: str,
    species: str,
    gene_symbol: str = "",
    gene_full_name: str = "",
    filter: str = "",
) -> None:
    """
    This function downloads the FASTA and Genbank files of a given Gene and Species.
    """

    # Provide email to NCBI to use the Entrez API
    Entrez.email = email

    # Build query.
    search_query: str = __generate_query_from_args(
        gene_symbol, gene_full_name, species, filter
    )

    # Assign gene
    gene: str = __assign_gene(gene_symbol, gene_full_name)

    # Search for the gene in the NCBI database
    with Entrez.esearch(db="nucleotide", term=search_query) as handle:
        search_results: SearchResult = dict(Entrez.read(handle))  # type: ignore

    # Check if there are any search results
    if int(search_results["Count"]) == 0:
        print("No results found for the given gene and species.")
        return

    # Fetch the first result and download the FASTA and Genbank files
    gb_result = Entrez.efetch(
        db="nucleotide", id=search_results["IdList"][0], rettype="gb", retmode="text"
    )

    fasta_result = Entrez.efetch(
        db="nucleotide", id=search_results["IdList"][0], rettype="fasta", retmode="text"
    )

    def __save_output(
        output_dir: str, gene: str, content: TextIO, file_ext: str
    ) -> None:
        # Save the files to the local system
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        file_name: str = f"{gene.lower().replace(' ','_')}-{species.lower().replace(' ','_')}.{file_ext.strip('.')}"
        with open(os.path.join("output", file_name), "w") as f:
            f.write(content.read())

    # Write genbank.
    __save_output("output", gene, gb_result, "gb")

    # Write fasta.
    __save_output("output", gene, fasta_result, "fasta")

    print(
        f"FASTA and Genbank files for {gene} in {species} are saved to the output directory."
    )


if __name__ == "__main__":

    # Example usage
    def main() -> int:

        download_gene_files(
            email="majerdaniel93@gmail.com",
            species="Drosophila melanogaster",
            gene_full_name="alcohol dehydrogenase",
            filter="mRNA",
        )
        return 0

    main()
