"""
NOT FINISHED YET!
"""

from Bio import Entrez, SeqIO
from Bio.Entrez import Parser
from typing import TypedDict
import os


class SearchResult(TypedDict):
    Count: str
    RetMax: str
    RetStart: str
    IdList: list[str]
    TranslationSet: list[dict[str, str]]
    TranslationStack: list[dict[str, str]]
    QueryTranslation: str

def download_gene_files(email: str, gene: str, species: str) -> None:
    """
    This function downloads the RefSeq FASTA and Genbank files of a given Gene and Species.
    """

    # Provide email to NCBI to use the Entrez API
    Entrez.email = email

    # Search for the gene in the NCBI database
    search_query: str = f"{gene}[Gene] AND {species}[Organism]"

    with Entrez.esearch(db="nucleotide", term=search_query) as handle:
        search_results: SearchResult = dict(Entrez.read(handle)) #type: ignore

    # Check if there are any search results
    if int(search_results["Count"]) == 0:
        print("No results found for the given gene and species.")
        return

    # Fetch the first result and download the RefSeq FASTA and Genbank files
    gb_result = Entrez.efetch(
        db="nucleotide", id=search_results["IdList"][0], rettype="gb", retmode="text"
    )

    fasta_result = Entrez.efetch(
        db="nucleotide", id=search_results["IdList"][0], rettype="fasta", retmode="text"
    )

    # Save the files to the local system
    if not os.path.exists("output"):
        os.mkdir("output")

    # Write genbank.
    gb_file_name: str = f"{gene}-{species.lower().replace(' ','_')}.gb"
    with open(os.path.join("output", gb_file_name), "w") as f:
        f.write(gb_result.read())

    # Write fasta.
    fasta_file_name: str = f"{gene}-{species.lower().replace(' ','_')}.fasta"
    with open(os.path.join("output", fasta_file_name), "w") as f:
        f.write(fasta_result.read())

    print(
        f"RefSeq FASTA and Genbank files for {gene} in {species} are saved to the output directory."
    )

# Example usage
def main() -> int:
    download_gene_files("majerdaniel93@gmail.com", "ADH", "Drosophila melanogaster")
    return 0

if __name__ == '__main__':
    main()
