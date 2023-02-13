"""
Module for downloading FASTA and Genbank files of the first result from querying the NCBI database.

Usage example:
    download_gene_files(
        email="majerdaniel93@gmail.com",
        species="Drosophila melanogaster",
        gene_full_name="alcohol dehydrogenase",
        filter="mRNA",
    )

"""

from Bio import Entrez, SeqIO
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


def __save_efetch_result(
    output_dir: str, file_name: str, file_extension: str, content: TextIO
) -> None:
    """
    Creates the output directory if doesn't exists,
    then saves the content TextIO type content argument.
    """

    # Save the files to the local system
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    with open(
        os.path.join(output_dir, f'{file_name}.{file_extension.lstrip(".")}'), "w"
    ) as f:
        f.write(content.read())


def download_gene_files(
    database: str, query: str, email: str, output_dir: str, file_name: str
) -> None:
    """
    This function downloads the FASTA and Genbank files of a given Gene and Species.
    """

    # Provide email to NCBI to use the Entrez API
    Entrez.email = email

    # Search for the gene in the NCBI database
    with Entrez.esearch(db=database, term=query) as handle:
        search_results: SearchResult = dict(Entrez.read(handle))  # type: ignore

    # Check if there are any search results
    if int(search_results["Count"]) == 0:
        print("No results found for the given gene and species.")
        return

    # Fetch the first result and download the FASTA and Genbank files
    gb_result = Entrez.efetch(
        db=database, id=search_results["IdList"][0], rettype="gb", retmode="text"
    )

    # print(list(SeqIO.parse(gb_result, 'genbank'))[0].annotations)

    fasta_result = Entrez.efetch(
        db=database, id=search_results["IdList"][0], rettype="fasta", retmode="text"
    )

    # Write genbank.
    __save_efetch_result(output_dir, file_name, "gb", gb_result)

    # Write fasta.
    __save_efetch_result(output_dir, file_name, "fasta", fasta_result)

    print(
        f"""
Files for the following query has been successfully downloaded to the {output_dir} directory:
- {query}
        """
    )


if __name__ == "__main__":

    # Example usage:
    download_gene_files(
        database="protein",
        query='alcohol dehydrogenase[Protein Name] AND "Drosophila melanogaster"[Organism] AND refseq[filter]',
        email="majerdaniel93@gmail.com",
        output_dir="../data",
        file_name="alcohol_dehydrogenase-drosophila_melanogaster",
    )
