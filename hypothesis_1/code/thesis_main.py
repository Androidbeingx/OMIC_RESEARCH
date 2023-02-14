"""
Main script of the first thesis.
"""

from download_ncbi_query import download_gene_files

# Global types
FileName = str
NCBIQuery = str


def download_fasta_and_genbank() -> None:

    queries: dict[FileName, NCBIQuery] = {
        "alcohol_dehydrogenase-drosophila_melanogaster": 'alcohol dehydrogenase[Protein Name] AND "Drosophila melanogaster"[Organism] AND refseq[filter]',
        "alcohol_dehydrogenase-drosophila_simulans": 'alcohol dehydrogenase[Protein Name] AND "Drosophila simulans"[Organism] AND refseq[filter]',
        "alcohol_dehydrogenase-drosophila_sechellia": 'alcohol dehydrogenase[Protein Name] AND "Drosophila sechellia"[Organism] AND refseq[filter]',
        "alcohol_dehydrogenase-drosophila_yakuba": 'alcohol dehydrogenase[Protein Name] AND "Drosophila yakuba"[Organism] AND refseq[filter]',
        "alcohol_dehydrogenase-drosophila_erecta": 'alcohol dehydrogenase[Protein Name] AND "Drosophila erecta"[Organism] AND refseq[filter]',
        "alcohol_dehydrogenase-drosophila_ananassae": 'alcohol dehydrogenase[Protein Name] AND "Drosophila ananassae"[Organism] AND refseq[filter]',
        "alcohol_dehydrogenase-drosophila_suzukii": 'alcohol dehydrogenase[Protein Name] AND "Drosophila suzukii"[Organism] AND refseq[filter]',
        "alcohol_dehydrogenase-drosophila_teissieri": 'alcohol dehydrogenase[Protein Name] AND "Drosophila teissieri"[Organism] AND refseq[filter]',
        "alcohol_dehydrogenase-homo_sapiens": 'ADH1B[Gene Name] AND "Homo sapiens"[Organism] AND refseq[filter]'
    }

    for file_name, ncbi_query in queries.items():

        download_gene_files(
            database="protein",
            query=ncbi_query,
            email="majerdaniel93@gmail.com",
            output_dir="../data",
            file_name=file_name,
        )


def main() -> None:

    download_fasta_and_genbank()

if __name__ == '__main__':
    main()
