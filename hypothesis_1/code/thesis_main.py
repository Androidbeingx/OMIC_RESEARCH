"""
Main script of the first thesis.
"""

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
AlignmentScores = dict[tuple[ProteinName, ProteinName], AlignmentResults]
# ==============================


def download_fasta_and_genbank(queries: dict[FileName, NCBIQuery]) -> list[str]:

    result_paths: list[str] = []
    for file_name, ncbi_query in queries.items():

        genbank_path, fasta_path = download_gene_files(
            database="protein",
            query=ncbi_query,
            email="majerdaniel93@gmail.com",
            output_dir="../data",
            file_name=file_name,
        )

        if genbank_path is not None and fasta_path is not None:
            result_paths.append(genbank_path)
            result_paths.append(fasta_path)

    return result_paths

def filter_fasta_paths(paths: list[str]) -> list[str]:
    return [path for path in paths if path.lower().endswith('.fasta')]


def get_alignments(file_paths: tuple[str, ...] | list[str]) -> dict:

    file_path_combinations: tuple[tuple[str, str]] = tuple(
        combination for combination in itertools.combinations(file_paths, 2)
    )

    alignments: AlignmentScores = {}

    for paths in file_path_combinations:
        protein_1: SeqRecord = SeqIO.read(paths[0], "fasta")
        protein_2: SeqRecord = SeqIO.read(paths[1], "fasta")

        protein_name_1: str = " ".join(protein_1.description.split(" ")[1:])
        protein_name_2: str = " ".join(protein_2.description.split(" ")[1:])

        alignment_result: AlignmentResults = use_blosum62_matrix(protein_1.seq, protein_2.seq)

        alignments[protein_name_1, protein_name_2] = alignment_result

    return alignments


def main() -> None:

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
    # fasta_paths: tuple[str, ...] = (
    #     "../data/genbanks-and-fastas/alcohol_dehydrogenase-homo_sapiens.fasta",
    #     "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_teissieri.fasta",
    #     "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_suzukii.fasta",
    #     "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_ananassae.fasta",
    #     "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_erecta.fasta",
    #     "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_yakuba.fasta",
    #     "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_sechellia.fasta",
    #     "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_simulans.fasta",
    #     "../data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_melanogaster.fasta"
    # )
    # ==============================================================================


    paths: list[str] = download_fasta_and_genbank(queries)
    fasta_paths: list[str] = filter_fasta_paths(paths)

    alignments: AlignmentScores = get_alignments(fasta_paths)
    print(alignments)


if __name__ == "__main__":
    main()
