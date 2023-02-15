"""
Main script of the first thesis.
"""

import itertools
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from download_ncbi_query import download_gene_files
from aligments_fasta import use_substitution_matrix


# Global types
# ==============================
FileName = str
NCBIQuery = str

ProteinName = str
AlignmentScore = float
AlignmentScores = dict[tuple[ProteinName, ProteinName], AlignmentScore]
# ==============================


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
        "alcohol_dehydrogenase-homo_sapiens": 'ADH1B[Gene Name] AND "Homo sapiens"[Organism] AND refseq[filter]',
    }

    for file_name, ncbi_query in queries.items():

        download_gene_files(
            database="protein",
            query=ncbi_query,
            email="majerdaniel93@gmail.com",
            output_dir="../data",
            file_name=file_name,
        )


def get_alignment_scores() -> dict:
    file_paths: tuple[str, ...] = (
        "../data/alcohol_dehydrogenase-homo_sapiens.fasta",
        "../data/alcohol_dehydrogenase-drosophila_teissieri.fasta",
        "../data/alcohol_dehydrogenase-drosophila_suzukii.fasta",
        "../data/alcohol_dehydrogenase-drosophila_ananassae.fasta",
        "../data/alcohol_dehydrogenase-drosophila_erecta.fasta",
        "../data/alcohol_dehydrogenase-drosophila_yakuba.fasta",
        "../data/alcohol_dehydrogenase-drosophila_sechellia.fasta",
        "../data/alcohol_dehydrogenase-drosophila_simulans.fasta",
        "../data/alcohol_dehydrogenase-drosophila_melanogaster.fasta",
    )

    file_path_combinations: tuple[tuple[str, str]] = tuple(
        combination for combination in itertools.combinations(file_paths, 2)
    )

    alignment_results: AlignmentScores = {}

    for paths in file_path_combinations:
        protein_1: SeqRecord = SeqIO.read(paths[0], 'fasta')
        protein_2: SeqRecord = SeqIO.read(paths[1], 'fasta')

        protein_name_1: str = " ".join(protein_1.description.split(' ')[1:])
        protein_name_2: str = " ".join(protein_2.description.split(' ')[1:])

        alignment_score: float = use_substitution_matrix(protein_1.seq, protein_2.seq)

        alignment_results[protein_name_1, protein_name_2] = alignment_score

    return alignment_results


def main() -> None:

    # download_fasta_and_genbank()
    alignment_scores: AlignmentScores = get_alignment_scores()


if __name__ == "__main__":
    main()
