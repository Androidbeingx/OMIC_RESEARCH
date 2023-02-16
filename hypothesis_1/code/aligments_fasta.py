"""
Module for parse aligments of extracted proteins to the research
"""

from Bio            import SeqIO
from Bio.Seq        import Seq
from Bio.SeqRecord  import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.Align      import PairwiseAligner, substitution_matrices, PairwiseAlignments

from Bio.Align.substitution_matrices import Array
from typing import TypedDict

class AlignmentResults(TypedDict):
    alignment: PairwiseAlignments
    score: float
    percent_identity: float


def read_fasta(filepath: str) -> str:
    """ 
    Reads a fasta file
    param: filepath,
    returns a sequence string
    """

    
    record: SeqRecord = SeqIO.read(filepath, 'fasta')

    return record.seq
  


# Substitution Matrix (BLOSUM62, etc.) for aligning proteins.
# BLAST uses BLOSUM62 in blastp.
# ---------------------------------------------------------------------
def use_blosum62_matrix(prot1: str, prot2: str) -> AlignmentResults:
    """
    Calculates the score of aligment between two proteins
    Param: two strings of proteins,
    Return int with the score
    """

    # Create Aligner
    aligner: PairwiseAligner = PairwiseAligner()

    # Get one of the matrices
    blosum62_matrix: list[str] | Array = substitution_matrices.load('BLOSUM62')

    # Put the matrix in the aligner.
    # This invalidates the match, mismatch and gap scores.
    aligner.substitution_matrix = blosum62_matrix

    alignment: PairwiseAlignments = aligner.align(prot1, prot2)

    score: float = alignment.score

    # Calc. the max score can be achived by the this alignment.
    # (Align the first seq. with itself.)
    max_score: float = aligner.score(prot1,prot1)

    # Calc. the alignment score as a percentage.
    percent_identity: float = (alignment.score / max_score) * 100

    alignment_results: AlignmentResults = {
            'alignment': alignment,
            'score': score,
            'percent_identity': percent_identity
            }

    return alignment_results

