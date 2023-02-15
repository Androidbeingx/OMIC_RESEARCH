"""
Module for parse aligments of extracted proteins to the research
"""

from Bio            import SeqIO
from Bio.Seq        import Seq
from Bio.SeqRecord  import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.Align      import PairwiseAligner, substitution_matrices

from Bio.Align.substitution_matrices import Array





def read_fasta(filepath: str) -> str:
    """ 
    Reads a fasta file
    param: filepath,
    returns a sequence string
    """

    
    record: SeqRecord = SeqIO.read(filepath, 'fasta')

    return record.seq
  

# in the main module as dict

# drosphila_ananassae       = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/alcohol_dehydrogenase-drosophila_ananassae.fasta')
# drosphila_erecta          = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/alcohol_dehydrogenase-drosophila_erecta.fasta')
# drosphila_melanogaster    = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/alcohol_dehydrogenase-drosophila_melanogaster.fasta')
# drosphila_sechellia       = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/alcohol_dehydrogenase-drosophila_sechellia.fasta')
# drosphila_simulans        = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/alcohol_dehydrogenase-drosophila_simulans.fasta')
# drosphila_suzukii         = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/alcohol_dehydrogenase-drosophila_suzukii.fasta')
# drosphila_teissieri       = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/alcohol_dehydrogenase-drosophila_teissieri.fasta')
# drosphila_yakuba          = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/alcohol_dehydrogenase-drosophila_yakuba.fasta')
# homo_sapiens              = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/alcohol_dehydrogenase-homo_sapiens.fasta')


# Substitution Matrix (BLOSUM62, etc.) for aligning proteins.
# BLAST uses BLOSUM62 in blastp.
# ---------------------------------------------------------------------
def use_substitution_matrix(prot1: str, prot2: str) -> float:
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

    # Align using the matrix
    score: float = aligner.score(prot1,prot2)
    return score

