'''
NOT FINISHED YET!
'''

from Bio import SeqIO
from Bio.Seq import Seq

# Read the FASTA file
records = list(SeqIO.parse("./output/adh-test.fasta", "fasta"))

# Translate the sequences to protein
for record in records:
    protein_seq = record.translate()
    print(protein_seq)
