from io import StringIO
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _Matrix
from Bio import Phylo

# Input your protein sequences in FASTA format
protein_sequences = """>seq111
MRYKAAA
>seq222
MRYKAKA
>seq333
MKYKAQA
>seq444
MKYKAKA"""

# Create an alignment object from the protein sequences
alignment = AlignIO.read(StringIO(protein_sequences), "fasta")

# Create a distance matrix using the alignment
calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignment)

# Construct the phylogenetic tree using the Neighbor Joining algorithm
constructor = DistanceTreeConstructor()

# Build tree.
tree = constructor.nj(distance_matrix)

# Format tree to a "phylo"xml format.
# If you want to use it on the web with D3.js,
# then you can convert the obtained xml into a json format, with another tool
# ------------------------------
# print(tree.format('phyloxml'))
# ------------------------------

# Draw the tree using the Phylo library
Phylo.draw_ascii(tree)
