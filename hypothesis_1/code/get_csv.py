"""
Module to make a table and csv with all the data
"""
from aligments_fasta import *
import pandas as pd
from pathlib import Path  

drosophila_ananassae       = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_ananassae.fasta')
drosophila_erecta          = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_erecta.fasta')
drosophila_melanogaster    = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_melanogaster.fasta')
drosophila_sechellia       = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_sechellia.fasta')
drosophila_simulans        = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_simulans.fasta')
drosophila_suzukii         = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_suzukii.fasta')
drosophila_teissieri       = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_teissieri.fasta')
drosophila_yakuba          = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_yakuba.fasta')
drosophila_guanche         = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_guanche.fasta')
drosophila_persimilis      = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_persimilis.fasta')
homo_sapiens              = read_fasta('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-homo_sapiens.fasta')

#GET ALIGMENTS
aligments = []
aligments.append(use_blosum62_matrix(homo_sapiens, drosophila_ananassae )) 
aligments.append(use_blosum62_matrix(homo_sapiens, drosophila_erecta)) 
aligments.append(use_blosum62_matrix(homo_sapiens, drosophila_melanogaster)) 
aligments.append(use_blosum62_matrix(homo_sapiens, drosophila_sechellia)) 
aligments.append(use_blosum62_matrix(homo_sapiens, drosophila_simulans)) 
aligments.append(use_blosum62_matrix(homo_sapiens, drosophila_suzukii)) 
aligments.append(use_blosum62_matrix(homo_sapiens, drosophila_teissieri)) 
aligments.append(use_blosum62_matrix(homo_sapiens, drosophila_yakuba ))
aligments.append(use_blosum62_matrix(homo_sapiens, drosophila_guanche )) 
aligments.append(use_blosum62_matrix(homo_sapiens, drosophila_persimilis  ))  


# GET ACCESSION AND DEFINITION
drosophila_ananassae_gb_path    = ('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_ananassae.gb')
drosophila_erecta_gb_path       = ('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_erecta.gb')
drosophila_melanogaster_gb_path = ('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_melanogaster.gb')
drosophila_sechellia_gb_path    = ('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_sechellia.gb')
drosophila_simulans_gb_path     = ('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_simulans.gb')
drosophila_suzukii_gb_path      = ('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_suzukii.gb')
drosophila_teissieri_gb_path    = ('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_teissieri.gb')
drosophila_yakuba_gb_path       = ('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_yakuba.gb')
drosophila_guanche_gb_path      = ('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_guanche.gb')
drosophila_persimilis_gb_path   = ('dawbio2_m14_uf2_pt1/hypothesis_1/data/genbanks-and-fastas/alcohol_dehydrogenase-drosophila_persimilis.gb')


def get_acc_def (gb_path: str) -> dict :
    """
    From SeqRecord gets id and description
    Param gb filepath
    Return dict with data
    """
    record: SeqRecord = SeqIO.read(gb_path, 'genbank')
    return {"accesion_num": record.id, "description_seq": record.description}

gb_paths = [drosophila_ananassae_gb_path, drosophila_erecta_gb_path, drosophila_melanogaster_gb_path, 
            drosophila_sechellia_gb_path, drosophila_simulans_gb_path, drosophila_suzukii_gb_path,  drosophila_teissieri_gb_path,
            drosophila_yakuba_gb_path, drosophila_guanche_gb_path, drosophila_persimilis_gb_path]
id_acc = []

for path in gb_paths:
    id_acc.append(get_acc_def(path))


#SET DATFRAME
df=pd.DataFrame(aligments)
names = [ "drosphila_ananassae","drosphila_erecta", "drosphila_melanogaster", "drosphila_sechellia", 
        "drosphila_simulans", "drosphila_suzukii", "drosphila_teissieri", "drosphila_yakuba", "drosphila_guanche ", "drosphila_persimilis" ]
df.insert(0, 'specie', names)



df1 = pd.DataFrame(id_acc)

result= pd.concat([df, df1], axis=1)

#falta las regex de ani
database = result.loc[:,["specie", "accesion_num", "description_seq","score", "percent_identity"]].sort_values(by="score", ascending=False).reset_index(drop=True)

#TO CSV RESULTS
filepath = Path('dawbio2_m14_uf2_pt1/hypothesis_1/data/csv/results.csv')  
filepath.parent.mkdir(parents=True, exist_ok=True)  
database.to_csv(filepath, sep=";")
