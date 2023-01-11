import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import Counter

#Nos fonctions
from rescript_pretraitement import preprocess_gff, preprocess_fasta
from rescript_codon_extraction import extract_codon_1bact
from rescript_mainv2 import (counter_and_proportion, give_df_codon_pption_per_phylum, 
                             top3_stacked_barplot, df_pour_PCA)

# Import des fichiers fasta et gff 
# Deux exemples start : Proteobacteria & Firmicutes
fasta_pw_proteobacteria="D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/proteobacteria.fasta"
gff_pw_proteobacteria="D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/proteobacteria.gff3"
gff_proteobacteria = preprocess_gff(gff_pw_proteobacteria)
fasta_proteobacteria = preprocess_fasta(fasta_pw_proteobacteria)

fasta_pw_firmicutes="D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/firmicutes.fasta"
gff_pw_firmicutes="D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/firmicutes.gff3"
gff_firmicutes = preprocess_gff(gff_pw_firmicutes)
fasta_firmicutes = preprocess_fasta(fasta_pw_firmicutes)

# Stockage des listes de codons start et stop
start_proteobacteria,stop_proteobacteria = extract_codon_1bact(fasta_proteobacteria, gff_proteobacteria)
start_firmicutes,stop_firmicutes = extract_codon_1bact(fasta_firmicutes, gff_firmicutes)

# Transformation en dataframe avec les proportions
df_proteobacteria=give_df_codon_pption_per_phylum(counter_and_proportion(start_proteobacteria))
df_firmicutes=give_df_codon_pption_per_phylum(counter_and_proportion(start_firmicutes))  
                                              
# Tracé des Stacked bar graphs
top3_stacked_barplot(df_proteobacteria,"Proteobacteria","start")
top3_stacked_barplot(df_firmicutes,"Frimicutes","start")

# Scatter plot de la PCA
df1=df_pour_PCA(df_proteobacteria,"Proteobacteria")
df2=df_pour_PCA(df_firmicutes, "Firmicutes")
# On fusionne les 2 dataframes pour faire la PCA
fusion=pd.concat([df1,df2])
# On visualise les 2 phyla avec 2 couleurs différentes sur le scatter plot
import seaborn as sns
from sklearn.decomposition import PCA
pca = PCA(n_components=2)
pcs = pca.fit_transform(fusion.drop('ID', axis=1).drop('phylum',axis=1))
sns.scatterplot(x=pcs[:,0], y=pcs[:,1], hue=fusion['phylum'])

