import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from collections import Counter
from sklearn.decomposition import PCA
#Nos fonctions
from rescript_pretraitement import preprocess_gff, preprocess_fasta
from rescript_codon_extraction import extract_codon_1bact
from stats_and_graphs import (counter_and_proportion, give_df_codon_pption_per_phylum, 
                             top3_stacked_barplot, df_pour_PCA, PCA_tous, df_pour_boxplots, boxplots_tous)

#%% Import des fichiers fasta et gff 
gff_proteobacteria = preprocess_gff("D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/proteobacteria.gff3")
fasta_proteobacteria = preprocess_fasta("D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/proteobacteria.fasta")

gff_firmicutes = preprocess_gff("D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/firmicutes.gff3")
fasta_firmicutes = preprocess_fasta("D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/firmicutes.fasta")

gff_actinobacteria = preprocess_gff("D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/actinobacteria.gff3")
fasta_actinobacteria = preprocess_fasta("D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/actinobacteria.fasta")

gff_cfb = preprocess_gff("D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/cfb.gff3")
fasta_cfb = preprocess_fasta("D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/cfb.fasta")

gff_fusobacteria = preprocess_gff("D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/fusobacteria.gff3")
fasta_fusobacteria = preprocess_fasta("D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/fusobacteria.fasta")

gff_spirochetes = preprocess_gff("D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/spirochetes.gff3")
fasta_spirochetes = preprocess_fasta("D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/spirochetes.fasta")

gff_cyanobacteria = preprocess_gff("D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/cyanobacteria.gff3")
fasta_cyanobacteria = preprocess_fasta("D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/cyanobacteria.fasta")

#%% Stockage des listes de codons start et stop
start_proteobacteria,stop_proteobacteria = extract_codon_1bact(fasta_proteobacteria, gff_proteobacteria)
start_firmicutes,stop_firmicutes = extract_codon_1bact(fasta_firmicutes, gff_firmicutes)
start_actinobacteria,stop_actinobacteria = extract_codon_1bact(fasta_actinobacteria, gff_actinobacteria)
start_cfb,stop_cfb = extract_codon_1bact(fasta_cfb, gff_cfb)
start_fusobacteria,stop_fusobacteria = extract_codon_1bact(fasta_fusobacteria, gff_fusobacteria)
start_spirochetes,stop_spirochetes = extract_codon_1bact(fasta_spirochetes, gff_spirochetes)
start_cyanobacteria,stop_cyanobacteria = extract_codon_1bact(fasta_cyanobacteria, gff_cyanobacteria)

#%% START : Transformation en dataframe avec les proportions
df_proteobacteria=give_df_codon_pption_per_phylum(counter_and_proportion(start_proteobacteria))
df_firmicutes=give_df_codon_pption_per_phylum(counter_and_proportion(start_firmicutes))  
df_actinobacteria=give_df_codon_pption_per_phylum(counter_and_proportion(start_actinobacteria))  
df_cfb=give_df_codon_pption_per_phylum(counter_and_proportion(start_cfb)) 
df_fusobacteria=give_df_codon_pption_per_phylum(counter_and_proportion(start_fusobacteria)) 
df_spirochetes=give_df_codon_pption_per_phylum(counter_and_proportion(start_spirochetes)) 
df_cyanobacteria=give_df_codon_pption_per_phylum(counter_and_proportion(start_cyanobacteria))

# GRAPHES
# Tracé des Stacked bar graphs
top3_stacked_barplot(df_proteobacteria,"Proteobacteria","start")
top3_stacked_barplot(df_firmicutes,"Frimicutes","start")
top3_stacked_barplot(df_actinobacteria,"Actinobacteria","start")
top3_stacked_barplot(df_cfb,"DFB","start")
top3_stacked_barplot(df_fusobacteria,"Fusobacteria","start")
top3_stacked_barplot(df_spirochetes,"Spirochetes","start")
top3_stacked_barplot(df_cyanobacteria,"Cyanobacteria","start")

# Scatter plot de la PCA
#PCA_tous(df_proteobacteria,df_firmicutes,df_actinobacteria)
df_fuso=df_fusobacteria.drop('NG_050724.1')
PCA_tous(df_proteobacteria,df_firmicutes,df_actinobacteria,df_cfb,df_fuso,df_spirochetes,df_cyanobacteria)

# Boxplots  pour le top3
#boxplots_tous(df_proteobacteria, df_firmicutes, df_actinobacteria)
boxplots_tous(df_proteobacteria,df_firmicutes,df_actinobacteria,df_cfb,df_fuso,df_spirochetes,df_cyanobacteria)

#%% STOP : Transformation en dataframe avec les proportions
dfP_proteobacteria=give_df_codon_pption_per_phylum(counter_and_proportion(stop_proteobacteria))
dfP_firmicutes=give_df_codon_pption_per_phylum(counter_and_proportion(stop_firmicutes))  
dfP_actinobacteria=give_df_codon_pption_per_phylum(counter_and_proportion(stop_actinobacteria))  
dfP_cfb=give_df_codon_pption_per_phylum(counter_and_proportion(stop_cfb)) 
dfP_fusobacteria=give_df_codon_pption_per_phylum(counter_and_proportion(stop_fusobacteria)) 
dfP_spirochetes=give_df_codon_pption_per_phylum(counter_and_proportion(stop_spirochetes)) 
dfP_cyanobacteria=give_df_codon_pption_per_phylum(counter_and_proportion(stop_cyanobacteria))
                                             
# GRAPHES
# Tracé des Stacked bar graphs
top3_stacked_barplot(dfP_proteobacteria,"Proteobacteria","stop")
top3_stacked_barplot(dfP_firmicutes,"Frimicutes","stop")
top3_stacked_barplot(dfP_actinobacteria,"Actinobacteria","stop")
top3_stacked_barplot(dfP_cfb,"DFB","stop")
top3_stacked_barplot(dfP_fusobacteria,"Fusobacteria","stop")
top3_stacked_barplot(dfP_spirochetes,"Spirochetes","stop")
top3_stacked_barplot(dfP_cyanobacteria,"Cyanobacteria","stop")

# Scatter plot de la PCA
#PCA_tous(df_proteobacteria,df_firmicutes,df_actinobacteria)
PCA_tous(dfP_proteobacteria,dfP_firmicutes,dfP_actinobacteria,dfP_cfb,dfP_fusobacteria,dfP_spirochetes,dfP_cyanobacteria)

# Boxplots  pour le top3
#boxplots_tous(df_proteobacteria, df_firmicutes, df_actinobacteria)
boxplots_tous(dfP_proteobacteria,dfP_firmicutes,dfP_actinobacteria,dfP_cfb,dfP_fusobacteria,dfP_spirochetes,dfP_cyanobacteria)




    



