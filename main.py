from data_preprocessing import preprocess_gff, preprocess_fasta
from codon_extraction import extract_codon_1bact
from stats_and_graphs import (counter_and_proportion, give_df_codon_pption_per_phylum, 
                             top3_stacked_barplot, PCA_tous, boxplots_tous)
import os

#%% Define global variables
# =========================
DATA_PATH = "demo_data"

FASTA_EXTENSION = ".fasta"
GFF_EXTENSION = ".gff3"
###

#%% Preprocessing data and transformation to dataframe
# ====================================================
list_phylum = [f.split(".")[0] for f in os.listdir(DATA_PATH) if f.endswith("fasta")]  # get list of phylum

dfs_start = []
dfs_stop = []

for phylum in list_phylum:
    #%% Import files: fasta and gff
    gff = preprocess_gff(os.path.join(DATA_PATH, phylum + GFF_EXTENSION))
    fasta = preprocess_fasta(os.path.join(DATA_PATH, phylum + FASTA_EXTENSION))

    #%% get start and stop codons
    start, stop = extract_codon_1bact(fasta, gff)

    #%% START : Transformation to dataframe within proportions
    dfstart = give_df_codon_pption_per_phylum(counter_and_proportion(start))
    if phylum == "fusobacteria":  # drop outlier for fuso (temporal solution)
        dfstart.drop('NG_050724.1', inplace=True) 
    dfs_start.append(dfstart)

    #%% STOP : Transformation to dataframe within proportions
    dfs_stop.append(give_df_codon_pption_per_phylum(counter_and_proportion(stop)))

#%% Graph and visualization
# =========================
# Start:
for df, phylum_name in zip(dfs_start, phylum):
    # Stacked bar graphs
    top3_stacked_barplot(df, phylum_name, "start")
# PCA plot
PCA_tous(dfs_start, phylum)
# Boxplots top3
boxplots_tous(dfs_start, phylum)

# Stop:
for df, phylum_name in zip(dfs_stop, phylum):
    # Stacked bar graphs
    top3_stacked_barplot(df, phylum_name, "stop")
# PCA plot
PCA_tous(dfs_stop, phylum)
# Boxplots top3
boxplots_tous(dfs_stop, phylum)