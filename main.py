from data_preprocessing import data_preprocessing
from codon_extraction import extract_codon
from stats_and_graphs import (top3_stacked_barplot, PCA_all, boxplots_all)
import os

#%% Define global variables
# =========================
DATA_PATH = "../data"

FASTA_EXTENSION = ".fasta"
GFF_EXTENSION = ".gff3"
###

#%% Preprocessing data and transformation to dataframe
# ====================================================
list_phylum = [f.split(".")[0] for f in os.listdir(DATA_PATH) if f.endswith("fasta")]  # get list of phylum

dfs_start = []
dfs_stop = []

for phylum_name in list_phylum:
    #%% Import files: fasta and gff
    fasta, gff = data_preprocessing(DATA_PATH, phylum_name, FASTA_EXTENSION, GFF_EXTENSION)

    #%% Transformation to dataframe within proportions
    dfstart, dfstop = extract_codon(fasta, gff)

    if phylum_name == "fusobacteria":  # drop outlier for fuso (temporal solution)
        dfstart.drop('NG_050724.1', inplace=True) 
    dfs_start.append(dfstart)

    dfs_stop.append(dfstop)

    #%% save dataframes
    dfstart.to_csv(f"{phylum_name}_start.csv", index=False)
    dfstop.to_csv(f"{phylum_name}_stop.csv", index=False)

#%% Graph and visualization
# =========================
# Start:
for df, phylum_name in zip(dfs_start, list_phylum):
    # Stacked bar graphs
    top3_stacked_barplot(df, phylum_name, "start")
# PCA plot
PCA_all(dfs_start, list_phylum)
# Boxplots top3
boxplots_all(dfs_start, list_phylum)

# Stop:
for df, phylum_name in zip(dfs_stop, list_phylum):
    # Stacked bar graphs
    top3_stacked_barplot(df, phylum_name, "stop")
# PCA plot
PCA_all(dfs_stop, list_phylum)
# Boxplots top3
boxplots_all(dfs_stop, list_phylum)