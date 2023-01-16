import pandas as pd
import os
from stats_and_graphs import (top3_stacked_barplot, PCA_all, boxplots_all)

#%% DEFINE GLOBAL VARIABLES
# =========================
DATA_PATH = "./data"

LIST_PHYLUM = ['Actinobacteria',
                'CFB',
                'Cyanobacteria',
                'Firmicutes',
                'Fusobacteria',
                'Proteobacteria',
                'Spirochetes']

#%% READ DATAFRAME
# ================
dfs_start = []
dfs_stop = []

for phylum_name in LIST_PHYLUM:
    #%% save dataframes
    dfs_start.append(pd.read_csv(os.path.join(DATA_PATH, f"{phylum_name.lower()}_start.csv")))
    dfs_stop.append(pd.read_csv(os.path.join(DATA_PATH, f"{phylum_name.lower()}_stop.csv")))

#%% GRAPH AND VISUALIZATION
# =========================
# Start:
for df, phylum_name in zip(dfs_start, LIST_PHYLUM):
    # Stacked bar graphs
    top3_stacked_barplot(df, phylum_name, "start")
# PCA plot
PCA_all(dfs_start, LIST_PHYLUM)
# Boxplots top3
boxplots_all(dfs_start, LIST_PHYLUM)

# Stop:
for df, phylum_name in zip(dfs_stop, LIST_PHYLUM):
    # Stacked bar graphs
    top3_stacked_barplot(df, phylum_name, "stop")
# PCA plot
PCA_all(dfs_stop, LIST_PHYLUM)
# Boxplots top3
boxplots_all(dfs_stop, LIST_PHYLUM)