import pandas as pd
import os
from stats_and_graphs import (top3_stacked_barplot, PCA_all, boxplots_all, plot_clustering, plot_classify, chi2_test)
from ml_tools import clustering, classify

#%% DEFINE GLOBAL VARIABLES
# =========================
DATA_PATH = "./data"

LIST_PHYLUM = ['Actinobacteria',
                'CFB',
                'Proteobacteria',
                'Firmicutes']
                # 'Cyanobacteria',
                # 'Fusobacteria',
                # 'Spirochetes']

LIST_CLASSIFICATION_METHOD = ['Linear',
                            'DecisionTree',
                            'KNeighbors',
                            'RandomForest']

#%% READ DATAFRAME
# ================
dfs_start = []
dfs_stop = []

for phylum_name in LIST_PHYLUM:

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

# #%% TESTS STATS
# # =============
# # Chi2
# p_val_start = chi2_test(dfs_start,LIST_PHYLUM)
# print(f"Chi-square test on start codons: p-value={p_val_start}")

# p_val_stop = chi2_test(dfs_stop,LIST_PHYLUM)
# print(f"Chi-square test on stop codons: p-value={p_val_stop}")

#%% CLUSTERING
# ============
# Start with Kmeans:
kmeans = clustering(dfs_start, n_cluster=len(dfs_start), method="kmeans")
plot_clustering(kmeans, dfs_start, LIST_PHYLUM, "kmeans", "start")
# Start with AgglomerativeHierachiqueClustering:
ac = clustering(dfs_start, n_cluster=len(dfs_start), method="ac")
plot_clustering(ac, dfs_start, LIST_PHYLUM, "agglomerativeclustering", "start")

# # Stop with Kmeans:
# kmeans_2 = clustering(dfs_stop, n_cluster=len(dfs_stop), method="kmeans")
# plot_clustering(kmeans_2, dfs_stop, LIST_PHYLUM, "kmeans", "stop")
# # Stop with AgglomerativeHierachiqueClustering:
# ac_2 = clustering(dfs_stop, n_cluster=len(dfs_stop), method="ac")
# plot_clustering(ac_2, dfs_stop, LIST_PHYLUM, "agglomerativeclustering", "stop")

#%% CLASSIFICATION
# ================
# Start:
for method in LIST_CLASSIFICATION_METHOD:

    print(f"</> {method} Classification")
    
    mat, global_score = classify(dfs_start, method, top3_only=False)
    
    plot_classify(mat, method, LIST_PHYLUM, "start")
    print(f"    Global score: {global_score}")
