import pandas as pd
import os
from stats_and_graphs import (top3_stacked_barplot, PCA_all, boxplots_all, plot_clustering, plot_classify)
from ml_tools import clustering, classify, calculate_tpr_fpr, get_all_roc_coordinates, ROC_curves_and_AOC_scores


#%% GLOBAL VARIABLES
# ==================
DATA_PATH = "./data"

LIST_PHYLUM = ['Actinobacteria',
                'CFB',
                'Proteobacteria',
                'Firmicutes']
                # 'Cyanobacteria',
                # 'Fusobacteria',
                # 'Spirochetes']

LIST_CLUSTERING_METHOD = ['Kmeans', 'AgglomerativeClustering']

LIST_CLASSIFICATION_METHOD = ['Linear',
                            'DecisionTree',
                            'KNeighbors',
                            'RandomForest',
                            'SVM',
                            'NeuralNetwork']

#%% READ DATAFRAME
# ================
dfs_start = []
dfs_stop = []

for phylum_name in LIST_PHYLUM:

    dfs_start.append(pd.read_csv(os.path.join(DATA_PATH, f"{phylum_name.lower()}_start.csv")))
    dfs_stop.append(pd.read_csv(os.path.join(DATA_PATH, f"{phylum_name.lower()}_stop.csv")))

#%% STATISTICAL ANALYSIS
# ======================
# Start:
for df, phylum_name in zip(dfs_start, LIST_PHYLUM):
    top3_stacked_barplot(df, phylum_name, "start")  # Stacked bar graphs

PCA_all(dfs_start, LIST_PHYLUM)  # PCA plot

boxplots_all(dfs_start, LIST_PHYLUM)  # Boxplots top3

# Stop:
for df, phylum_name in zip(dfs_stop, LIST_PHYLUM):
    top3_stacked_barplot(df, phylum_name, "stop")  # Stacked bar graphs

PCA_all(dfs_stop, LIST_PHYLUM)  # PCA plot

boxplots_all(dfs_stop, LIST_PHYLUM)  # Boxplots top3

#%% CLUSTERING
# ============
# Only for tart:
for method in LIST_CLUSTERING_METHOD:
    clust = clustering(dfs_start, n_cluster=len(dfs_start), method=method, top3_only=True)
    
    plot_clustering(clust, dfs_start, LIST_PHYLUM, method, "start")

#%% CLASSIFICATION
# ================
# Only for start:
for method in LIST_CLASSIFICATION_METHOD:

    print(f"</> {method} Classification")
    
    mat, global_score = classify(dfs_start, method, test_size=0.2, top3_only=True)
    
    plot_classify(mat, method, LIST_PHYLUM, "start")
    print(f"    Global score: {global_score}")
    
    
 #%% ROC plot and ROC AUC calculation
# ================

ROC_curves_and_AOC_scores(dfs_start)
