import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import confusion_matrix
from scipy.cluster.hierarchy import dendrogram


#%% TOP3
def top3_stacked_barplot(dataframe, phylum, startoustop):
    """
    Parameters
    ----------
    dataframe : pandas DataFrame
            Output of give_df_codon_pption_per_phylum

    phylum : str,
            phylum name

    startoustop : str
            start or stop according to the nature of the codons
    """
    #On cherche le top3
    sum = dataframe.sum()
    sum = sum[:-1]
    for i in range(len(sum)):
        sum[i]=(float(sum[i])/len(sum))
    df_sum = pd.DataFrame(index=list(sum.index),data=list(sum.values))
    df_sum.columns=['Somme']
    top=df_sum.sort_values(by='Somme', ascending=False)[:3]
    
    #on crée un sous dataframe pour seulement ce top3 avec la catégorie 'others'
    sous_df = dataframe[list(top.index)]
    others = 1 -(sous_df.sum(axis=1))
    sous_df = sous_df.assign(Others=list(others))
    sous_df = sous_df.assign(ID=list(dataframe.index))
    
    #on crée une figure du top3
    fig=sous_df.plot(x = 'ID', kind = 'barh', stacked = True, 
                       title = f"Stacked Bar Graph of {phylum} {startoustop} codons", 
                       mark_right = True, figsize=(8,10), fontsize=10)
    fig.legend(ncol=4,bbox_to_anchor =(0.4,-0.1),loc="lower center")
    plt.yticks([])  # Hide ID number to avoid overlapping lables
    plt.savefig(f"figs/barplot_{phylum}_{startoustop}_codons.png")
    plt.show()

#%% PCA
def df_pour_PCA(dataframe, nom_phylum):
    """
    Parameters
    ----------
    dataframe : panda dataframe, sortie de give_df_codon_pption_per_phylum
    nom_phylum : str, nom du phylum
    Returns
    -------
    df_bis : dataframe avec colonne info nom du phylum
        Data frame des proportions des codons
    """
    phylum=[]
    for i in range(len(dataframe)):
        phylum.append(nom_phylum) 
    df_bis=dataframe.assign(phylum=phylum)
    return df_bis

def PCA_all(dfs, pls):
    # dfs: list of dataframes
    # pls: list of phylum names

    dfs = [df_pour_PCA(df, name) for df, name in zip(dfs, pls)]

    #On fusionne
    fusion=pd.concat(dfs, ignore_index=True)
    fusion = fusion.fillna(0)

    #On fait la PCA
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(fusion.drop('ID', axis=1).drop('phylum',axis=1))

    #On trace le scatter plot
    fig=sns.scatterplot(x=pcs[:,0], y=pcs[:,1], hue=fusion['phylum']) 
    plt.xlabel('PC 1 (%.2f%%)' % (pca.explained_variance_ratio_[0]*100))
    plt.ylabel('PC 2 (%.2f%%)' % (pca.explained_variance_ratio_[1]*100))
    fig.legend(fontsize=8)
    plt.savefig("figs/PCA.png")
    plt.show()

#%% BOXPLOTS
def df_pour_boxplots(dataframe, nom_phylum):
    sum = dataframe.sum()
    sum = sum[:-1]
    for i in range(len(sum)):
        sum[i]=(float(sum[i])/len(sum))
    df_sum = pd.DataFrame(index=list(sum.index),data=list(sum.values))
    df_sum.columns=['Somme']
    top=df_sum.sort_values(by='Somme', ascending=False)[:3]
    df_boxs=dataframe[list(top.index)]
    df_boxs=pd.melt(df_boxs)
    df_boxs = df_pour_PCA(df_boxs, nom_phylum)
    return df_boxs

def boxplots_all(dfs, pls):

    boxs = [df_pour_boxplots(df, name) for df, name in zip(dfs, pls)]

    # On fusionne pour faire le boxplot
    boxs_fusion = pd.concat(boxs)
    boxs_fusion = boxs_fusion.fillna(0)
    fig=sns.boxplot(x='variable', y='value', data=boxs_fusion, hue='phylum')
    plt.xlabel('Top 3 codons')
    plt.ylabel('Proportion (%)')
    fig.legend(fontsize=8)
    plt.savefig(f"figs/condons_proportion_boxplot.png")
    plt.show()

def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)

def plot_clustering(opt, dfs, pls, method, startoustop):
    
    # Confusion Matrix
    target = []
    for i in range(len(dfs)):
        target += [i for j in range(len(dfs[i]))]

    mat = confusion_matrix(target, opt.labels_)

    sns.heatmap(mat.T, square=True, annot=True, fmt='d', cbar=False,
                xticklabels=pls,
                yticklabels=range(len(pls)))

    plt.xlabel('Phylum')
    plt.xticks(rotation=15, ha='right', fontsize=5)
    plt.ylabel('Group')
    plt.title(f"Confusion Matrix of {startoustop} codons using {method} method")
    plt.savefig(f"figs/CM_clustering_{startoustop}_codons_{method}.png")
    plt.show()

    if method == "ac" or method == "agglomerativeclustering":
        # Hierarchical Clustering Dendrogram
        plot_dendrogram(opt, truncate_mode="level", p=3)

        plt.xlabel("Number of points in node (or index of point if no parenthesis).")
        plt.title(f"Hierarchical Clustering Dendrogram of {startoustop} codons using {method} method")
        plt.savefig(f"figs/dendrogram_{startoustop}_codons.png")
        plt.show()

def plot_classify(mat, method, list_phylum, startoustop):

    sns.heatmap(mat.T, square=True, annot=True, fmt='d', cbar=False,
                xticklabels=list_phylum,
                yticklabels=list_phylum)
    plt.xticks(rotation=15, ha='right', fontsize=5)
    plt.xlabel('True label')
    plt.ylabel('Prediction')
    plt.yticks(rotation=15, ha='right', fontsize=5)
    plt.title(f"Confusion Matrix of {startoustop} codons using {method} classification")
    plt.savefig(f"figs/CM_classification_{startoustop}_codons_{method}.png")
    plt.show()

    score = mat.diagonal()/mat.sum(axis=1)
    print(f"</> Scores of prediction based on {startoustop} codon:")

    for s, phylum in zip(score, list_phylum):
        print(f"    Phylum {phylum}: {s}")
