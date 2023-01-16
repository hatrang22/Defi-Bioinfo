from collections import Counter
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

#%% Count codons & Proportions
def counter_and_proportion(result_codons):
    """
    Parameters
    ----------
    result_codons : dict, sortie de extract_codon_1bact
    Returns
    -------
    d : dict, sous la forme de counter avec les proportions calculées
    """
    d={}
    for cle, liste in result_codons.items():
            count=Counter(liste)
            pption=count
            for i in pption:
                pption[i]=count[i]/len(liste)
            d[cle]=pption          
    return d

#%% DATAFRAME TRANSFORMATION
def give_df_codon_pption_per_phylum(dictionnary_propotions:dict):
    """
    Parameters
    ----------
    dictionnary_propotions : dict
    Returns
    -------
    df : panda dataframe
        Data frame des proportions des codons
    """
    #on extrait liste de tous les codons du phylum
    tous_codons=[]
    for counter in dictionnary_propotions.values():
        tous_codons=tous_codons+(list(counter.keys())) 
    unique=np.unique(np.array(tous_codons))
    #on produit panda data frame des proportions de chaque codons dans un phylum
    df = pd.DataFrame(index=unique)
    for organism,counter in dictionnary_propotions.items(): #pour chaque organisme
        for key,value in counter.items(): #pour chaque codon (key)
                df.loc[key,organism]=value
    
    #on remplace les nan par des valeurs 0          
    df = df.fillna(0)
    df=df.T
    df=df.assign(ID=list(df.index))
    return df

#%% TOP3
def top3_stacked_barplot(dataframe, phylum, startoustop):
    """
    Parameters
    ----------
    dataframe : panda dataframe, sortie de give_df_codon_pption_per_phylum
    phylum : str, nom du phylum
    startoustop : str, start ou stop selon la nature des codons
    
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
    #plt.figure(figsize=(300,150))
    fig=sous_df.plot(x = 'ID', kind = 'barh', stacked = True, 
                       title = 'Stacked Bar Graph of '+phylum+' '+startoustop +' codons', 
                       mark_right = True, figsize=(8,10), fontsize=10)
    fig.legend(ncol=4,bbox_to_anchor =(0.4,-0.1),loc="lower center")
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

def PCA_tous(dfs, pls):
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
    fig.legend(ncol=4,bbox_to_anchor =(0.45,-0.35),loc="lower center")
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

def boxplots_tous(dfs, pls):

    boxs = [df_pour_boxplots(df, name) for df, name in zip(dfs, pls)]

    # On fusionne pour faire le boxplot
    boxs_fusion = pd.concat(boxs)
    boxs_fusion = boxs_fusion.fillna(0)
    fig=sns.boxplot(x='variable', y='value', data=boxs_fusion, hue='phylum')
    plt.xlabel('Top 3 codons')
    plt.ylabel('Proportion (%)')
    fig.legend(ncol=4,bbox_to_anchor =(0.45,-0.35),loc="lower center")
    plt.show()