#%% preprocess
import pandas as pd
from Bio import SeqIO

from rescript_pretraitement import preprocess_gff, preprocess_fasta

fasta_pw_proteobacteria="D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/proteobacteria.fasta"
gff_pw_proteobacteria="D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/proteobacteria.gff3"
gff_proteobacteria = preprocess_gff(gff_pw_proteobacteria)
fasta_proteobacteria = preprocess_fasta(fasta_pw_proteobacteria)

fasta_pw_firmicutes="D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/firmicutes.fasta"
gff_pw_firmicutes="D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/firmicutes.gff3"
gff_firmicutes = preprocess_gff(gff_pw_firmicutes)
fasta_firmicutes = preprocess_fasta(fasta_pw_firmicutes)

# fasta_pw_actinobacteria="D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/actinobacteria.fasta"
# gff_pw_actinobacteria="D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/actinobacteria.gff3"
# gff_actinobacteria = preprocess_gff(gff_pw_actinobacteria)
# fasta_actinobacteria = preprocess_fasta(fasta_pw_actinobacteria)

# fasta_pw_cfb="D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/cfb.fasta"
# gff_pw_cfb="D:/Mes Documents/INSA/5A/UF7 - Défi Bioinfo/Defi-Bioinfo-main/demo_data/cfb.gff3"
# gff_cfb = preprocess_gff(gff_pw_cfb)
# fasta_cfb = preprocess_fasta(fasta_pw_cfb)
#%% extract codon
from rescript_codon_extraction import extract_codon_1bact

start_proteobacteria,stop_proteobacteria = extract_codon_1bact(fasta_proteobacteria, gff_proteobacteria)
start_firmicutes,stop_firmicutes = extract_codon_1bact(fasta_firmicutes, gff_firmicutes)

# start_actinobacteria,stop_actinobacteria = extract_codon_1bact(fasta_actinobacteria, gff_actinobacteria)
# start_cfb,stop_cfb = extract_codon_1bact(fasta_cfb, gff_cfb)
#%% Count codons & Proportions
from collections import Counter
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

#test    
pption_proteobacteria = counter_and_proportion(start_proteobacteria)
pption_firmicutes = counter_and_proportion(start_firmicutes)
#%% DATAFRAME TRANSFORMATION
import numpy as np
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

df_proteobacteria=give_df_codon_pption_per_phylum(pption_proteobacteria)
df_firmicutes=give_df_codon_pption_per_phylum(pption_firmicutes)

# Anciens graphes impossibles à lire à cause des couleurs et trop de codons
# fig_proteobacteria=df_proteobacteria.plot(x = 'ID', kind = 'barh', stacked = True, 
#                     title = 'Stacked Bar Graph of Proteobacteria Start Codons', 
#                     mark_right = True)
# fig_proteobacteria.legend(ncol=8,bbox_to_anchor =(0.4,-0.7),loc="lower center")

# fig_firmicutes=df_firmicutes.plot(x = 'ID', kind = 'barh', stacked = True, 
#                     title = 'Stacked Bar Graph of Firmicutes Start Codons', 
#                     mark_right = True)
# fig_firmicutes.legend(ncol=8,bbox_to_anchor =(0.4,-0.7),loc="lower center")

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
    fig=sous_df.plot(x = 'ID', kind = 'barh', stacked = True, 
                       title = 'Stacked Bar Graph of '+phylum+' '+startoustop +' codons', 
                       mark_right = True)
    fig.legend(ncol=4,bbox_to_anchor =(0.4,-0.2),loc="lower center")

#test
top3_stacked_barplot(df_proteobacteria,"Proteobacteria","start")
top3_stacked_barplot(df_firmicutes,"Frimicutes","start")

# import scipy.stats as stats
# def student_test(proportions_codons1,proportions_codons2,nom_phylum1,nom_phylum2):
#     if sorted(proportions_codons1.keys())==sorted(proportions_codons2.keys()):
#         Student=stats.ttest_ind(list(proportions_codons1.values()),list(proportions_codons2.values()))
#         if Student.pvalue < 0.05:
#             return print("Test de Student entre", nom_phylum1, "et", nom_phylum2,
#                   "=> moyennes non égales", Student.pvalue)
#         else:
#             return print("Test de Student entre", nom_phylum1, "et", nom_phylum2,
#                   "=> moyennes égales", Student.pvalue)
        
# #test
# student_test(counter_and_proportion(start_firmicutes), 
#              counter_and_proportion(start_proteobacteria),
#              "Firmicutes", "Proteobacteria")

#%% PCA
import seaborn as sns
from sklearn.decomposition import PCA

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
    import plotly.express as px
    from sklearn.decomposition import PCA
    
    df = px.data.iris()
    X = df[['sepal_length', 'sepal_width', 'petal_length', 'petal_width']]
    
    pca = PCA(n_components=2)
    components = pca.fit_transform(X)
    
    fig = px.scatter(components, x=0, y=1, color=df['species'])
    fig.show()    
        
        
    df_bis=dataframe.assign(phylum=phylum)
    return df_bis













