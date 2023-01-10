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
def counter_and_proportion(result_codons,startoustop):
    """
    Argument startoustop : indicate a string "start" or "stop"
    """
    d_start={}
    d_stop={}
    for cle, liste in result_codons.items():
        
        if startoustop == "start":
            count_start=Counter(liste)
            pption_start=count_start
            for i in pption_start:
                pption_start[i]=count_start[i]/len(liste)
            d_start[cle]=pption_start          
        if startoustop == "stop":
            count_stop=Counter(liste)
            pption_stop=count_stop
            for i in count_stop:
                pption_stop[i]=count_stop[i]/len(result_codons)
            d_stop[cle]=pption_stop
        
    if startoustop == "start":
        df=d_start
    if startoustop == "stop":
        df=d_stop
    return df

#test    
pption_proteobacteria = counter_and_proportion(start_proteobacteria, "start")

#%% Transform df proportions into data.frame of codons proportions
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def give_df_codon_pption_per_phylum(dictionnary_propotions:dict):
    """
    Parameters
    ----------
    all_codons_phylum : list
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
        
    #on produit panda data frame des proportions de chaque codons dans un phylum
    df = pd.DataFrame(index=[tous_codons])
    for organism,counter in dictionnary_propotions.items(): #pour chaque organisme
        for key,value in counter.items(): #pour chaque codon (key)
                df.loc[key,organism]=value
    
    #on remplace les nan par des valeurs 0          
    df = df.fillna(0)
    return df

dataframe= give_df_codon_pption_per_phylum(pption_proteobacteria)

#%%STACKED BARPLOTS

dataframe.plot(
    x = 'index',
    kind = 'barh',
    stacked = True,
    title = 'Stacked Bar Graph',
    mark_right = True)

# Ci-dessous un tracé de barplots non empilés, il y en a un par espèce. Mais tout ça ne nous sert à rien!
# for cle, counter in pption_proteobacteria.items():
#     y_pos = range(len(list(counter.keys())))
#     plt.bar(y_pos, sorted(list(counter.values())))
#     plt.xticks(y_pos, sorted(list(counter.keys())),fontsize=5,rotation=50)
#     plt.show()

#%% TEST SHAPIRO
import numpy as np
from scipy.stats import shapiro
def shapiro_test(proportions_codons,nom_phylum):
    Shapiro=shapiro(list(proportions_codons.values()))
    if Shapiro.pvalue < 0.05:
        return print("Test de Shapiro de", nom_phylum, "=> distribution non normale", Shapiro.pvalue)
    else:
        return print("Test de Shapiro de", nom_phylum, "=> distribution normale", Shapiro.pvalue)

#tests
shapiro_test(counter_and_proportion(start_proteobacteria, "start"),"Proteobacteria")
shapiro_test(counter_and_proportion(start_firmicutes, "start"),"Firmicutes")

#%% TEST STUDENT 
import scipy.stats as stats
def student_test(proportions_codons1,proportions_codons2,nom_phylum1,nom_phylum2):
    if sorted(proportions_codons1.keys())==sorted(proportions_codons2.keys()):
        Student=stats.ttest_ind(list(proportions_codons1.values()),list(proportions_codons2.values()))
        if Student.pvalue < 0.05:
            return print("Test de Student entre", nom_phylum1, "et", nom_phylum2,
                  "=> moyennes non égales", Student.pvalue)
        else:
            return print("Test de Student entre", nom_phylum1, "et", nom_phylum2,
                  "=> moyennes égales", Student.pvalue)
        
#test
student_test(counter_and_proportion(start_firmicutes, "start"), 
             counter_and_proportion(start_proteobacteria, "start"),
             "Firmicutes", "Proteobacteria")

#%% PCA

    
    
    
