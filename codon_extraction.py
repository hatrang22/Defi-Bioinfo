from tqdm import tqdm
import pandas as pd
from collections import Counter
import numpy as np

#%%rev comp
def rev_comp_st(seq):
    """
    This function returns a reverse complement
    of a DNA strand
    """
    for i in range(len(seq)):

        # complement strand
        seq[i] = seq[i].replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
        seq[i] = seq[i].upper()
         
        # reverse strand
        seq[i] = seq[i][::-1]
        
#%% extract codon
def extract_codon_1bact(fasta, gff):

    d_start={} #dictionnaire avec ID, stop et start
    d_stop={}

    for cle, seq in tqdm(fasta.items(), desc="    Extracting codons"):

        # On extrait les codons start et stop pour un ID
        start_plus,stop_plus,start_moins,stop_moins = ([] for i in range(4))
        for i in gff.index:
            if gff.loc[i,'ID']==cle:
                if gff.loc[i,'strand']=="+":
                    start_plus.append(seq[int(gff.loc[i,'start'])-1 : int(gff.loc[i,'start'])+2]) #i[3] correspond à l'indice du début de la CSD, on extrait les premiers éléments de la CSD (du fasta, attention pas la même indexation entre python et fasta)
                    
                    if int(gff.loc[i,'stop']) > len(seq)-1:
                        continue
                    else:
                        stop_plus.append(seq[int(gff.loc[i,'stop'])-3 : int(gff.loc[i,'stop'])])  #i[4] correspond à l'indice de fin d'une CSD
                    # if seq[int(gff.loc[i,'stop'])-3 : int(gff.loc[i,'stop'])] == '':
                    #     print(i, "Stop+")
                        
                if gff.loc[i,'strand']=="-": #attention au sens de lecture du brin '-' !!
                    start_moins.append(seq[int(gff.loc[i,'start'])-1 : int(gff.loc[i,'start'])+2])
                    
                    if int(gff.loc[i,'stop']) > len(seq)-1:
                        continue
                    else:
                        stop_moins.append(seq[int(gff.loc[i,'stop'])-3 : int(gff.loc[i,'stop'])])
                    # if seq[int(gff.loc[i,'stop'])-3 : int(gff.loc[i,'stop'])] == '':
                    #     print(i, "Stop-")

        rev_comp_st(stop_moins)
        rev_comp_st(start_moins)
        
        # On fusionne les codons strand + et -
        start = start_plus + stop_moins
        stop = stop_plus + start_moins
            
        d_start[cle]=start
        d_stop[cle]=stop

    return d_start,d_stop

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

def extract_codon(fasta, gff):
    #%% get start and stop codons
    start, stop = extract_codon_1bact(fasta, gff)

    dfstart = give_df_codon_pption_per_phylum(counter_and_proportion(start))
    dfstop = give_df_codon_pption_per_phylum(counter_and_proportion(stop))

    return dfstart, dfstop
