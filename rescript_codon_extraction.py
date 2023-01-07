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
def extract_codon_1bact(fasta: dict, gff):
    #On extrait les codons start et stop
    start_plus,stop_plus,start_moins,stop_moins = ([] for i in range(4))
    
    for i in gff.index:
        for cle,seq in fasta.items() :
            if gff.loc[i,'ID']==cle:
                if gff.loc[i,'strand']=="+":
                    start_plus.append(seq[int(gff.loc[i,'start'])-1 : int(gff.loc[i,'start'])+2]) #i[3] correspond à l'indice du début de la CSD, on extrait les premiers éléments de la CSD (du fasta, attention pas la même indexation entre python et fasta)
                    stop_plus.append(seq[int(gff.loc[i,'stop'])-3 : int(gff.loc[i,'stop'])])  #i[4] correspond à l'indice de fin d'une CSD
                    
                if gff.loc[i,'strand']=="-": #attention au sens de lecture du brin '-' !!
                    start_moins.append(seq[int(gff.loc[i,'start'])-1 : int(gff.loc[i,'start'])+2])
                    stop_moins.append(seq[int(gff.loc[i,'stop'])-3 : int(gff.loc[i,'stop'])])

    rev_comp_st(stop_moins)
    rev_comp_st(start_moins)

    # On fusionne les codons strand + et -
    start = start_plus + stop_moins
    stop = stop_plus + start_moins

    return start, stop
