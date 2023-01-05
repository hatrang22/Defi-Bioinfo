def extract_codon_phylum(dseq: dict):
    """
    Parameters
    ----------
    dseq: dict
            Dictionary consisting of phylum infomation (accession number, fasta and gff).

    Returns
    -------
    ret: dict
        Dictionary with 3 keys:

            - accession_number: list of accession numbers
            - start: list of corresponding start codons
            - stop: list of corresponding stop codons
    """
    
    ret = {"accession_number": [], "start": [], "stop": []}

    for i, acc_num in enumerate(dseq["accession_number"]):
        
        ret["accession_number"].append(acc_num)

        start, stop = extract_codon_1bact(dseq["fasta"][i], dseq["gff"][i])

        ret["start"].append(start)
        ret["stop"].append(stop)

    return ret

def extract_codon_1bact(fasta: str, gff: list):
    #%%On extrait les codons start et stop
    start_plus,stop_plus,start_moins,stop_moins = ([] for i in range(4))
    for i in gff:
        if i[2]=="CDS":
            if i[6]=="+":
                start_plus.append(fasta[int(i[3])-1 : int(i[3])+2]) #i[3] correspond à l'indice du début de la CDS, on extrait les premiers éléments de la CDS (du fasta, attention pas la même indexation entre python et fasta)
                stop_plus.append(fasta[int(i[4])-3 : int(i[4])])  #i[4] correspond à l'indice de fin d'une CDS
                
            if i[6]=="-": #attention au sens de lecture du brin '-' !!
                start_moins.append(fasta[int(i[3])-1 : int(i[3])+2])
                stop_moins.append(fasta[int(i[4])-3 : int(i[4])])

    rev_comp_st(stop_moins)
    rev_comp_st(start_moins)

    #%% On fusionne les codons strand + et -
    start = start_plus + stop_moins
    stop = stop_plus + start_moins

    return start, stop

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
