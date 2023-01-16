import pandas as pd
from Bio import SeqIO

#%% df_gff
def preprocess_gff(gff_datapath) :

    print(f"</> Preprocessing gff at {gff_datapath}")

    df_gff = pd.read_csv(gff_datapath, sep='\t', skiprows = 5, header= None,
                         usecols = list(range(0,8)),
                         names = ['ID', 'RefSeq', 'type', 'start', 'stop', 'dot1', 'strand', 'dot2']
                         )
                         
    df_gff = df_gff[df_gff['type'] == 'CDS'] # filter 'CDS'

    return(df_gff)

#%% df_fasta
def preprocess_fasta(fasta_datapath):

    print(f"</> Preprocessing fasta at {fasta_datapath}")

    fasta_dictio = {}

    fasta_sequences = SeqIO.parse(open(fasta_datapath),'fasta')

    for fasta in fasta_sequences:

        identity, sequence, description = str(fasta.id), str(fasta.seq), str(fasta.description)

        if 'plasmid' in description :
            continue

        else :
            fasta_dictio[identity] = sequence
            
    return fasta_dictio