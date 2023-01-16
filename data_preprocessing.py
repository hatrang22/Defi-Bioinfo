import pandas as pd
from Bio import SeqIO
import os
from codon_extraction import extract_codon

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

def data_preprocessing(raw_data_path, list_phylum, save_dir="data", fasta_ext=".fasta", gff_ext=".gff3"):
    """
    Read fasta/gff files, extract condons, then convert into DataFrame and save results.

    Parameters
    ----------
    raw_data_path : str
        Path to raw data (fasta and gff files).

    list_phylum : list of str
        List of studied phylum names.

    save_dir : str, default 'data'
        Path to save output dataframes.

    fasta_ext : str, default '.fasta'
        Extension of fasta files.

    gff_ext : str, default '.gff3'
        Extension of gff files.
    """

    dfs_start = []
    dfs_stop = []

    for phylum_name in list_phylum:

        phylum_name = phylum_name.lower()

        #%% Import files: fasta and gff
        gff = preprocess_gff(os.path.join(raw_data_path, phylum_name + gff_ext))
        fasta = preprocess_fasta(os.path.join(raw_data_path, phylum_name + fasta_ext))

        #%% Transformation to dataframe within proportions
        dfstart, dfstop = extract_codon(fasta, gff)

        if phylum_name == "fusobacteria":  # drop outlier for fuso (temporal solution)
            dfstart.drop('NG_050724.1', inplace=True) 
        dfs_start.append(dfstart)

        dfs_stop.append(dfstop)

        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        #%% save dataframes
        dfstart.to_csv(os.path.join(save_dir, f"{phylum_name}_start.csv"), index=False)
        dfstop.to_csv(os.path.join(save_dir, f"{phylum_name}_stop.csv"), index=False)


### main preprocessing ###
##########################

if __name__ == '__main__':

    RAW_DATA_PATH = "../raw_data"  # Path to raw data (fasta and gff)

    LIST_PHYLUM = ['Actinobacteria',
                'CFB',
                'Cyanobacteria',
                'Firmicutes',
                'Fusobacteria',
                'Proteobacteria',
                'Spirochetes']  # All studied phylums

    data_preprocessing(RAW_DATA_PATH, LIST_PHYLUM)
