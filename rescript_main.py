#%% preprocess
import pandas as pd
from Bio import SeqIO

from rescript_pretraitement import preprocess_gff, preprocess_fasta

fasta_pw="C:/Users/Clementine/OneDrive/Documents/ENSAT/3a/UF7/Projet Info/Phylum nucleotide way/actino.fasta"
gff_pw="C:/Users/Clementine/OneDrive/Documents/ENSAT/3a/UF7/Projet Info/Phylum nucleotide way/actino.gff3"

df_gff = preprocess_gff(gff_pw)
fasta = preprocess_fasta(fasta_pw)

#%% extract codon
from rescript_codon_extraction import extract_codon_1bact

start,stop = extract_codon_1bact(fasta, df_gff)

#%%Count codon
from collections import Counter

count_start=Counter(start)
count_stop=Counter(stop)