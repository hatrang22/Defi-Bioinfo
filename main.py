from data_preprocessing import read_phylum
from extract_codon import extract_codon_phylum
from plot import plot_codon_repartition

####################
# DEFINE PHYLUM
####################
fp = "demo_data/actinobacteria10_genome_fasta.tar"
gp = "demo_data/actinobacteria10_genome_gff.tar"

####################
# PREPROCESSING DATA
####################
phylum = read_phylum(fp, gp)

####################
# EXTRACT CODON
####################
codon_phylum = extract_codon_phylum(phylum)

####################
# ANALYSE
####################
# TODO
# Example: plot codon repartition for the first bacterie
plot_codon_repartition(codon_phylum["accession_number"][0], codon_phylum["start"][0], codon_phylum["stop"][0])