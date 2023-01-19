# Statistical exploration of codon usage

This project is part of "Défi-Bioinformatique" at INSA/ENSAT Toulouse, France, consisting of statistical exploration of start and stop codons in prokaryotes. The objective is to know the distribution of these codons according to the bacterial phyla and to allow the development of a prediction tool that, from an unknown genome, would assign the corresponding phylum based on the proportions of its start and stop codons.

Firstly, the data collected from [NCBI](https://www.ncbi.nlm.nih.gov/nuccore) are preprocessed using the following command, that allows us to read the fasta/gff files, extract codons and convert raw data into DataFrame (assume that this process has already been done in this git and the generated DataFrames are available at `./data`).
```
$ python data_preprocessing.py
</> Preprocessing gff at ../raw_data\actinobacteria.gff3
</> Preprocessing fasta at ../raw_data\actinobacteria.fasta
    Extracting codons: 100%|███████████████████████████████████████████████████████████| 75/75 [03:13<00:00,  2.58s/it]
</> Preprocessing gff at ../raw_data\cfb.gff3
</> Preprocessing fasta at ../raw_data\cfb.fasta
    Extracting codons: 100%|███████████████████████████████████████████████████████████| 75/75 [02:45<00:00,  2.21s/it]
</> Preprocessing gff at ../raw_data\cyanobacteria.gff3
</> Preprocessing fasta at ../raw_data\cyanobacteria.fasta
    Extracting codons: 100%|███████████████████████████████████████████████████████████| 75/75 [03:06<00:00,  2.49s/it]
</> Preprocessing gff at ../raw_data\firmicutes.gff3
</> Preprocessing fasta at ../raw_data\firmicutes.fasta
    Extracting codons: 100%|███████████████████████████████████████████████████████████| 75/75 [02:07<00:00,  1.70s/it]
</> Preprocessing gff at ../raw_data\fusobacteria.gff3
</> Preprocessing fasta at ../raw_data\fusobacteria.fasta
    Extracting codons: 100%|███████████████████████████████████████████████████████████| 75/75 [01:32<00:00,  1.23s/it]
</> Preprocessing gff at ../raw_data\proteobacteria.gff3
</> Preprocessing fasta at ../raw_data\proteobacteria.fasta
    Extracting codons: 100%|███████████████████████████████████████████████████████████| 75/75 [02:15<00:00,  1.81s/it]
</> Preprocessing gff at ../raw_data\spirochetes.gff3
</> Preprocessing fasta at ../raw_data\spirochetes.fasta
    Extracting codons: 100%|███████████████████████████████████████████████████████████| 75/75 [01:11<00:00,  1.05it/s]
```

The main script is written in `main.py` involving statistical visualization, cluster analysis and classification. Run this code in a terminal with the following command:
```
$ python main.py
================
</> Linear Classification
</> Scores of prediction based on start codon:
    Phylum Actinobacteria: 0.9354838709677419
    Phylum CFB: 1.0
    Phylum Proteobacteria: 0.0
    Phylum Firmicutes: 0.7619047619047619
    Global score: 0.6565656565656566
</> DecisionTree Classification
</> Scores of prediction based on start codon:
    Phylum Actinobacteria: 0.8709677419354839
    Phylum CFB: 0.8
    Phylum Proteobacteria: 0.5185185185185185
    Phylum Firmicutes: 0.7142857142857143
    Global score: 0.7272727272727273
</> KNeighbors Classification
</> Scores of prediction based on start codon:
    Phylum Actinobacteria: 0.9354838709677419
    Phylum CFB: 0.9
    Phylum Proteobacteria: 0.8148148148148148
    Phylum Firmicutes: 0.8095238095238095
    Global score: 0.8686868686868687
</> RandomForest Classification
</> Scores of prediction based on start codon:
    Phylum Actinobacteria: 0.9354838709677419
    Phylum CFB: 0.85
    Phylum Proteobacteria: 0.8148148148148148
    Phylum Firmicutes: 0.8571428571428571
    Global score: 0.8686868686868687
```
