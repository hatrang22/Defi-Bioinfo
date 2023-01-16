# Statistical exploration of codon usage

This project is part of "Défi-Bioinformatique" at INSA/ENSAT Toulouse, France, consisting of statistical exploration of start and stop codons in prokaryotes. The objective is to know the distribution of these codons according to the bacterial phyla and to allow the development of a prediction tool that, from an unknown genome, would assign the corresponding phylum based on the proportions of its start and stop codons.

The main script is written in `main.py`. Run this code in a terminal with the following command:
```
$ python main.py
</> Preprocessing gff at demo_data\fusobacteria.gff3
</> Preprocessing fasta at demo_data\fusobacteria.fasta
    Extracting codon: 100%|████████████████████████████████████████████████████████████| 75/75 [01:28<00:00,  1.18s/it]
</> Preprocessing gff at demo_data\spirochetes.gff3
</> Preprocessing fasta at demo_data\spirochetes.fasta
    Extracting codon: 100%|████████████████████████████████████████████████████████████| 75/75 [01:09<00:00,  1.08it/s]
```

