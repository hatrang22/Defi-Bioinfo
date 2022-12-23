import  gzip
import os
import tarfile
import warnings

def read_phylum(fasta_path, gff_path):
    """
    Parameters
    ----------
    fasta_path: str
            Path to fasta tar file.

    gff_path: str
            Path to gff tar file.

    Returns
    -------
    ret: dict
        Dictionary with 3 keys:

            - accession_number: list of accession numbers
            - fasta: list of corresponding fastas
            - gff: list of corresponding gffs

    """

    fasta_dir = os.path.splitext(fasta_path)[0] # remove extension file
    with tarfile.open(fasta_path) as f:
        f.extractall(fasta_dir)
    fasta_dir = os.path.join(fasta_dir, os.listdir(fasta_dir)[0]) # files extracted in fasta_dir

    gff_dir = os.path.splitext(gff_path)[0] # remove extension file
    with tarfile.open(gff_path) as f:
        f.extractall(gff_dir)
    gff_dir = os.path.join(gff_dir, os.listdir(gff_dir)[0]) # files extracted in gffr_dir

    list_fasta = [f for f in os.listdir(fasta_dir) if f.endswith(".gz")]
    list_fasta.sort()  # list of sorted fasta names

    list_gff = [f for f in os.listdir(gff_dir) if f.endswith(".gz")]
    list_gff.sort()  # list of sorted gff names

    ret = {"accession_number": [], "fasta": [], "gff": []}  # returned dict

    for fasta_fn, gff_fn in zip(list_fasta, list_gff):

        fan = fasta_fn.split("_")
        fan = "_".join(fan[:2])  # get accession number for fasta

        gan = gff_fn.split("_")
        gan = "_".join(gan[:2])  # get accession number for gff

        if fan == gan:  # check if accession numbers are consistent

            fasta, gff = preprocessing_fasta_gff(os.path.join(fasta_dir, fasta_fn), os.path.join(gff_dir, gff_fn))

            ret["fasta"].append(fasta)

            ret["gff"].append(gff)

            ret["accession_number"].append(fan)

        else:
            warnings.warn(f"Inconsistent accession number between fasta file ({fan}) and gff file ({gan})")
            
    return ret

def preprocessing_fasta_gff(fasta_gz, gff_gz):

    with gzip.open(fasta_gz, "rt") as gzip_fasta:

        fasta= []
        for f in gzip_fasta.readlines():
            if f[0].isalpha():
                fasta.append(f)
            elif fasta:
                break
        fasta = "".join(fasta) # concatener list de string en str
        fasta = fasta.replace("\n","") # replace \n

    with gzip.open(gff_gz, "rt") as gzip_gff:

        gff = [f.split("\t") for f in gzip_gff.readlines() if f[0].isalpha()] # Conserver uniquement les lignes commencant par une lettre

    return fasta, gff