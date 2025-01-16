import sys
import os

# assign fasta file and msa folder from command line
# sys.argv[0] is the script name itself
fasta_file = sys.argv[1]
msa_folder = sys.argv[2]

# read fasta file
with open(fasta_file, "r") as f:
    fasta = f.read().lstrip(">").split(">")

# convert fasta file to dictionary (store index, header and protein sequence)
fasta_dict = {}
for i in fasta:
    header, protein_seq = i.rstrip("\n").split("\n")
    index = str(fasta.index(i))
    fasta_dict[index] = [header, protein_seq]

# make a list of a3m files ans sort them by file name (index)
a3m_files = [f for f in os.listdir(msa_folder) if f.endswith('.a3m')]
a3m_files.sort()

# go through fasta dictionary and a3m files and compare index and protein sequence
for k, v in fasta_dict.items():

    for file in a3m_files:

        # check if a3m file name and fasta index matches
        file_index = file.split(".")[0]
        if k == file_index:

            with open(os.path.join(msa_folder, file), "r") as f:
                a3m_content = f.read().lstrip(">").split(">")
            
            # get protein sequence for the probe from a3m file
            for i in a3m_content:
                if i.startswith("101\n"):
                    header, protein_seq = i.rstrip("\n").split("\n")
                    break
            
            # rename a3m files with headers from fasta file if index and protein sequence match
            if v[1] == protein_seq:
                os.rename(os.path.join(msa_folder, file), os.path.join(msa_folder, v[0] + ".a3m"))
            else:
                os.rename(os.path.join(msa_folder, file), os.path.join(msa_folder, v[0] + ".a3m"))
                print(f"SEQS_DONT_MATCH:\t{msa_folder}/{v[0]}.a3m\t{protein_seq}\t{v[1]}")