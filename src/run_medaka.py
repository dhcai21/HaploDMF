import os
import sys
import numpy as np
import pickle
from Bio import SeqIO

file_path = sys.argv[1]
file_prefix = sys.argv[2]
flag = sys.argv[3]

if flag == '1':
    os.system(f'medaka_consensus -i {file_path}/medaka/fastx/cluster.fastq -d {file_path}/medaka/fastx/consensus.fasta -o {file_path}/medaka/medaka0')
    os.system(f'cp {file_path}/medaka/medaka0/consensus.fasta {file_path}/{file_prefix}_haplotypes.fasta')
else:
    file = f'{file_path}/{file_prefix}_consensus.fasta'
    num = 0
    for seq in SeqIO.parse(file,'fasta'):
        num = num + 1
    for i in range(num):
        os.system(f"rm -rf {file_path}/medaka/medaka_{i}")
        os.system(f'medaka_consensus -i {file_path}/medaka/fastx/cluster_{i}.fastq -d {file_path}/medaka/fastx/consensus_{i}.fasta -o {file_path}/medaka/medaka_{i}')
    os.system(f"cat {file_path}/medaka/medaka_{0}/consensus.fasta > {file_path}/{file_prefix}_haplotypes.fasta")
    for i in range(1,num):
        os.system(f"cat {file_path}/medaka/medaka_{i}/consensus.fasta >> {file_path}/{file_prefix}_haplotypes.fasta")

exit()
