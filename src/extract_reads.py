import os
import sys
import numpy as np
import pickle
from Bio import SeqIO

file_path = sys.argv[1]
file_prefix = sys.argv[2]
flag = sys.argv[3]

os.system(f'rm -rf {file_path}/medaka/fastx')
os.system(f'mkdir -p {file_path}/medaka/fastx')

if flag == '1':
    os.system(f'samtools fastq {file_path}/alignment/{file_prefix}.bam > {file_path}/medaka/fastx/cluster.fastq')
    os.system(f'cp {file_path}/{file_prefix}_consensus.fasta {file_path}/medaka/fastx/consensus.fasta')
else:
    file = f'{file_path}/{file_prefix}_consensus.fasta'
    num = 0
    for seq in SeqIO.parse(file,'fasta'):
        num = num + 1
    f = open(f'{file_path}/{file_prefix}_consensus.fasta','r')
    for i in range(num):
        os.system(f'samtools fastq {file_path}/clusters/cluster_{i}.bam > {file_path}/medaka/fastx/cluster_{i}.fastq')
        f_o = open(f'{file_path}/medaka/fastx/consensus_{i}.fasta','w')
        f_o.write(f.readline())
        f_o.write(f.readline())
        f_o.close()

exit()
        
    
