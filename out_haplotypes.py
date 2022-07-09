import sys
import pickle
import numpy as np
import pysam
import pandas as pd
import os
## load clusters and initialize
##
file_prefix = sys.argv[1]
file_bam= sys.argv[2]
file_path=sys.argv[3]
file_ref= sys.argv[4]
file_clusters = f"{file_path}/{file_prefix}_clusters.pickle"
file_out = f"{file_path}/{file_prefix}_consensus.fasta" #input
os.system(f"samtools index {file_bam}")
os.system(f"mkdir -p {file_path}/clusters")
s_pos = int(sys.argv[5])
e_pos = int(sys.argv[6])


f = open(file_clusters,'rb')
data = pickle.load(f)
f.close()
reads = data['Reads_all']
label = np.array(data['label'])
reads = np.array(reads)
haplo_cluster = []
haplo_fre = []
i = 0
while True:
    flag = label == i
    count = np.count_nonzero(flag)
    if count == 0:
        break
    haplo_cluster.append(flag)
    haplo_fre.append(count/len(label))
    i = i+1

haplo_fre = np.array(haplo_fre)
index = list(np.argsort(-haplo_fre))
haplo_fre = haplo_fre[index]
haplo_cluster = [haplo_cluster[i] for i in index]
del data
def out_reads(reads,haplo):
    R = []
    for h in haplo:
        r_temp = reads[h]
        R.append(r_temp)
    return R

R = out_reads(reads,haplo_cluster)


attributes = ['pos', 'reads_all','deletions','A', 'C', 'G', 'T']
chara = 'ACGT'

f_out = open(file_out,'w')

for i in range(len(R)):
    bam = pysam.AlignmentFile(file_bam)
    qnames = R[i]
    out_name = file_path+"/clusters/cluster_"+str(i)+".bam"
    obam = pysam.AlignmentFile(out_name, "wb", template= bam)
    for b in bam.fetch(until_eof=True):
        if b.query_name in qnames:
            obam.write(b)
    obam.close()
    bam.close()
    file_cluster_bam = out_name
    file_cluster_bam_sorted = file_cluster_bam[0:-4]+'_sorted.bam'
    os.system(f"samtools sort {file_cluster_bam} -o {file_cluster_bam_sorted}")
    os.system(f"samtools index {file_cluster_bam_sorted}")
    os.system(f"python ./src/count_frequency.py {file_cluster_bam_sorted} {file_path}/clusters/cluster_{i}_acgt.txt")
    seq_temp = []
    file_acgt = file_path+"/clusters/cluster_"+str(i) + '_acgt.txt'
    df = pd.read_csv(file_acgt, sep='\t')
    data = df[attributes]
    mat = data.values
    flag = np.logical_and(mat[:,0]>=s_pos,mat[:,0]<=e_pos)
    mat = mat[flag,:]
    for j in range(len(mat)):
        pos = mat[j, 0] - 1
        if mat[j,2]/mat[j,1]>=0.5:
            continue
        index = np.argsort(-mat[j,3:])[0]
        base = chara[index]
        seq_temp.append(base)
    h = ''.join(seq_temp)
    f_out.write('>haplotype_'+str(i)+'_length_'+str(len(h))+'_abundance_'+str(haplo_fre[i])+'\n')
    f_out.write(h+'\n')

f_out.close()
exit()
