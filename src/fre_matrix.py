import sys
import pysam
import numpy as np
import pickle
import time
import concurrent.futures
from tqdm import tqdm
import multiprocessing as mp
mp.set_start_method("fork")
### load data
file_in=sys.argv[1]  
file_snv = sys.argv[2] 
file_prefix = sys.argv[3] 
cond_pro = float(sys.argv[4])
smallest_snv = int(sys.argv[5])
num_read_1 = int(sys.argv[6])
num_read_2 = int(sys.argv[7])
gap = int(sys.argv[8])
fre_most_base = float(sys.argv[9])
thread = int(sys.argv[10])
input_snv = sys.argv[11]


bamfile = pysam.AlignmentFile(file_in)
ref_name = bamfile.references[0]
f_snv = open(file_snv,'r')
snv_sites = f_snv.readline().split()
snv_sites = np.array([int(i)-1 for i in snv_sites])
f_snv.close()
ref_length = bamfile.get_reference_length(ref_name)

### data initialization
def error_rate(qual):
    er = 10 ** (-qual / 10)
    return er

Reads_all = []
seq_mat = []
for read in bamfile.fetch():
    align_index = np.array(read.get_aligned_pairs(matches_only=True))
    align_seq = np.array(list(read.query_sequence))
    read_seq = np.array(['0']*ref_length)
#    read_seq[align_index[0,1]:align_index[-1,1]+1]='-'
    read_seq[list(align_index[:,1])]= align_seq[list(align_index[:,0])]
    if not (read_seq[snv_sites] == '0').all():
        Reads_all.append(read.query_name)
        seq_mat.append(read_seq[snv_sites])
    

seq_mat = np.array(seq_mat)

num_snv = len(snv_sites)
num_reads = len(Reads_all)

Coun = []
for i in range(len(snv_sites)):
    chara = ['A','C','G','T']
    temp = []
    for j in chara:
        temp.append(np.count_nonzero(seq_mat[:,i]==j))
    Coun.append(temp)


Coun = np.array(Coun)

####  calculate nucleotide frequency and rank
def rank_nucl(Coun):
	res = []
	for i in Coun:
		res.append(list(np.sort(i)/sum(i)))
	return res


nucl_fre = np.array(rank_nucl(Coun))

if input_snv =='None':
    print("strat maximum conditional probability")

###   find out the candicate site to be verified >=fre_most_base  
to_be_verified_sites = np.array([i for i in range(num_snv)])
to_be_verified_sites = to_be_verified_sites[nucl_fre[:,-1]>=fre_most_base]

###  calculate the 1-order  probability for the SNV at each candicate site
subsitution_rate = {i:nucl_fre[i,2]/sum(nucl_fre[i,2:]) for i in to_be_verified_sites}

###  reserve the top-2 bases  
chara = ['A','C','G','T']
snv_base = []
major_base = []
for i in range(len(Coun)):
    temp = np.argsort(Coun[i])
    snv_base.append(chara[temp[-2]])
    major_base.append(chara[temp[-1]])

snv_base = np.array(snv_base)
major_base = np.array(major_base)

seq_snv_flag = []
seq_major_flag = []
for i in range(len(seq_mat)):
    temp = seq_mat[i]==snv_base
    seq_snv_flag.append(temp)
    temp = seq_mat[i]==major_base
    seq_major_flag.append(temp)

seq_snv_flag = np.array(seq_snv_flag)
seq_major_flag = np.array(seq_major_flag)

###   calculate the conditional probability
def cond_pro_1_order(c_site,g_site,ssf=seq_snv_flag,smf=seq_major_flag,v_min = num_read_1):
    #c_site: current site
    #g_site: given site
    #ssf: seq_snv_flag
    #smf: seq_major_flag
    current_snv_flag = ssf[:,c_site]
    current_major_flag = smf[:,c_site]
    given_flag = ssf[:,g_site]
    num_snv = np.count_nonzero(current_snv_flag[given_flag])
    num_major = np.count_nonzero(current_major_flag[given_flag])
    denominator  = num_snv + num_major
    numerator = num_snv
    if denominator > v_min:
        pro = numerator/(denominator+1e-16)
    else:
        pro = 0
    return pro

def cond_pro_extend(c_site,g_site,cg_flag,ssf,smf):
    #c_site: current site
    #g_site: given site
    #ssf: seq_snv_flag
    #smf: seq_major_flag
    #current_given_flag
    current_snv_flag = ssf[:,c_site]
    current_major_flag = smf[:,c_site]
    given_flag = np.logical_and(ssf[:,g_site],cg_flag)
    num_snv = np.count_nonzero(current_snv_flag[given_flag])
    num_major = np.count_nonzero(current_major_flag[given_flag])
    denominator = num_snv + num_major
    numerator = num_snv
    pro = numerator/(denominator+1e-16)
    return pro,given_flag,denominator

def verify_snv_site(site_id,to_be_verified_sites=to_be_verified_sites,seq_snv_flag=seq_snv_flag,seq_major_flag=seq_major_flag,v_min_1=num_read_1,v_min_2=num_read_2,hd=gap):
    temp_flag = to_be_verified_sites!=site_id
    given_sites = to_be_verified_sites[temp_flag]
    con_pro = []
    for j in given_sites:
        con_pro.append(cond_pro_1_order(site_id,j))
    con_pro = np.array(con_pro)
    index_sort = np.argsort(-con_pro)
    current_pro = subsitution_rate[site_id]
    current_given_flag = np.array([True]*num_reads)
    for k in index_sort:
        given_site = given_sites[k]
        if abs(snv_sites[given_site]-snv_sites[site_id]) < hd:
            continue
        new_pro,new_given_flag,new_v = cond_pro_extend(site_id,given_site,current_given_flag,seq_snv_flag,seq_major_flag)
        if new_pro>current_pro and new_v>=v_min_2:
            current_pro = new_pro
            current_given_flag = new_given_flag
            if current_pro>=cond_pro:
                break
        else:
            break
    return current_pro

if input_snv=='None':
    query_list = [i for i in to_be_verified_sites]
    with concurrent.futures.ProcessPoolExecutor(thread) as executor:
        snv_max_pro = list(tqdm(executor.map(verify_snv_site, query_list), total=len(query_list)))
    snv_max_pro = np.array(snv_max_pro)
    fake_snv_sites = to_be_verified_sites[snv_max_pro<cond_pro]
    del snv_max_pro
    final_snv = [i for i in range(num_snv) if i not in fake_snv_sites]
    ###  update snv_sites file
    f = open(file_snv,'w')
    for i in final_snv:
        f.write(str(snv_sites[i]+1)+'\t')
    ### check the satisfaction of SNV sites's number
    if len(final_snv)<smallest_snv:
        print("No enough SNV sites.\nexit")
        f.write("\nA small number of SNV sites indicates one haplotype\nexit")
        f.close()
        exit()
    f.close()
else:
    final_snv = [i for i in range(num_snv)]


###  update some variables
Coun = Coun[final_snv]
seq_mat = seq_mat[:,final_snv]
snv_base = snv_base[final_snv]
major_base = major_base[final_snv]



data = {'seq_mat':seq_mat,'Reads_all':Reads_all,'Coun':Coun,'snv_base':snv_base,'major_base':major_base}
f_mat = open(file_prefix+'_matrix.pickle','wb')
pickle.dump(data,f_mat,protocol=pickle.HIGHEST_PROTOCOL)
f_mat.close()
del data
exit()
