import sys
import pickle
import numpy as np
import torch.nn.functional as F
import random
import torch
from torch.utils.tensorboard import SummaryWriter
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from collections import Counter
from torch.utils.data import DataLoader
import torch.nn as nn
from tqdm import tqdm
import pandas as pd
from typing import Optional
import os
import time

sys.path.append('./src/')
import dmf_model
import dl_dataset
import utils

# Device configuration
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(device)

# Seed
seed = 21
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)

## weighted loss function
def weighted_mse_loss(input, target):
    weight = 1/(target+1e-16)
    return torch.sum(torch.multiply(weight,(input - target) ** 2))

def snv_loss(input,target):
    up = torch.matmul(input,target.T)
    temp1 = torch.unsqueeze(torch.linalg.vector_norm(input,dim=1),1)
    temp2 = torch.unsqueeze(torch.linalg.vector_norm(target,dim=1),0)
    down = torch.matmul(temp1,temp2)
    res = torch.sum((-torch.divide(up,down)+1)*0.5)
    return res

def fre_cal(Coun,major_base,snv_base):
    chara = {'A':0,'C':1,'G':2,'T':3}
    res_mat = []
    for i in range(len(major_base)):
        top1 = Coun[i,chara[major_base[i]]]
        top2 = Coun[i,chara[snv_base[i]]]
        top12 = top1+top2
        res_mat.append([top1/top12,top2/top12])
    res_mat = np.array(res_mat)
    return res_mat.T

def fre_mat_create(seq_mat,fre_flag,major_base,snv_base):
    m,n = seq_mat.shape
    res_fre_mat = np.zeros([m,n])
    major_fre = fre_flag[0]
    snv_fre = fre_flag[1]
    for i in range(m):
        seq = seq_mat[i]
        flag_1 = seq == major_base
        flag_2 = seq == snv_base
        res_fre_mat[i] = flag_1 * major_fre + flag_2*snv_fre
    return res_fre_mat


if __name__ == "__main__":
    args = utils.cmd()
    # load data
    print(args)
    train_path = args.i
    f = open(train_path,'rb')
    data = pickle.load(f)
    seq_mat = data['seq_mat']
    major_base = data['major_base']
    snv_base = data['snv_base']
    Coun = data['Coun']
    fre_flag = fre_cal(Coun,major_base,snv_base)
    fre_mat = fre_mat_create(seq_mat,fre_flag,major_base,snv_base)
    Reads_all = data['Reads_all']
    train_data = eval(f"dl_dataset.fre_coding")(fre_mat)
    #network parameter
    batch_size = int(args.batch_size)
    # DON'T shuffle the testing datasets
    train_loader = DataLoader(train_data, batch_size=batch_size, shuffle=True, num_workers=1)
    del train_data

# train function
def train(args, train_loader, fre_mat,read_mid_dim = 128, snv_mid_dim=256, final_vector = 32):
    # create a log directory to save the logging
    local_time = time.strftime("%m-%d_%H-%M", time.localtime())
    # determine log file name and `mkdir`
    log_file_name = args.p
    os.system(f"mkdir -p {args.o}/log/{log_file_name}")
    writer = SummaryWriter(f"{args.o}/log/{log_file_name}/tensorboard")
    # load model
    model = eval(f"dmf_model.dfm")(
        fre_mat=fre_mat, read_mid_dim=read_mid_dim, snv_mid_dim= snv_mid_dim, final_vector=final_vector
        ).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.l)
    total_step = len(train_loader)
    n_epochs = args.n
    lw = args.w
    # the time of starting training
    train_s_time = time.time()
    _,snv_num = fre_mat.shape
    ref_index = torch.LongTensor([i for i in range(snv_num)])
    ref_index = ref_index.to(device)
    for epoch in range(n_epochs):  # times of feeding data to the network
        train_losses = []
        # Train the model
        model.train()
        for (datas, fre) in train_loader:#tqdm(train_loader):
            datas = datas.to(device)
            fre = fre.float()
            fre = fre.to(device)
            # Forward pass
            # loss for predict fre
            _,snv_vec,outputs = model(datas[:,0],datas[:,1])
            loss1 = weighted_mse_loss(outputs,fre)
            # ref_snv output
            ref_mat = model.snv_embed(ref_index)
            ref_snv = F.relu(model.snv_layer_0(ref_mat))
            ref_snv = model.snv_layer_1(ref_snv)
            # loss of distances between snv sites
            loss2 = snv_loss(snv_vec,ref_snv)
            # Backward and optimize
            loss = (1-lw)*loss1 + lw/(snv_num-1)*loss2
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            train_losses.append(loss.item())
        # calculate average loss over an epoch
        train_loss = np.average(train_losses)
        writer.add_scalars('loss', {'train': train_loss}, epoch)
        # print training/validation statistics
        epoch_len = len(str(n_epochs))
        print(f'Epoch [{epoch + 1:>{epoch_len}}/{n_epochs:>{epoch_len}}] ' +
              f'train_loss: {train_loss:.5f} ')
    # the time of ending training
    train_e_time = time.time()
    print(f"Training time: {(train_e_time - train_s_time) / 60 :.2f} min.")
    # save the model to the log file
    torch.save(model.state_dict(), f"{args.o}/log/{log_file_name}/{log_file_name}.pt")
    # the time of ending training
    train_e_time = time.time()
    print(f"Training time: {(train_e_time - train_s_time) / 60 :.2f} min.")
    # save the model to the log file
    torch.save(model.state_dict(), f"{args.o}/log/{log_file_name}/{log_file_name}.pt")
    return model,log_file_name

def trans(label,num_haplo):
    cluster_index = {i:[] for i in range(num_haplo)}
    for i in range(len(label)):
        cluster_index[label[i]].append(i)
    for key in cluster_index.keys():
        cluster_index[key] = np.array(cluster_index[key])
    return cluster_index

### training
_,log_file_name = train(args = args, train_loader=train_loader, fre_mat = fre_mat, read_mid_dim = 64, snv_mid_dim=128, final_vector = 32)
model = eval(f"dmf_model.dfm")(fre_mat=fre_mat, read_mid_dim=64, snv_mid_dim= 128, final_vector=32)
del _
del train_loader
model.load_state_dict(torch.load(f"{args.o}/log/{log_file_name}/{log_file_name}.pt",map_location=torch.device('cpu')))

####  output learned vector
read_num,snv_num = fre_mat.shape
read_snv_once = []
for i in range(read_num):
    read_snv_once.append([i,0])

read_snv_once = torch.tensor(read_snv_once)
read_snv_once = read_snv_once#.to(device)
read_vec,snv_vec,_ = model(read_snv_once[:,0],read_snv_once[:,1])

feature = read_vec.cpu().detach().numpy()

del model
del fre_mat

##### normalize vector
normalize_feature = []
for i in feature:
    normalize_feature.append(i/np.sqrt(np.sum(i**2)))

normalize_feature = np.array(normalize_feature)

#### determine the number of haplotypes
seq_flag = (seq_mat == major_base)+(seq_mat == snv_base)
def call_consensus(seq_mat,cluster_index,depth=5):
    consensus = {}
    consensus_flag = {}
    cluster_label = {}
    for key in cluster_index.keys():
        sub_seq_mat = seq_mat[cluster_index[key]]
        seq = np.array(['0']*snv_num)
        for i in range(snv_num):
            temp = Counter(sub_seq_mat[:, i])
            if '-' in temp.keys():
                temp['-']=0
            elif '0' in temp.keys():
                temp['0']=0
            max_value = max(temp.values())
            values_list = np.array(list(temp.values()))
            if np.sum(values_list)>=depth:
                nucl = max(temp, key=temp.get)
                seq[i] = nucl
        temp = seq != '0'
        if np.count_nonzero(seq_flag)==0:
            cluster_label[key] = False
        else:
            cluster_label[key] = True
        consensus[key] = seq
        consensus_flag[key] = temp
    return consensus,consensus_flag,cluster_label


def max_error_correction(con_seq,con_flag,clus_label,seq_mat,seq_flag,cluster_index):
    total_mec={}
    for key in con_seq.keys():
        sub_seq_mat = seq_mat[cluster_index[key]]
        sub_seq_flag = seq_flag[cluster_index[key]]
        num_temp,_ = sub_seq_mat.shape
        consensus = con_seq[key]
        consensus_flag = con_flag[key]
        if clus_label[key]:
            temp_mec = sub_seq_mat!=consensus
            info = []
            for i in range(num_temp):
                compare_flag = np.logical_and(sub_seq_flag[i],consensus_flag)
                edit_flag = temp_mec[i][compare_flag]
                edit = np.count_nonzero(edit_flag)
                info.append(edit)
            info = np.array(info)
            total_mec[key] = np.sum(info)
    return total_mec



def find_clusters(args,normalize_feature,seq_mat,seq_flag):
    upper_clusters = args.largest_num
    algorithm = args.algorithm
    a,b = seq_mat.shape
    mec_all = []
    final_num = 2
    step = 1
    label_list = [[0]*a]
    cluster_index = {}
    cluster_index[0] = np.array([i for i in range(len(seq_mat))])
    con_seq,con_flag,clus_label = call_consensus(seq_mat,cluster_index)
    total_mec = max_error_correction(con_seq,con_flag,clus_label,seq_mat,seq_flag,cluster_index)
    mec_all.append(sum(total_mec.values()))
    for num_haplo in range(2,upper_clusters):
        print(num_haplo)
        if algorithm == 'ward':
            clustering = AgglomerativeClustering(n_clusters=num_haplo).fit(normalize_feature)
        elif algorithm == 'kmeans':
            clustering = KMeans(n_clusters=num_haplo,random_state=0).fit(normalize_feature)
        else:
            clustering = AgglomerativeClustering(n_clusters=num_haplo).fit(normalize_feature)
        label = clustering.labels_
        label_list.append(label)
        cluster_index = trans(label,num_haplo)
        con_seq,con_flag,clus_label = call_consensus(seq_mat,cluster_index)
        total_mec = max_error_correction(con_seq,con_flag,clus_label,seq_mat,seq_flag,cluster_index)
        mec_all.append(sum(total_mec.values()))
        thres = mec_all[-1]/mec_all[-2]
        final_label = label_list[-2]
        final_num = num_haplo-1
        if thres >args.c1:
            if algorithm == 'ward':
                clustering = AgglomerativeClustering(n_clusters=num_haplo+step).fit(normalize_feature)
            elif algorithm == 'KMeans':
                clustering = KMeans(n_clusters=num_haplo+step,random_state=0).fit(normalize_feature)
            else:
                clustering = AgglomerativeClustering(n_clusters=num_haplo+step).fit(normalize_feature)
            label = clustering.labels_
            cluster_index = trans(label,num_haplo+step)
            con_seq,con_flag,clus_label = call_consensus(seq_mat,cluster_index)
            total_mec = max_error_correction(con_seq,con_flag,clus_label,seq_mat,seq_flag,cluster_index)
            ref_mec = sum(total_mec.values())
            final_thres = ref_mec/mec_all[-1]
            if final_thres > args.c2:
                if final_num == 1:
                    print("only one cluster")
                break
    return mec_all,final_num,final_label




#### clusters via vector
print(log_file_name)

### hierarhical 
mec_list,final_num,final_label = find_clusters(args,normalize_feature,seq_mat,seq_flag)
print(mec_list)

### save clusters for calling consensus
data = {'label':final_label,'Reads_all':Reads_all,'log_file_name':log_file_name}
f= open(f"{args.o}/{args.p}_clusters.pickle",'wb')

pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
f.close()
exit()
