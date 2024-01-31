from torch.utils import data
import numpy as np
import pickle
import random


class fre_coding(data.Dataset):
    def __init__(self,fre_mat,type='train'):
        self.type = type
        self.fre_mat = fre_mat
        self.read_id,self.snv_id = np.where(self.fre_mat!=0)
        self.fre_len = len(self.read_id)
        _,self.snv_num = fre_mat.shape
    
    def __len__(self):
        return len(self.read_id)
    
    def __getitem__(self, index):
        read_snv = np.array([self.read_id[index],self.snv_id[index],], dtype=int)
        fre = self.fre_mat[self.read_id[index],self.snv_id[index]]
        return read_snv,fre

