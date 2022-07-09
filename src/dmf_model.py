import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable
import numpy as np

class dfm(nn.Module):
    def __init__(self,fre_mat, read_mid_dim,snv_mid_dim, final_vector=32, dropout=0):
        super(dfm, self).__init__()
        self.read_num, self.snv_num = fre_mat.shape
        self.read_hidden = read_mid_dim
        self.snv_hidden = snv_mid_dim
        self.final_vector = final_vector
        ## embedding
        self.read_embed = nn.Embedding.from_pretrained(torch.FloatTensor(fre_mat),freeze=True)
        self.snv_embed = nn.Embedding.from_pretrained(torch.FloatTensor(fre_mat.T),freeze=True)
        ## network structure
        self.read_layer_0 = nn.Linear(self.snv_num,self.read_hidden)
        self.read_layer_1 = nn.Linear(self.read_hidden,self.final_vector)
        self.snv_layer_0 = nn.Linear(self.read_num,self.snv_hidden)
        self.snv_layer_1 = nn.Linear(self.snv_hidden,self.final_vector)
    def forward(self, read_id, snv_id):
        ## encoder layer
        read_mat = self.read_embed(read_id)
        snv_mat = self.snv_embed(snv_id)
        read_hidden = F.relu(self.read_layer_0(read_mat))
        snv_hidden = F.relu(self.snv_layer_0(snv_mat))
        read_vector = self.read_layer_1(read_hidden)
        snv_vector = self.snv_layer_1(snv_hidden)
        y = F.cosine_similarity(read_vector,snv_vector)
        y = torch.clamp(y,1e-6)
        return read_vector,snv_vector,y
