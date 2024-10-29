import random
import torch
import torch.nn.functional as F
import torch.nn as nn
import time
from multiprocessing import Pool
from model.layers import GAT_gate
import os
from ops.os_operation import mkdir
import shutil
import  numpy as np
from data_processing.Prepare_Input import Prepare_Input
from model.GNN_Model import GNN_Model
import torch
from ops.train_utils import count_parameters,initialize_model
from data_processing.collate_fn import collate_fn
from data_processing.Single_Dataset import Single_Dataset
from torch.utils.data import DataLoader
from ops.argparser import argparser
import os
import time
from prepare_learning_data import AverageMeter
from ops.argparser import argparser
from model.GNN_Model import GNN_Model
import torch
import os
from prepare_learning_data import collate_fn
from ops.train_utils import count_parameters,initialize_model
from data_processing.Single_Dataset import Single_Dataset
from torch.utils.data import DataLoader
import torch.nn as nn
from model.GNN_Model import GNN_Model
import torch
from ops.train_utils import count_parameters,initialize_model
from data_processing.Single_Dataset import Single_Dataset
from torch.utils.data import DataLoader
from ops.argparser import argparser
import time

def train_GNN(model,train_dataloader,optimizer,loss_fn,device):

    model.train()

    Loss = AverageMeter()
    Accu1 = AverageMeter()
    end_time = time.time()
    
    for batch_idx, sample in enumerate(train_dataloader):
        # start = time.time()
        b = time.time()
        H, A1, A2, V, Atom_count, Y = sample
        batch_size = H.size(0)
        H, A1, A2, Y, V = H.cuda(), A1.cuda(), A2.cuda(), Y.cuda(), V.cuda()
        # end = time.time()
        # print("loop end, train start: ", end - start)
        pred = model.train_model((H, A1, A2, V, Atom_count), device)
        # start = time.time()
        # print("end train: ", start - end)

        loss = loss_fn(pred, Y)
        # end = time.time()
        # print("loss end: ", end - start)
        optimizer.zero_grad()
        # start = time.time()
        # print("optimizer end:", start - end)
        #torch.nn.utils.clip_grad_norm_(model.parameters(), params['clip'])
        torch.nn.utils.clip_grad_value_(model.parameters(), clip_value = 1.0)
        #loss_list.append(loss)

        loss.backward()

        optimizer.step()

        Accu1.update(pred, batch_size)
        Loss.update(loss.item(), batch_size)
        # end = time.time()

        # print("batch end: ",batch_idx, time.time() - b)


    return Loss.avg, Accu1.avg#, loss_list