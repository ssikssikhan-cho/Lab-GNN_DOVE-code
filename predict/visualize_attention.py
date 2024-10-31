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
from predict.predict_single_input import init_model

def Get_Attention(dataloader,device,model):
    Final_atten1 = []
    Final_atten2 = []
    with torch.no_grad():
        for batch_idx, sample in enumerate(dataloader):
            H, A1, A2, V, Atom_count = sample
            batch_size = H.size(0)
            H, A1, A2, V = H.to(device), A1.to(device), A2.to(device), V.to(device)
            atten1,atten2= model.eval_model_attention((H, A1, A2, V, Atom_count), device)
            atten1 = atten1.detach().cpu().numpy()
            atten2 = atten2.detach().cpu().numpy()
            Final_atten1 += list(atten1)
            Final_atten2 += list(atten2)
    return Final_atten1,Final_atten2

def visualize_attention(input_path,params):
    #create saving path
    save_path=os.path.join(os.getcwd(),"Predict_Result")
    mkdir(save_path)
    save_path = os.path.join(save_path, "Visulize_Target")
    mkdir(save_path)
    save_path = os.path.join(save_path, "Fold_"+str(params['fold'])+"_Result")
    mkdir(save_path)

    input_path=os.path.abspath(input_path)
    split_name=os.path.split(input_path)[1]
    original_pdb_name=split_name
    if ".pdb" in split_name:
        split_name=split_name[:-4]
    save_path=os.path.join(save_path,split_name)
    mkdir(save_path)

    #load model
    fold_choice = params['fold']
    model_path = os.path.join(os.getcwd(), "best_model")
    model_path = os.path.join(model_path, "fold" + str(fold_choice))
    model_path = os.path.join(model_path, "checkpoint.pth.tar")
    model, device = init_model(model_path, params)

    structure_path = os.path.join(save_path, "Input.pdb")
    shutil.copy(input_path, structure_path)
    input_file = Prepare_Input(structure_path)
    list_npz = [input_file]
    dataset = Single_Dataset(list_npz)
    dataloader = DataLoader(dataset, 1, shuffle=False,
                            num_workers=params['num_workers'],
                            drop_last=False, collate_fn=collate_fn)

    Final_atten1,Final_atten2 = Get_Attention(dataloader, device, model)
    tmp_save_path1 = os.path.join(save_path, "attention1.npy")
    tmp_save_path2 = os.path.join(save_path, "attention2.npy")
    np.save(tmp_save_path1, Final_atten1)
    np.save(tmp_save_path2, Final_atten2)

    atom_count = input_file[1]
    print(f'Atom Count : {atom_count}')