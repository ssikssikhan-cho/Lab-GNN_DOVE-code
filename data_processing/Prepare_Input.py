import os
from data_processing.Extract_Interface import Extract_Interface
from rdkit.Chem.rdmolfiles import MolFromPDBFile
from data_processing.Feature_Processing import get_atom_feature
import numpy as np
from rdkit.Chem.rdmolops import GetAdjacencyMatrix
from scipy.spatial import distance_matrix


def Prepare_Input(structure_path):
    # extract the interface region
    interface_path = Extract_Interface(structure_path)
<<<<<<< HEAD
    protein_mol = MolFromPDBFile(interface_path, sanitize=False)
    

    if protein_mol is None:
        raise ValueError("fail to create from PDB-file")
    atom_count = protein_mol.GetNumAtoms()
    
    # extract atom feature 
    atom_features = get_atom_feature(protein_mol)


    # get protain adj matrix
    conformer = protein_mol.GetConformers()[0]
    positions = np.array(conformer.GetPositions())
    adj_matrix = GetAdjacencyMatrix(protein_mol) + np.eye(atom_count)
    distance_mat = distance_matrix(positions, positions)
    # combine analysis

    valid=np.ones((atom_count,))

    root_path = os.path.dirname(interface_path)
    input_file = os.path.join(root_path, "Input.npz")
    np.savez(input_file, H=atom_features, A1=adj_matrix, A2=distance_mat, V=valid)
    
#protain_mol = MolFromPDBFile(interface_path, sanitize=False)
#atom_count = protain_mol.GetNumAtoms()
    
    # extract atom feature 
#atom_features = get_atom_feature(mol)


    # get protain adj matrix
#conformer = mol.GetConformers()[0]
#positions = np.array(conformer.GetPositions())
#adj_matrix = GetAdjacencyMatrix(mol) + np.eye(atom_count)
#distance_mat = distance_matrix(positions, positions)
    # combine analysis

#H = np.concatenate([receptor_feature, ligand_feature], 0)
#agg_adj1 = np.zeros((receptor_count + ligand_count, receptor_count + ligand_count))
#agg_adj1[:receptor_count, :receptor_count] = adj1
#agg_adj1[receptor_count:, receptor_count:] = adj2  # array without r-l interaction
#dm = distance_matrix(d1, d2)
#agg_adj2 = np.copy(agg_adj1)
#agg_adj2[:receptor_count, receptor_count:] = np.copy(dm)
#agg_adj2[receptor_count:, :receptor_count] = np.copy(np.transpose(dm))  # with interaction array
    # node indice for aggregation
#valid = np.zeros((receptor_count + ligand_count,))
#valid[:receptor_count] = 1
#input_file=os.path.join(root_path,"Input.npz")
    # sample = {
    #     'H': H.tolist(),
    #     'A1': agg_adj1.tolist(),
    #     'A2': agg_adj2.tolist(),
    #     'V': valid,
    #     'key': structure_path,
    # }
#np.savez(input_file,  H=H, A1=agg_adj1, A2=agg_adj2, V=valid)

    return input_file

    # 모든 행렬을 하나의 구조로 결합
    valid_mask = np.ones(atom_count)
    output_path = os.path.join(os.path.split(structure_path)[0], "Input.npz")
    np.savez(output_path, H=atom_features, A1=adj_matrix, A2=distance_mat, V=valid_mask)

    return output_path



    interface_path = Extract_Interface(structure_path)
    mol = MolFromPDBFile(interface_path, sanitize=False)
    atom_count = mol.GetNumAtoms()

    # 원자 특징 추출
    atom_features = get_atom_feature(mol)

    # 인접 행렬 및 거리 행렬 계산
    conformer = mol.GetConformers()[0]
    positions = np.array(conformer.GetPositions())
    adj_matrix = GetAdjacencyMatrix(mol) + np.eye(atom_count)
    distance_mat = distance_matrix(positions, positions)

    # 모든 행렬을 하나의 구조로 결합
    valid_mask = np.ones(atom_count)
    output_path = os.path.join(os.path.split(structure_path)[0], "Input.npz")
    np.savez(output_path, H=atom_features, A1=adj_matrix, A2=distance_mat, V=valid_mask)

    return output_path
