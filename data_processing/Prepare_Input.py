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
    
    return input_file
