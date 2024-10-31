import os
from ops.Timer_Control import set_timeout,after_timeout
RESIDUE_Forbidden_SET={"FAD"}

def Extract_Interface(pdb_path):
    """
    specially for 2 docking models
    :param pdb_path:docking model path
    :rcount: receptor atom numbers
    :return:
    extract a receptor and ligand, meanwhile, write two files of the receptor interface part, ligand interface part
    """
    atom_list=[]
    atom_info_list=[]
    with open(pdb_path,'r') as file:
        line = file.readline()               # call readline()
        while line[0:4]!='ATOM':
            line=file.readline()
        while line :
            if line.startswith('ATOM'):
                x,y,z=float(line[30:38]), float(line[38:46]), float(line[46:54])
                atom_type=line[13:16].strip()
                atom_list.append(line)
                atom_info_list.append([x,y,z,atom_type])
            line=file.readline()
    print("Extracting atoms : %d" %len(atom_list))
    final_atoms=Form_interface(atom_info_list,atom_list)
    interface_path=Write_Interface(final_atoms, pdb_path, ".interface")
    print(interface_path)
    return interface_path

@set_timeout(100000, after_timeout)

def Form_interface(atom_info_list, atom_list ,cut_off=10):
    cut_off=cut_off**2
    interface_index = set()
    for item1 in range(len(atom_info_list)):
        for item2 in range(item1+1,len(atom_info_list)):
            dist_squar=sum((atom_info_list[item1][k] - atom_info_list[item2][k]) ** 2 for k in range(3))
            if dist_squar <= cut_off:
                interface_index.add(item1)
                interface_index.add(item2)               
    final_atoms=[atom_list[index] for index in interface_index]
    print("After filtering the interface region atom : %d" % len(final_atoms))
    return final_atoms


def Write_Interface(line_list,pdb_path,ext_file):
    new_path=pdb_path[:-4]+ext_file
    with open(new_path,'w') as file:
        for line in line_list:
            #check residue in the common residue or not. If not, no write for this residue
            residue_type = line[17:20]
            if residue_type in RESIDUE_Forbidden_SET:
                continue
            file.write(line)
    return new_path
