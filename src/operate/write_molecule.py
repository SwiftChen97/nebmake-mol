import os
import sys
import config
from src.center import cal_center

# get mass center of each state
def center_info(fin_path):
    """read POSCAR and return all of the mass center of molecules and neighbor atoms.

    Returns:
        _type_: center and neighbors.
    """
    center_list = []
    for state in [os.path.join(config.base_dir, sys.argv[1]), os.path.join(config.base_dir, fin_path)]:
        with open(state, 'r', encoding = 'utf-8') as file_object:
            info = file_object.readlines()[: config.total]
        parameter = [info, config.atom_mass, float(sys.argv[4])]
        center = cal_center(*parameter)
        center_list.append(center[0])
    return center_list

def write_text(folder, center_vec, state, seq, name):
    """
    Args:
        folder: the folder to store the files.
        vec: the vector of center of mass.
        state: the state before or after phase transition.
    """
    file_path = os.path.join(folder, f'{name}_{str(seq)}.dat')
    
    with open(file_path, 'w', encoding = 'utf-8') as file_object:
        center = '%.10f %.10f %.10f\n' %(center_vec[0], center_vec[1], center_vec[2])
        file_object.write(center)

        for line in state[seq][1]:
            line = '%.10f %.10f %.10f\n' %(line[0], line[1], line[2])
            file_object.write(line)

def run(path=sys.argv[2], folder='center_neighbor'):
    # create directories to store the files
    write_folder = os.path.join(config.base_dir, folder)
    if not os.path.exists(write_folder):
        os.makedirs(write_folder)
    
    # extract the center of mass of each state
    initial = center_info(path)[0]
    final = center_info(path)[1]
    
    # write files
    for ele in range(len(initial)):
        initial_center = initial[ele][0]
        final_center = final[ele][0]
    
        write_text(write_folder, initial_center, initial, ele, name = 'initial')
        write_text(write_folder, final_center, final, ele, name = 'final')