import os
import sys
import shutil
import config
from src.operate import write_molecule, linear_interpolate
from src.interp_process import interp
from src.operate.combine_poscar import extract_organic, extract_inorganic, inorganic_part, organic_part

def rot_mol():
    ini_path = os.path.join(config.base_dir, 'interp_result', 'interp_1.vasp')
    fin_path = os.path.join(config.base_dir, 'rot_fin.vasp')
    new_ini = os.path.join(config.base_dir, 'rot_ini.vasp')
    extract_inorganic(ini_path, fin_path)
    shutil.copyfile(ini_path, new_ini)
    rot = extract_organic()[3]

    with open(fin_path, 'a', encoding = 'utf-8') as file_object:
        for coor in rot[config.idx_seq[0][0]: config.idx_seq[0][1]]:
            file_object.write(f'{coor[0]:.6f} {coor[1]:.6f} {coor[2]:.6f}\n')
        for coor in rot[config.idx_seq[1][0]: config.idx_seq[1][1]]:
            file_object.write(f'{coor[0]:.6f} {coor[1]:.6f} {coor[2]:.6f}\n')
        for coor in rot[config.idx_seq[2][0]: config.idx_seq[2][1]]:
            file_object.write(f'{coor[0]:.6f} {coor[1]:.6f} {coor[2]:.6f}\n')
    
    
    shutil.rmtree(os.path.join(config.base_dir, 'center_neighbor'))
    shutil.rmtree(os.path.join(config.base_dir, 'interp_result'))
    write_molecule.run(fin_path, 'rot-mol')
    linear_interpolate.run(new_ini, fin_path, 'rot-result', 'rot')

    folder_path = os.path.join(config.base_dir, 'rot-result')
    inorganic_part(folder_path, 'rot')
    organic_part(folder_path, 'rot-mol', 'rot')

    for num in range(int(sys.argv[3])):
        direc = os.path.join(config.base_dir, f'0{num}')
        shutil.rmtree(direc)