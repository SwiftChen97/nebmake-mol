import os
import sys
from src.center import cal_center

organic_info = {'C':12.011, 'N':14.007, 'H':1.008}
base_dir =  os.getcwd()

with open(os.path.join(base_dir, sys.argv[1]), 'r', encoding = 'utf-8') as file_object:
    info = file_object.readlines()
    info = [line.strip() for line in info]
    values = info[6].split()
    total = sum(int(num) for num in values) + 8
    info = info[: total]
atom_mass = {ele: organic_info[ele] for ele in info[5].split() if ele in organic_info}
result = cal_center(info, atom_mass, float(sys.argv[4]))
matrix, imatrix, molecule_number, idx_seq, max_a, min_a = result[1], result[2], result[3], result[4], result[5], result[6]

keys = info[5].split()
name_num = dict(zip(keys, values))
# basic information of interpolate
inorganic = ['Pb', 'I']
organic = ['C', 'N', 'H']
ele_name = f"{'  '.join(inorganic)}  {'  '.join(organic)}\n"
ele_num = f"{'  '.join([name_num[ele] for ele in ele_name.strip().split()])}\n"