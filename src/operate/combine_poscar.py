"""
integrate the inorganic atoms and organic atoms into one POSCAR
"""
import re
import os
import shutil
import sys
import config
from src.interp_process import interp

# extract the coordinates of inorganic atoms
def extract_inorganic(ini_path, fin_path):
	with open(ini_path, 'r', encoding = 'utf-8') as file_object:
		data = file_object.readlines()
		ele = data[5].strip().split()
		num = [int(k) for k in data[6].strip().split()]

		head = data[:8]
		head[5] = config.ele_name
		head[6] = config.ele_num
		coord = data[8:]

		# remove the element symbol at the end
		for seq in range(len(coord)):
			coord[seq] = coord[seq].strip()

			tail = re.search('[A-Z][a-z]', coord[seq]) \
					or re.search('[A-Z]', coord[seq])

			if not tail:
				coord[seq] = '{}\n'.format(coord[seq])
			else:
				coord[seq] = '{}\n'.format(coord[seq].strip(tail.group()))

		part_pos = []
		for scop in num:
			part = list(coord[:scop])
			part_pos.append(part)
			del coord[0: scop]
		element = config.inorganic

		lead_iodine = []
		for seq in range(len(ele)):
			for atom in element:
				if ele[seq] != atom:
					continue
				lead_iodine.append(part_pos[seq])

		with open(fin_path, 'w', encoding = 'utf-8') as file_obj:
			for info in head:
				file_obj.write(info)

			for tp in lead_iodine:
				for coor in tp:
					file_obj.write(coor)

def extract_organic(mol_path='center_neighbor'):
	list_first = []
	list_second = []
	list_third = []
	for file in range(config.molecule_number):
		# read every molecule information from molecule file
		folder = os.path.join(config.base_dir, mol_path)

		ini_path = os.path.join(folder, f'initial_{str(file)}.dat')
		fin_path = os.path.join(folder, f'final_{str(file)}.dat')

		first, second, third, rot = interp(ini_path, fin_path)
		list_first.append(first)
		list_second.append(second)
		list_third.append(third)

	list_first = list(map(list, zip(*list_first)))
	list_second = list(map(list, zip(*list_second)))
	list_third = list(map(list, zip(*list_third)))
	return list_first, list_second, list_third, rot

# write lattice constant and the coordinates of lead and iodine.
def inorganic_part(folder, label='interp'):
	for image in range(int(sys.argv[3])):
		# find structure file path and get content
		file_path = os.path.join(folder, f'{label}_{str(image + 1)}.vasp')
		extract_inorganic(file_path, file_path)
		

def organic_part(wfolder,path='center_neighbor', label='interp'):
	l1, l2, l3, rot = extract_organic(path)
	# write organic section
	for image in range(int(sys.argv[3])):
		file_path = os.path.join(wfolder, f'{label}_{str(image + 1)}.vasp')

		with open(file_path, 'a', encoding = 'utf-8') as file_object:
			for j in l1[image]:
				for k in range(len(j)):
					file_object.write('{}\n'.format(j[k]))
			for j in l2[image]:
				for k in range(len(j)):
					file_object.write('{}\n'.format(j[k]))
			for j in l3[image]:
				for k in range(len(j)):
					file_object.write('{}\n'.format(j[k]))

def run():
	folder_path = os.path.join(config.base_dir, 'interp_result')
	inorganic_part(folder_path)
	organic_part(folder_path)

	#construct NEB input directories
	for num in range(int(sys.argv[3])):
		new_direc = os.path.join(config.base_dir, f'0{num}')
		if not os.path.exists(new_direc):
			os.makedirs(new_direc)
		file_path = os.path.join(folder_path, f'interp_{str(num + 1)}.vasp')
		new_fpath = os.path.join(new_direc, 'POSCAR')
		shutil.copyfile(file_path, new_fpath)