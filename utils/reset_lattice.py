import os
import sys
import numpy as np
import config

def delete_symbol(file_path):
	with open(file_path, 'r') as file_object:
		info = file_object.readlines()
	head = info[:8]
	coord = info[8:]

	for seq in range(len(coord)):
		coord[seq] = np.array(list(map(float, coord[seq].split())))
	return head, coord

def average(coord):
	return sum(coord)/ len(coord)

def reset_final(header, coord, relative, final_path):
	file_object = open(final_path, 'w')
	
	# write basic information and lattice constant.
	for line in header:
		file_object.write(line)
	# write coordinate of each atom.
	for seq in range(len(coord)):
		coord[seq] = (coord[seq] + relative).tolist()
		coord_str = ' '.join(map(str, ['%.6f' %i for i in coord[seq]]))
		coord_str = '{}\n'.format(coord_str)
		file_object.write(coord_str)

def run():
	folder = os.path.join(config.base_dir, 'interp_result')

	for num in range(2, int(sys.argv[3]) + 1):
		ini_path = os.path.join(folder, 'interp_1.vasp')
		fin_path = os.path.join(folder, f'interp_{num}.vasp')

		relative = average(delete_symbol(ini_path)[1]) - \
					average(delete_symbol(fin_path)[1])

		para = list(delete_symbol(fin_path))
		para.append(relative)
		para.append(fin_path)
		print(para[2])
		reset_final(*para)