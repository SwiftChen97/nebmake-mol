import numpy as np
from src.pboundry import neighbors

def cal_center(info, atom_mass, length):
	"""This function returns the center of mass of the structure.

	Args:
		info (_type_): readlines POSCAR file. (list)
		atom_mass (_type_): element name and mass. (dict)
		length (_type_): the distance in where find all atoms of one molecule. (float)

	Returns:
		_type_: center_list, matrix and imatrix. (list, ndarray, ndarray)
	"""
# calculate basis vector of 3 directions
	matrix = np.array([list(map(float, info[2].split())), 
						list(map(float, info[3].split())), 
						list(map(float, info[4].split()))]) \
					* float(info[1])
	imatrix = np.linalg.inv(matrix)

	atom_type = info[5].split()
	atom_num = [int(i) for i in info[6].split()]
	coord = info[8: ]
	coord = [list(map(float, cord.split())) for cord in coord]
	coord = [np.dot(lst, matrix) for lst in coord]
# extract the coordinates of atomic coordinate of organic molecule
	organic_unit = []
	for number in atom_num:
		each = [np.array(i) for i in coord[:number]]
		organic_unit.append(each)
		del coord[0: number]
	for symbol in range(len(atom_type) - 1, -1, -1):
		if atom_type[symbol] in atom_mass.keys():
			continue
		del organic_unit[symbol]

# find which atom have the least quantities and set the coordinate of this atom as standard
	numbers = [len(i) for i in organic_unit]
	idx = numbers.index(min(numbers))
	std_atom = organic_unit[idx]

# find the neighbor atom of the standard atom and transform back to direct coordinate
	center_list = []
	sort = []
	for xyz in std_atom:
		std_atomass = atom_mass.get(list(atom_mass.keys())[idx])
		xyz = np.dot(xyz, imatrix)
		molecule = []

# the mass and coordinate of standard atom 
		molecule_vecmass = xyz * std_atomass
		total_mass = std_atomass
		for other in range(len(organic_unit)):
			para = (xyz, organic_unit[other], 
					matrix, imatrix, length)
			nbors = neighbors(*para)
			sort.append(len(nbors))
			molecule.extend(nbors)
# get the value of coordinate multiply atomic mass of each molecule
			ele_symbol = list(atom_mass.keys())[other]
			atomass = atom_mass.get(ele_symbol)
			vec_mass = np.dot(nbors, atomass)
			molecule_vecmass += sum(vec_mass)
			total_mass += atomass * len(nbors)
# calculate coordinate of mass center
		center = molecule_vecmass / total_mass
		center_neighbor = [center, molecule]
		center_list.append(center_neighbor)

	new_sort = sort[0: len(organic_unit)]
	ini_num = 0
	idx_aran = [0]

	for index,ele in enumerate(new_sort):
		if index == len(new_sort) - 1:
			ini_num += ele
			idx_aran.append(ini_num)
		
		else:
			ini_num +=ele
			idx_aran.extend([ini_num, ini_num])

	idx_aran = [idx_aran[i:i+2] for i in range(0,len(idx_aran),2)]
	max_seq = new_sort.index(max(new_sort))
	min_seq = new_sort.index(min(new_sort))


	atom_more = idx_aran[max_seq][1]
	atom_less = idx_aran[min_seq][0]

	return center_list, matrix, imatrix, len(std_atom), idx_aran, atom_more, atom_less