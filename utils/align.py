import operator
import numpy as np
import itertools
import config
from src.operate.write_molecule import center_info
from src.molecule_code import dof_quant, coor_recurr

def coor_move1(coord,interp,arg_fin):
	for i in range(3):
		for j in range(interp[i]):
			nj = arg_fin[i][j]
			coord[nj][i] += 1.0

def atom_align(coord_ini, coord_fin, coor_abc, count_a_bool, fixed_a = None):
	if fixed_a is None:
		fixed_a = [0.,0.,0.]
	the_n = len(coord_fin)
	nlist = np.arange(the_n)

	nums = itertools.permutations(nlist)
	ncoms = list(nums)
	ncoms = np.array(ncoms)

	xinterp, yinterp, zinterp = np.meshgrid(nlist,nlist,nlist,indexing='ij')
	interps = [[xc, yc, zc]
	           for xc, yc, zc in zip(xinterp.flatten(), yinterp.flatten(),
	                                 zinterp.flatten())]
	interps = np.array(interps)

	co_ini = []
	for cini in coord_ini:
		co_ini_t = np.dot(cini,coor_abc)
		co_ini.append(co_ini_t)
	co_ini = np.array(co_ini)

	arg_fin = []
	for i in range(3):
		arg_fin_t = coord_fin[:, i].argsort()
		arg_fin.append(arg_fin_t)
	arg_fin = np.array(arg_fin)

	interp_set = []
	dis_set = []
	a_set = []

	for ncom in ncoms:
		dis_set_t = []
		a_set_t = []
		for interp in interps:
			coord_fin_t = coord_fin.copy()
			coor_move1(coord_fin_t,interp,arg_fin)

			co_fin = []
			for cfin in coord_fin:
				co_fin_t = np.dot(cfin,coor_abc)
				co_fin.append(co_fin_t)
			co_fin = np.array(co_fin)

			drs = []
			for i in nlist:
				ni = ncom[i]
				dr = co_fin[ni] - co_ini[i]
				drs.append(dr)
			drs = np.array(drs)

			if count_a_bool:
				vector_a = np.sum(drs, axis = 0)/(-the_n)
			else:
				vector_a = np.array([0,0,0])

			for i in nlist:
				drs[i] = drs[i] + vector_a

			dis_t = np.sum(drs**2)

			dis_set_t.append(dis_t)
			a_set_t.append(vector_a)

		dis_set_t = np.array(dis_set_t)
		arg_dis_t = dis_set_t.argsort()
		n0 = arg_dis_t[0]

		interp_set.append(n0)
		dis_set.append(dis_set_t[n0])
		a_set.append(a_set_t[n0])

	dis_set = np.array(dis_set)
	arg_dis = dis_set.argsort()
	n00 = arg_dis[0]

	best_com = ncoms[n00]
	best_interp = interps[interp_set[n00]]
	best_a = a_set[n00]

	co_fin_out = coord_fin.copy()
	coor_move1(co_fin_out,best_interp,arg_fin)
	co_fin_out = co_fin_out[best_com]

	imatrix = np.linalg.inv(coor_abc)
	best_a = np.dot(best_a, imatrix)

	if count_a_bool:
		imatrix = np.linalg.inv(coor_abc)
		best_a = np.dot(best_a, imatrix)
	else:
		best_a = np.array(fixed_a)

	for i in nlist:
		co_fin_out[i] = co_fin_out[i] + best_a

	return co_fin_out, best_a, best_com

def inorganic():
	pass

def align_one_hydrogen(center_atom):
	coord_list = []
	for mol in center_atom:
		dis_seq = []
		hydrogen = mol[1][config.hydrogen_begin: config.hydrogen_end]
		carbon = np.dot(mol[1][config.carbon_begin], config.matrix)
		nitrogen = mol[1][config.nitrogen_begin]
		for idx, h in enumerate(hydrogen, 0):
			h = np.dot(h, config.matrix)
			dis = np.sqrt(np.sum((h - carbon) ** 2))
			dis_seq.append([dis, idx])
		dis_seq.sort(key=operator.itemgetter(0))
		new_seq = [list[1] for list in dis_seq]
		new_hydrogen = [hydrogen[seq] for seq in new_seq]
		new_hydrogen.insert(config.carbon_begin, np.dot(carbon, config.imatrix))
		new_hydrogen.insert(config.nitrogen_begin, nitrogen )
		coord_mol = new_hydrogen.copy()
		coord_list.append((mol[0],coord_mol))
	return coord_list

def organic():
	molecule_ini = center_info()[0]
	molecule_fin = center_info()[1]
	center_ini = np.array([mol[0] for mol in molecule_ini])
	center_fin = np.array([mol[0] for mol in molecule_fin])

	# adjust the sequence of molecule based on the initial structure
	# align the center of the molecules in the initial and final states.
	result = atom_align(center_ini, center_fin, config.matrix, True) # FALSE, input fix_a
	seq = result[-1]
	center_fin = result[0]
	new_fin = [molecule_fin[num] for num in seq]

	for seq in range(len(new_fin)):
		new_fin[seq][0] = center_fin[seq]

	ini_mol = align_one_hydrogen(molecule_ini)
	fin_mol = align_one_hydrogen(new_fin)

	for seq in range(len(ini_mol)):
		n_atom = len(ini_mol[seq][1])

		ini_center = ini_mol[seq][0]
		fin_center = fin_mol[seq][0]

		# the final molecule which center has been aligned with the initial molecule
		ini_coord = np.resize(ini_mol[seq][1], (n_atom,3))
		fin_coord = np.resize(fin_mol[seq][1], (n_atom,3))

		ini_para = (ini_center, ini_coord, n_atom)
		fin_para = (fin_center, fin_coord, n_atom)
		# get initial structure and degree of freedom of the final structure
		struc_data_i = np.array(dof_quant(*ini_para)[-1])
		theta_f, phi_f, gamma_f, struc_data_f = dof_quant(*fin_para)
		struc_data_f = np.array(struc_data_f)
		
		# extract other atoms which don't need to be aligned
		fin_carbon = struc_data_f[0]
		fin_hydro = struc_data_f[1]
		fin_nitrogen = struc_data_f[-1]

		fin_strt = atom_align(struc_data_i[2: 7], struc_data_f[2: 7], config.matrix, False)[0]
		# combing the aligned atoms and other atoms
		fin_strt = list(fin_strt)
		fin_strt.insert(0, fin_carbon)
		fin_strt.insert(1, fin_hydro)
		fin_strt.append(fin_nitrogen)

		# transfer the interpolated structure data to reccurence structure
		recur_para = (fin_center, theta_f, phi_f, gamma_f, fin_strt, n_atom)
		fin_mol = coor_recurr(*recur_para)

		# TO DO:
		# 1. write coordinates of inorganic atoms to POSCAR
		# 2. align inorganic atoms
