"""
interpolate
"""
import sys
import pandas as pd
import numpy as np
import config
from src.molecule_code import dof_quant, coor_recurr

def load_data(file):
	cols = ["x", "y", "z"]
	cols = list(map(str, cols))
	df = pd.read_csv(file, sep=" ", header = None)
	df.columns = cols
	
	mass_center = np.dot(df.iloc[0], config.matrix)
	coor_data  = np.array([])
	
	x_list = df['x'][1:] 
	y_list = df['y'][1:]
	z_list = df['z'][1:]
	n_atom = len(x_list)
	
	for x,y,z in zip (x_list, y_list, z_list):
		coor_data = np.append(coor_data, [x,y,z])
	coor_data = np.resize(coor_data,(n_atom,3))
	coor_data = np.array([np.dot(vec, config.matrix) for vec in coor_data])
	
	return mass_center, coor_data, n_atom

def judge_dof_quan(file_i, file_f):
	phi_i = dof_quant(*load_data(file_i))[1]
	phi_f = dof_quant(*load_data(file_f))[1]
	if phi_i > phi_f:
		return dof_quant(*load_data(file_i), ref = - np.pi / 2), \
		dof_quant(*load_data(file_f), ref = - np.pi / 2)
	else:
		return dof_quant(*load_data(file_i)), dof_quant(*load_data(file_f))

def interp(file_i, file_f):
	theta_i, phi_i, gamma_i, struc_data_i = judge_dof_quan(file_i, file_f)[0]
	theta_f, phi_f, gamma_f = judge_dof_quan(file_i, file_f)[1][: 3]

	# interploate x, y, z degree of freedom
	interp_x = np.linspace(load_data(file_i)[0][0], load_data(file_f)[0][0], int(sys.argv[3]))
	interp_y = np.linspace(load_data(file_i)[0][1], load_data(file_f)[0][1], int(sys.argv[3]))
	interp_z = np.linspace(load_data(file_i)[0][2], load_data(file_f)[0][2], int(sys.argv[3]))

	interp_mc = []
	for i in range(len(interp_x)):
		coordinate = [interp_x[i], interp_y[i], interp_z[i]]
		interp_mc.append(coordinate)

	# consider clockwise or counterclockwise closer
	if abs(theta_i - theta_f) > np.pi:
		if theta_i > theta_f:
			theta_i = theta_i - 2 * np.pi
		else:
			theta_f = theta_f - 2 * np.pi

	if abs(phi_i - phi_f) > np.pi / 2:
		if phi_i > phi_f:
			phi_i = phi_i - np.pi
		else:
			phi_f = phi_f - np.pi

	# interpolate theta, phi, gamma (the rotation degree of freedom of molecule)
	interp_theta = np.linspace(theta_i, theta_f, int(sys.argv[3]))
	interp_phi = np.linspace(phi_i, phi_f, int(sys.argv[3]))
	interp_gamma = np.linspace(gamma_i, gamma_f, int(sys.argv[3]))

	#To Do:
	# 1. interpolate the structure data

	first_list = []
	second_list = []
	third_list = []

	if phi_i > phi_f:
		for i in range(len(interp_x)):
			coor_list = coor_recurr(interp_mc[i], 
								interp_theta[i], 
								interp_phi[i],
								interp_gamma[i], 
								struc_data_i,load_data(file_i)[-1], ref = - np.pi / 2)

			# reproduce the H atom
			first = list(coor_list[config.idx_seq[0][0]: config.idx_seq[0][1]])
			for a in range(len(first)):
				first[a] = list(np.dot(first[a], config.imatrix))
				first[a] = ' '.join(['%.6f' %k for k in first[a]])

			# reproduce the C atom
			second = list(coor_list[config.idx_seq[1][0]: config.idx_seq[1][1]])
			for a in range(len(second)):
				second[a] = list(np.dot(second[a], config.imatrix))
				second[a] = ' '.join(['%.6f' %k for k in second[a]])

			# reproduce the N atom
			third = list(coor_list[config.idx_seq[2][0]: config.idx_seq[2][1]])
			for a in range(len(third)):
				third[a] = list(np.dot(third[a], config.imatrix))
				third[a] = ' '.join(['%.6f' %k for k in third[a]])

			first_list.append(first)
			second_list.append(second)
			third_list.append(third)
	else:
		for i in range(len(interp_x)):
			coor_list = coor_recurr(interp_mc[i], 
								interp_theta[i], 
								interp_phi[i],
								interp_gamma[i], 
								struc_data_i,load_data(file_i)[-1])

			# reproduce the H atom
			first = list(coor_list[config.idx_seq[0][0]: config.idx_seq[0][1]])
			for a in range(len(first)):
				first[a] = list(np.dot(first[a], config.imatrix))
				first[a] = ' '.join(['%.6f' %k for k in first[a]])

			# reproduce the C atom
			second = list(coor_list[config.idx_seq[1][0]: config.idx_seq[1][1]])
			for a in range(len(second)):
				second[a] = list(np.dot(second[a], config.imatrix))
				second[a] = ' '.join(['%.6f' %k for k in second[a]])

			# reproduce the N atom
			third = list(coor_list[config.idx_seq[2][0]: config.idx_seq[2][1]])
			for a in range(len(third)):
				third[a] = list(np.dot(third[a], config.imatrix))
				third[a] = ' '.join(['%.6f' %k for k in third[a]])

			first_list.append(first)
			second_list.append(second)
			third_list.append(third)
	return first_list, second_list, third_list