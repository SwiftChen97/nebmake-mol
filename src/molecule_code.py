import pandas as pd
import numpy as np
import config
from src.coordinates_convert import cartesian_spherical_converter

def count_mass_center(coor_data, mass_data):
	x = np.dot(coor_data[:,0], mass_data)
	y = np.dot(coor_data[:,1], mass_data)
	z = np.dot(coor_data[:,2], mass_data)
	tot_m = np.sum(mass_data)

	x /= tot_m
	y /= tot_m
	z /= tot_m
	return np.array([x, y, z])

def deg_calcu(vector_c, vector_e1, vector_e2):
	cos1 = np.dot(vector_c, vector_e1) / np.linalg.norm(vector_c)
	cos2 = np.dot(vector_c, vector_e2) / np.linalg.norm(vector_c)
	cos1 = np.around(cos1, 4)
	gamma = np.arccos(cos1)
	
	if cos2 < 0:
		gamma = 2 * np.pi -gamma
	return gamma

def dof_quant(mass_center, coor_data, n_atom, ref = np.pi / 2):
	print(coor_data)
	vector_a = coor_data[config.min_a + 1] - mass_center
	vector_b = coor_data[config.max_a - 1] - mass_center
	r1, theta, phi = cartesian_spherical_converter(vector_a, 'spherical')
	
	vector_e1 = np.array([1.0, theta, phi + ref]) 
	xe1, ye1, ze1 = cartesian_spherical_converter(vector_e1, 'cartesian')
	vector_e1 = np.array([xe1, ye1, ze1])
	
	vector_e2 = np.cross(vector_a, vector_e1) / np.linalg.norm(vector_a)
	
	vector_c = vector_b - vector_a * np.dot(vector_b, vector_a) / np.dot(
		vector_a,vector_a)
	
	gamma = deg_calcu(vector_c, vector_e1, vector_e2)
	
	vector_g1 = vector_a / np.linalg.norm(vector_a)
	vector_g2 = vector_c / np.linalg.norm(vector_c)
	vector_g3 = np.cross(vector_g1, vector_g2)

	struc_data = pd.DataFrame(columns=('l1', 'l2', 'l3'))
	for i in range(n_atom):
		vector_i = coor_data[i] - mass_center
		l1 = np.dot(vector_i, vector_g1)
		l2 = np.dot(vector_i, vector_g2)
		l3 = np.dot(vector_i, vector_g3)
		struc_data.loc[i] = [l1,l2,l3]
		
	return theta, phi, gamma, struc_data

def coor_recurr(mass_center, theta, phi, gamma, struc_data, n_atom, ref = np.pi / 2):
	vector_g1 = np.array([1.0, theta, phi])
	xg1, yg1, zg1 = cartesian_spherical_converter(vector_g1,'cartesian')
	vector_g1 = np.array([xg1, yg1, zg1])
	
	vector_e1 = np.array([1.0, theta, phi + ref])
	xe1, ye1, ze1 = cartesian_spherical_converter(vector_e1,'cartesian')
	vector_e1 = np.array([xe1, ye1, ze1])
	
	vector_e2 = np.cross(vector_g1, vector_e1)
	vector_g2 = vector_e1* np.cos(gamma) + vector_e2* np.sin(gamma)
	
	vector_g3 = np.cross(vector_g1, vector_g2)
	
	coor_data = np.array([])
	for i in range(n_atom):
		coor_i = mass_center 
		coor_i += struc_data['l1'][i]* vector_g1 
		coor_i += struc_data['l2'][i]* vector_g2
		coor_i += struc_data['l3'][i]* vector_g3
		coor_data = np.append(coor_data, coor_i)
	
	coor_data = np.resize(coor_data, (n_atom, 3))
	return coor_data