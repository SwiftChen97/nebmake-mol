'''
interpolate between two mass center and corresponding spherical coordinates
'''
from src.center import cal_center
from src.pboundry import adjust_centroid
from src.pboundry import find_neighbor
from src.coordinates_convert import cartesian_spherical_converter
import numpy as np

def sphere_interp(ele_i, ele_f, i_vec, f_vec, matrix, center_interp):
    ele_i = np.dot(find_neighbor(i_vec, ele_i), matrix)
    ele_f = np.dot(find_neighbor(f_vec, ele_f), matrix)

    r_ele_i, theta_ele_i, phi_ele_i = cartesian_spherical_converter(ele_i - np.dot(i_vec, matrix), 'spherical')
    r_ele_f, theta_ele_f, phi_ele_f = cartesian_spherical_converter(ele_f - np.dot(f_vec, matrix), 'spherical')
    interp_theta = np.linspace(theta_ele_i, theta_ele_f, 6)
    interp_phi = np.linspace(phi_ele_i, phi_ele_f, 6)

    interp_result = []
    for a in range(len(interp_theta)):
        spherical_ele = [r_ele_i, interp_theta[a], interp_phi[a]]
        x, y, z = cartesian_spherical_converter(spherical_ele, 'cartesian')
        direct_ele = np.dot(np.array([x, y, z]), imatrix)
        direct_ele = direct_ele + center_interp[a]
        interp_result.append(direct_ele)
    
    return interp_result 


centroid, matrix, imatrix = cal_center()
initial = centroid[0]
final = centroid[1]

# Linear interpolation for mass center
carbon = []
for idx in range(len(initial)):
    
    coord_interp = []
    center_interp = []

    initial_vec = initial[idx][0]
    final_vec = final[idx][0]
    # adjust the position of mass center based on periodic condition for continuous interpolation
    initial_vec, final_vec = adjust_centroid(initial_vec, final_vec)

    for seq in range(3):
        interp = np.linspace(initial_vec[seq], final_vec[seq], 6)
        coord_interp.append(interp)


    for num in range(len(coord_interp[0])):
        coordinate = [coord_interp[0][num], coord_interp[1][num], coord_interp[2][num]]
        center_interp.append(coordinate)

    # get every atom of organic molecule
    atom_pos_i = initial[idx][1]
    atom_pos_f = final[idx][1]
    for i in range(len(atom_pos_i)):
        atom_pos_i[i], atom_pos_f[i] = adjust_centroid(atom_pos_i[i], atom_pos_f[i])
    
    carbon_i = atom_pos_i[6]
    carbon_f = atom_pos_f[6]

    carbon_coord = sphere_interp(carbon_i, carbon_f, initial_vec, final_vec, matrix, center_interp)    
    carbon.append(carbon_coord)

carbon = list(map(list, zip(*carbon)))

