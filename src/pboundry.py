import numpy as np

def find_neighbor(ary1, ary2):
    if ary2[0] - ary1[0] > 0.5:
        ary2[0] -= 1
    elif ary2[0] - ary1[0] < -0.5:
        ary2[0] += 1

    if ary2[1] - ary1[1] > 0.5:
        ary2[1] -= 1
    elif ary2[1] - ary1[1] < -0.5:
        ary2[1] += 1

    if ary2[2] - ary1[2] > 0.5:
        ary2[2] -= 1
    elif ary2[2] - ary1[2] < -0.5:
        ary2[2] += 1
    return ary2

def neighbors(ini, coord_list, matrix, imatrix, length):
    '''
    find the neigbors of one atom based on the length between the atom and its neighbors

    Args:
        ini: the coordinate of the atom (Direct coordinate).
        coord_list: the coordinate of other atoms (Direct coordinate).
        matrix: the matrix of lattice basis vectors.
        imatrix: the inverse matrix of lattice basis vectors.
        length: the length distance the atom and its neighbors.

    Returns:
        List of neighbors
    '''
    neighbors = []
    for ary in coord_list:
        ary = np.dot(ary, imatrix)
        oary = find_neighbor(ini, ary)
        distance = np.sqrt(((np.dot(oary, matrix) - 
                    np.dot(ini, matrix)) ** 2).sum())
        if distance > length:
            continue
        neighbors.append(ary)
    return neighbors