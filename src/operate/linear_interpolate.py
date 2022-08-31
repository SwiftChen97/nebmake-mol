import os
import sys
import numpy as np
import pymatgen as pmg
import config

def interpolate2(start_structure, end_structure, ximages=np.linspace(0,1,11),
                interpolate_lattices=False, pbc=True, autosort_tol=0):
    # sourcery skip: low-code-quality
    """
    Interpolate between this structure and end_structure. Useful for
    construction of NEB inputs.

    Args:
        end_structure (Structure): structure to interpolate between this
            structure and end.
        nimages (int): No. of interpolation images. Defaults to 10 images.
        interpolate_lattices (bool): Whether to interpolate the lattices.
            Interpolates the lengths and angles (rather than the matrix)
            so orientation may be affected.
        pbc (bool): Whether to use periodic boundary conditions to find
            the shortest path between endpoints.
        autosort_tol (float): A distance tolerance in angstrom in
            which to automatically sort end_structure to match to the
            closest points in this particular structure. This is usually
            what you want in a NEB calculation. 0 implies no sorting.
            Otherwise, a 0.5 value usually works pretty well.

    Returns:
        List of interpolated structures. The starting and ending
        structures included as the first and last structures respectively.
        A total of (nimages + 1) structures are returned.
    """
    # Check length of structures
    if len(start_structure) != len(end_structure):
        raise ValueError("Structures have different lengths!")

    if not (interpolate_lattices or start_structure.lattice == end_structure.lattice):
        raise ValueError("Structures with different lattices!")

    # Check that both structures have the same species
    for i in range(len(start_structure)):
        if start_structure[i].species_and_occu != end_structure[i].species_and_occu:
            raise ValueError("Different species!\nStructure 1:\n" +
                             str(start_structure) + "\nStructure 2\n" +
                             str(end_structure))

    start_coords = np.array(start_structure.frac_coords)
    end_coords = np.array(end_structure.frac_coords)

    if autosort_tol:
        dist_matrix = start_structure.lattice.get_all_distances(start_coords,
                                                     end_coords)
        site_mappings = collections.defaultdict(list)
        unmapped_start_ind = []
        for i, row in enumerate(dist_matrix):
            ind = np.where(row < autosort_tol)[0]
            if len(ind) == 1:
                site_mappings[i].append(ind[0])
            else:
                unmapped_start_ind.append(i)

        if len(unmapped_start_ind) > 1:
            raise ValueError("Unable to reliably match structures "
                             "with auto_sort_tol = %f. unmapped indices "
                             "= %s" % (autosort_tol, unmapped_start_ind))

        sorted_end_coords = np.zeros_like(end_coords)
        matched = []
        for i, j in site_mappings.items():
            if len(j) > 1:
                raise ValueError("Unable to reliably match structures "
                                 "with auto_sort_tol = %f. More than one "
                                 "site match!" % autosort_tol)
            sorted_end_coords[i] = end_coords[j[0]]
            matched.append(j[0])

        if len(unmapped_start_ind) == 1:
            i = unmapped_start_ind[0]
            j = list(set(range(len(start_coords))).difference(matched))[0]
            sorted_end_coords[i] = end_coords[j]

        end_coords = sorted_end_coords

    vec = end_coords - start_coords
    if pbc:
        vec -= np.round(vec)
    sp = start_structure.species_and_occu
    structs = []

    if interpolate_lattices:
        # interpolate lattice matrices using polar decomposition
        from scipy.linalg import polar
        # u is unitary (rotation), p is stretch
        u, p = polar(np.dot(end_structure.lattice.matrix.T,
                            np.linalg.inv(start_structure.lattice.matrix.T)))
        lvec = p - np.identity(3)
        lstart = start_structure.lattice.matrix.T

    for x in ximages:
        if interpolate_lattices:
            l_a = np.dot(np.identity(3) + x* lvec, lstart).T
            l = pmg.Lattice(l_a)
        else:
            l = start_structure.lattice
        fcoords = start_coords + x  * vec
        structs.append(start_structure.__class__(l, sp, fcoords,
                                      site_properties=start_structure.site_properties))
    return structs

def get_Q(struct1, struct2):
    disp = []
    mass = []
    for isite in enumerate(struct1.sites):
        disp.append(isite[1].distance(struct2.sites[isite[0]]))
        mass.append(isite[1].specie.atomic_mass)
    disp = np.array(disp)
    mass = np.array(mass)
    for i, j in zip(mass, disp):
        print(i, j, i*j**2/np.sum(mass*disp**2)*100)
    # formula for the configuration coordinate diagram
    return np.sqrt(np.sum(mass*disp**2))


def run(path1=sys.argv[1], path2=sys.argv[2], end='interp_result', label='interp'):
    # create folder to store the interpolated structures
    folder_path = os.path.join(config.base_dir, end)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    q0 = pmg.Structure.from_file(os.path.join(config.base_dir, path1))
    q1 = pmg.Structure.from_file(os.path.join(config.base_dir, path2))
    dQ = get_Q(q0,q1)
    print(dQ)

    # how many structures to interpolate.
    xx = np.linspace(0, 1, int(sys.argv[3]))

    interp=interpolate2(q0, q1, ximages = xx, pbc = True, interpolate_lattices = True)

    for itr,istr in enumerate(interp):
        file_path = os.path.join(folder_path, f"{label}_{itr + 1}.vasp")
        istr.to(fmt = "poscar",filename = file_path)