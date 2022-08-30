from src.operate import linear_interpolate
from utils import mv_boundry
from utils import reset_lattice
from src.operate import write_molecule, combine_poscar

def run():
    # consider the move direction of organic molecule
    func_dict = {
        '1': {'title': "consider the direction of movement of molecule", 'func': mv_boundry.run},
        '2': {'title': "extract molecule", 'func': write_molecule.run},
        '3': {'title': "interp inorganic atoms", 'func': linear_interpolate.run},
        '4': {'title': "interp organic atoms and combine all atoms together", 'func': combine_poscar.run},
        '5': {'title': "adjust lattice of interpolated structure", 'func': reset_lattice.run},
    }

    for choice in range(5):
        choice_dict = func_dict.get(str(choice + 1))
        choice_dict['func']()