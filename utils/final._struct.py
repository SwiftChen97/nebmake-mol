import os
import sys
dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(dir)
import config
from utils.reset_lattice import reset_final, delete_symbol, average


relative = average(delete_symbol(config.pos_ini)[1]) - \
			average(delete_symbol(config.pos_fin)[1])
print(relative)

para = list(delete_symbol(config.pos_fin))
para.append(relative)
para.append(config.pos_fin)

if __name__ == '__main__':
    reset_final(*para)

