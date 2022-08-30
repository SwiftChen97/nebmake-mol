import os
import sys
import config

def run():
    with open(os.path.join(config.base_dir, sys.argv[1]), 'r', encoding = 'utf-8') as file_obj1:
        data = file_obj1.readlines()
        ini_coord = data[8: config.total]
        head = data[:8]
        end = data[config.total + 1: ]
    
    with open(os.path.join(config.base_dir, sys.argv[2]), 'r', encoding = 'utf-8') as file_obj2:
        fin_coord = file_obj2.readlines()[8: config.total]
    
    with open(os.path.join(config.base_dir, sys.argv[1]), 'w', encoding = 'utf-8') as file_object:
        # write head info and atomic coordinates
        for info in head:
            file_object.write(info)

        for seq in range(len(ini_coord)):
            ini_coord[seq] = list(map(float, ini_coord[seq].strip().split()))
            fin_coord[seq] = list(map(float, fin_coord[seq].strip().split()))

            for num in range(3):
                if ini_coord[seq][num] - fin_coord[seq][num] > 0.5:
                    ini_coord[seq][num] -= 1.0
                elif ini_coord[seq][num] - fin_coord[seq][num] < -0.5:
                    ini_coord[seq][num] += 1.0

            # write coordinates info
            line = '  '.join(['%.10f' %k for k in ini_coord[seq]])
            file_object.write('{}\n'.format(line))
        
        file_object.write('\n')
        for end_line in end:
            file_object.write(end_line)