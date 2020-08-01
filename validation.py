import sys
import os
import numpy as np

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path))
print(current_path)

from diffusion import Props, Equation
from sgrid import Sgrid, save_files_collection_to_file

time = 15.  # calculation time in sec
time_step = 0.1  # time step in sec
y_coord_in = 0.1  # length (m) is the length of the inner boundary
y_coord_out = 1.1  # length (m) is the length of the outer boundary
len_x = 20.  # length of grid block in x direction
len_z = 10.  # length of grid_block in z direction
grid_block_n = 100  # gridBlockN is a number of grid blocks
conc_left = 2.0  # outer concentration (density) (kg/m3)
conc_init = 10.0  # outer concentration (kg/m3)
diffusivity = 5.e-2  # diffusion coefficient (m2/sec)
it_accuracy = 1.e-20  # accuracy for iterative procedures

params = {'time': float(time),
          'time_step': float(time_step),
          'y_coord_in': float(y_coord_in),
          'y_coord_out': float(y_coord_out),
          'len_x': float(len_x),
          'len_z': float(len_z),
          'grid_block_n': int(grid_block_n),
          'conc_left': float(conc_left),
          'conc_init': float(conc_init),
          'diffusivity': float(diffusivity),
          'density': float(conc_left),
          'it_accuracy': float(it_accuracy)}

props = Props(params)
equation = Equation(params)
equation.calculate()


points_dims = np.array([2, params['grid_block_n'] + 1, 2], dtype=np.int32)
points_origin = np.array([0., params['y_coord_in'], 0.], dtype=np.float)
len_y = (params['y_coord_out'] - params['y_coord_in']) / params['grid_block_n']
spacing = np.array([params['len_x'], len_y, params['len_z']], dtype=np.float)

sgrid = Sgrid(points_dims, points_origin, spacing)

concs_dict = dict()
file_name = 'diffusion.pvd'
files_names = list()
files_descriptions = list()

for i in range(len(equation.times)):
    sgrid.cells_arrays = {'conc': np.array(equation.concs[i])}
    files_names.append('diffusion_' + str(i) + '.vtu')
    files_descriptions.append(str(round(equation.times[i], 3)))
    sgrid.save(files_names[i])

save_files_collection_to_file(file_name, files_names, files_descriptions)
