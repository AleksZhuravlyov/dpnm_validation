# MIT License
# 
# Copyright (c) 2020 Aleksandr Zhuravlyov and Zakhar Lanets
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import sys
import os
import numpy as np

from foamfile import FoamFile

import json

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path))
print(current_path)

from diffusion import Props, Equation
from sgrid import Sgrid, save_files_collection_to_file

import convection

# Diffusion

time_equil = 0.005
time = 0.002  # calculation time in sec
time_step = 0.000002  # time step in sec
y_coord_in = 0.002  # length (m) is the length of the inner boundary
y_coord_out = 0.006  # length (m) is the length of the outer boundary
len_x = 0.01  # length of grid block in x direction
len_z = 0.01  # length of grid_block in z direction
grid_block_d_n = 100  # grid_block_n is a number of diffusion grid blocks
conc_left = 2.0  # outer concentration (density) (kg/m3)
conc_init = 2.001  # outer concentration (kg/m3)
density = 2.0  # density (kg/m3)
diffusivity = 5.e-3  # diffusion coefficient (m2/sec)
it_accuracy = 1.e-20  # accuracy for iterative procedures

params = {'time_equil': float(time_equil),
          'time': float(time),
          'time_step': float(time_step),
          'y_coord_in': float(y_coord_in),
          'y_coord_out': float(y_coord_out),
          'len_x': float(len_x),
          'len_z': float(len_z),
          'grid_block_d_n': int(grid_block_d_n),
          'conc_left': float(conc_left),
          'conc_init': float(conc_init),
          'diffusivity': float(diffusivity),
          'density': float(density),
          'it_accuracy': float(it_accuracy)}

equation = Equation(params)
equation.calculate()

points_dims = np.array([2, params['grid_block_d_n'] + 1, 2], dtype=np.int32)
points_origin = np.array([0., params['y_coord_in'], 0.], dtype=np.float)
len_y = (params['y_coord_out'] - params['y_coord_in']) / params[
    'grid_block_d_n']
spacing = np.array([params['len_x'], len_y, params['len_z']], dtype=np.float)
sgrid = Sgrid(points_dims, points_origin, spacing)

os.system('rm -r inOut/diffusion/*.vtu')
os.system('rm -r inOut/diffusion/*.pvd')
concs_dict = dict()
file_name = 'inOut/diffusion/collection.pvd'
files_names = list()
files_descriptions = list()
for i in range(len(equation.times)):
    sgrid.cells_arrays = {'conc': np.array(equation.concs[i])}
    files_names.append(str(round(equation.times[i], 8)) + '.vtu')
    files_descriptions.append(str(round(equation.times[i], 8)))
    sgrid.save('inOut/diffusion/' + files_names[i])
save_files_collection_to_file(file_name, files_names, files_descriptions)

print(equation.velocities)

# Convection

grid_block_x_n = 30  # grid_block_x_n is a x number of convection grid blocks
grid_block_y_n = 10  # grid_block_y_n is a y number of convection grid blocks
viscosity = 2.e-2 / params['density']  # viscosity (m2/sec)
p_inlet = 300007.0 / params['density']  # inlet pressure (m2/sec2)
p_outlet = 300000.0 / params['density']  # outlet pressure (m2/sec2)

params['grid_block_x_n'] = int(grid_block_x_n)
params['grid_block_y_n'] = int(grid_block_y_n)
params['viscosity'] = float(viscosity)
params['diffusivity'] = float(diffusivity)
params['p_inlet'] = float(p_inlet)
params['p_outlet'] = float(p_outlet)

# controlDict

file_name = 'inOut/convection/system/controlDict'
with FoamFile(file_name) as f:
    foam_content = f.read()
foam_content['endTime'] = params['time'] + params['time_equil']
foam_content['deltaT'] = params['time_step']
with FoamFile(file_name, "w", foam_class="dictionary") as f:
    f.write(foam_content)

# transportProperties

file_name = 'inOut/convection/constant/transportProperties'
with FoamFile(file_name) as f:
    foam_content = f.read()
foam_content['nu'] = '[0 2 -1 0 0 0 0] ' + str(params['viscosity'])
with FoamFile(file_name, "w", foam_class="dictionary") as f:
    f.write(foam_content)

# blockMeshDict

file_name = 'inOut/convection/system/blockMeshDict'
with FoamFile(file_name) as f:
    foam_content = f.read()

x_min = 0.0
x_max = params['len_x']
y_min = 0.0
y_max = params['y_coord_in']
z_min = 0.0
z_max = params['len_z']

foam_content['vertices'][0] = [x_min, y_min, z_min]
foam_content['vertices'][1] = [x_max, y_min, z_min]
foam_content['vertices'][2] = [x_min, y_max, z_min]
foam_content['vertices'][3] = [x_max, y_max, z_min]
foam_content['vertices'][4] = [x_min, y_min, z_max]
foam_content['vertices'][5] = [x_max, y_min, z_max]
foam_content['vertices'][6] = [x_min, y_max, z_max]
foam_content['vertices'][7] = [x_max, y_max, z_max]

foam_content['blocks'][2] = [params['grid_block_x_n'],
                             params['grid_block_y_n'], 1]

with FoamFile(file_name, "w", foam_class="dictionary") as f:
    f.write(foam_content)

# p

file_name = 'inOut/convection/0/p'
with FoamFile(file_name) as f:
    foam_content = f.read()

foam_content[
    'boundaryField']['inlet']['value'] = 'uniform ' + str(params['p_inlet'])
foam_content[
    'boundaryField']['outlet']['value'] = 'uniform ' + str(params['p_outlet'])

with FoamFile(file_name, "w", foam_class="dictionary") as f:
    f.write(foam_content)

os.system('blockMesh -case inOut/convection')

os.system('rm -r inOut/convection/[1-9]*')
os.system('rm -r inOut/convection/0.*')

equil_steps_n = int(params['time_equil'] / params['time_step'])
results = convection.calculate(equation.velocities, equation.times,
                               equil_steps_n)

os.system('open -a ParaView-5.8.0.app inOut/diffusion/collection.pvd ')
os.system('open -a ParaView-5.8.0.app inOut/convection/system/controlDict.foam')

print(equation.velocities)
print(results['U_release_av'])
print(results['U_outlet_av'])
print()
print(results['U_outlet_equil_av'])
print()

np.save('U_release_av', np.array(results['U_release_av']))
np.save('U_outlet_av', np.array(results['U_outlet_av']))
np.save('times', np.array(equation.times))


with open('params.json', 'w') as f:
    json.dump(params, f)
