import sys
import os
import numpy as np

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path))
print(current_path)

from diffusion import Props, Equation

# calculation time in sec
time = 1.
# time step in sec
time_step = 0.1
# length (m) is the length of the inner boundary
x_coord_in = 0.1
# length (m) is the length of the outer boundary
x_coord_out = 1.1
# length of grid block in y direction
len_y = 1.
# length of grid_block in z direction
len_z = 1.
# gridBlockN is a number of grid blocks
grid_block_n = 5
# outer concentration (kg/m3)
conc_left = 2.0
# outer concentration (kg/m3)
conc_right = 10.0
# diffusion coefficient (m2/sec)
diffusivity = 5.e-2
# iterative_accuracy is the accuracy for iterative procedures
it_accuracy = 1.e-20

params = {"time": float(time),
          "time_step": float(time_step),
          "x_coord_in": float(x_coord_in),
          "x_coord_out": float(x_coord_out),
          "len_y": float(len_y),
          "len_z": float(len_z),
          "grid_block_n": int(grid_block_n),
          "conc_left": float(conc_left),
          "conc_right": float(conc_right),
          "diffusivity": float(diffusivity),
          "it_accuracy": float(it_accuracy)}

print(params)
props = Props(params)
equation_diffusion = Equation(params)
# props.print_params()
