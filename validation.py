import sys
import os
import numpy as np

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path))
print(current_path)

from diffusion import Props, Equation

time = 1.  # calculation time in sec
time_step = 0.1  # time step in sec
x_coord_in = 0.1  # length (m) is the length of the inner boundary
x_coord_out = 1.1  # length (m) is the length of the outer boundary
len_y = 1.  # length of grid block in y direction
len_z = 1.  # length of grid_block in z direction
grid_block_n = 5  # gridBlockN is a number of grid blocks
conc_left = 2.0  # outer concentration (density) (kg/m3)
conc_init = 10.0  # outer concentration (kg/m3)
diffusivity = 5.e-2  # diffusion coefficient (m2/sec)
it_accuracy = 1.e-20  # accuracy for iterative procedures

params = {"time": float(time),
          "time_step": float(time_step),
          "x_coord_in": float(x_coord_in),
          "x_coord_out": float(x_coord_out),
          "len_y": float(len_y),
          "len_z": float(len_z),
          "grid_block_n": int(grid_block_n),
          "conc_left": float(conc_left),
          "conc_init": float(conc_init),
          "diffusivity": float(diffusivity),
          "density": float(conc_left),
          "it_accuracy": float(it_accuracy)}

props = Props(params)
equation = Equation(params)
equation.calculate()

print()
print(equation.times)
print()
print(equation.velocities)
print()
print(equation.concs)


