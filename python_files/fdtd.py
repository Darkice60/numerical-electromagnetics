# Darkice60
# fdtd file for numerical-electromagnetics

"""
TEST 1
E = (0,0,0)
H = (1,2,3)
Expected: no change

TEST 2
E = (z,0,y)
H = (y,z,x)
Expected:
Hx decreases by dt/mu0
Hy decreases by dt/mu0
Hz unchanged
"""

# imports
import numpy as np
import math
import matplotlib.pyplot as plt
from typing import Final

## Constants
EPSILON_0: Final = 8.8541878128e-12
MU_0: Final = 4 * math.pi * 1e-7
C: Final = 299_792_458

# create vector function at time t
def create_func(arr):
    def func(x, y, z):
        return np.array([eval(expression, {"x": x, "y": y, "z":z, "math": math, "np": np}) for expression in arr])
    return func

def init_comp(arr, func, x_half, y_half, z_half, dim):
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            for k in range(arr.shape[2]):
                arr[i, j, k] = func((i + x_half) * delta_x, (j + y_half) * delta_y, (k + z_half) * delta_z)[dim]
    return arr


# get the function E^n and H^(n-1/2) from the user and make them into a function
print("THE EQUATIONS MUST SATISFY GAUSS' LAW OF ELECTRICITY AND MAGNETISM")
e_n = input("Enter the Electric field (E) at time n (Seperate by commas)\n").split(",")
e_n_func = create_func(e_n)
h_n_min_half = input("Enter the Magnetic Field Strength (H) at time n - 1/2 (Seperate by commas)\n").split(",")
h_n_min_half_func = create_func(h_n_min_half)
# ask the user for the size of the rect domain
domain_x = float(input("Size of x domain?\n"))
domain_y = float(input("Size of y domain?\n"))
domain_z = float(input("Size of z domain?\n"))
# ask user for number of voxels
N_x = int(input("Voxels in x?\n"))
N_y = int(input("Voxels in y?\n"))
N_z = int(input("Voxels in z?\n"))
# ask user for total sim time
total_time = float(input("Total simulation time?\n"))
# calculate Delta x, y, z to approx E and h and attain Delta t
delta_x = domain_x/N_x
delta_y = domain_y/N_y
delta_z = domain_z/N_z
# by the Courant-Friedrichs-Lewy Condition, delta t must be less tha or equal to
# 1/(C * math.sqrt((1/(delta_x)**2) + (1/(delta_y)**2) + (1/(delta_z)**2)))
# we multiply by 0.95 to create a 5% safety barrier
delta_t = 1/(C * math.sqrt((1/(delta_x)**2) + (1/(delta_y)**2) + (1/(delta_z)**2))) * 0.95

e_x = np.zeros((N_x, N_y + 1, N_z + 1), dtype=float)
e_y = np.zeros((N_x + 1, N_y, N_z + 1), dtype=float)
e_z = np.zeros((N_x + 1, N_y + 1, N_z), dtype=float)

h_x = np.zeros((N_x + 1, N_y, N_z), dtype=float)
h_y = np.zeros((N_x, N_y + 1, N_z), dtype=float)
h_z = np.zeros((N_x, N_y, N_z + 1), dtype=float)

old_h_x = init_comp(h_x, h_n_min_half_func, 0, 0.5, 0.5, 0)
e_y = init_comp(e_y, e_n_func, 0, 0.5, 0, 1)
e_z = init_comp(e_z, e_n_func, 0, 0, 0.5, 2)

for i in range(h_x.shape[0]):
    for j in range(h_x.shape[1]):
        for k in range(h_x.shape[2]):
            h_x[i, j, k] = old_h_x[i, j, k] - (delta_t/MU_0) * ((e_z[i, j+1, k] - e_z[i, j, k])/delta_y - (e_y[i, j, k + 1] - e_y[i, j, k])/delta_z)
old_h_x = h_x.copy()


print(h_x)