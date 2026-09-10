# Darkice60
# fdtd file for numerical-electromagnetics

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
# ask user for calculation resolution of the rect domain
res_x = int(input("Resolution in x?\n"))
res_y = int(input("Resolution in y?\n"))
res_z = int(input("Resolution in z?\n"))
# ask user for total sim time
total_time = float(input("Total simulation time?\n"))
# calculate Delta x, y, z to approx E and h and attain Delta t
delta_x = domain_x/(res_x-1)
delta_y = domain_y/(res_y-1)
delta_z = domain_z/(res_z-1)
# calculate the number of voxels across x, y, z
N_x = math.ceil(domain_x / delta_x)
N_y = math.ceil(domain_y / delta_y)
N_z = math.ceil(domain_z / delta_z)
# by the Courant-Friedrichs-Lewy Condition, delta t must be less tha or equal to
# 1/(C * math.sqrt((1/(delta_x)**2) + (1/(delta_y)**2) + (1/(delta_z)**2)))
# we multiply by 0.95 to create a 5% safety barrier
delta_t = 1/(C * math.sqrt((1/(delta_x)**2) + (1/(delta_y)**2) + (1/(delta_z)**2))) * 0.95

e_x = np.zeros((N_x, N_y + 1, N_z + 1), dtype=float)
e_y = np.zeros((N_x + 1, N_y, N_z + 1), dtype=float)
e_z = np.zeros((N_x + 1, N_y + 1, N_z), dtype=float)

h_x = np.zeros((N_x + 1, N_y, N_z), dtype=float)
h_y = np.zeros((N_x, N_y + 1, N_z), dtype=float)
h_z = np.zeros((N_x + 1, N_y, N_z), dtype=float)

