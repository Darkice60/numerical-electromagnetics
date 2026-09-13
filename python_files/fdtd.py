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

# initialize an np array of the values of a function across a givem domain
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

#calculate the number of steps we will take in time
num_steps = math.ceil(total_time / delta_t)

# declare the nparrays for each compnenet of the elctric field and magnetic field with size based on the Yee Lattice
e_x = np.zeros((N_x, N_y + 1, N_z + 1), dtype=float)
e_y = np.zeros((N_x + 1, N_y, N_z + 1), dtype=float)
e_z = np.zeros((N_x + 1, N_y + 1, N_z), dtype=float)

h_x = np.zeros((N_x + 1, N_y, N_z), dtype=float)
h_y = np.zeros((N_x, N_y + 1, N_z), dtype=float)
h_z = np.zeros((N_x, N_y, N_z + 1), dtype=float)

# initialize the value of H^(n-1/2) and E^n through init_comp and the user enetered intila function
# the 0.5 addtions are again based on the Yee lattice's staggering
# also copy the array into "old" versions for calculations
init_comp(h_x, h_n_min_half_func, 0, 0.5, 0.5, 0)
init_comp(h_y, h_n_min_half_func, 0.5, 0, 0.5, 1)
init_comp(h_z, h_n_min_half_func, 0.5, 0.5, 0, 2)
old_h_x = h_x.copy()
old_h_y = h_y.copy()
old_h_z = h_z.copy()
init_comp(e_x, e_n_func, 0.5, 0, 0, 0)
init_comp(e_y, e_n_func, 0, 0.5, 0, 1)
init_comp(e_z, e_n_func, 0, 0, 0.5, 2)
old_e_x = e_x.copy()
old_e_y = e_y.copy()
old_e_z = e_z.copy()

# ideclare the nparrays for each compnenet of the elctric field and magnetic field with size based on the Yee Lattice
# while having the time at which the field is calculated as the primary index
e_x_hist = np.zeros((num_steps, N_x, N_y + 1, N_z + 1), dtype=float)
e_y_hist = np.zeros((num_steps, N_x + 1, N_y, N_z + 1), dtype=float)
e_z_hist = np.zeros((num_steps, N_x + 1, N_y + 1, N_z), dtype=float)

h_x_hist = np.zeros((num_steps, N_x + 1, N_y, N_z), dtype=float)
h_y_hist = np.zeros((num_steps, N_x, N_y + 1, N_z), dtype=float)
h_z_hist = np.zeros((num_steps, N_x, N_y, N_z + 1), dtype=float)

# for lop that runs for the number of steps needed to reach the total time form t = 0
for n in range(num_steps):
    ## SUMMARY FOR HOW ALL 6 LOOPS WORK
    # essentially this is just yee's method but long, in python and inefficient, and extremely hs student
    # basically we update the h fields and e field using numerical approximations of in matter where J_f = 0, D = Epsilon_0 E, and H = MU_0 B
    # namely we use faraday's law of induction and The ampere-maxwell law, which yields 6 equations that with rearragment yield favourable rsults
    # the array indcies are once again based of of the Yee lattice
    # (∂B/∂t)_x = -(∇×E)_x -> (H_x)^(n+1/2)[i, j, k] = (H_x)^(n-1/2)[i, j, k] - (Δt/μ₀) * ([(E_z)^n[i, j + 1, k] - (E_z)^n[i, j, k]]/Δy - [(E_y)^n[i, j, k + 1] - (E_y)^n[i, j, k]]/Δz)
    # (∂B/∂t)_y = -(∇×E)_y -> (H_y)^(n+1/2)[i, j, k] = (H_y)^(n-1/2)[i, j, k] - (Δt/μ₀) * ([(E_x)^n[i, j, k + 1] - (E_x)^n[i, j, k]]/Δz - [(E_z)^n[i + 1, j, k] - (E_z)^n[i, j, k]]/Δx)
    # (∂B/∂t)_z = -(∇×E)_z -> (H_z)^(n+1/2)[i, j, k] = (H_z)^(n-1/2)[i, j, k] - (Δt/μ₀) * ([(E_y)^n[i + 1, j, k] - (E_y)^n[i, j, k]]/Δx - [(E_x)^n[i, j + 1, k] - (E_x)^n[i, j, k]]/Δy)
    # (∂E/∂t)_x = [(∇×H)_x]/ε₀ -> (E_x)^(n+1)[i, j, k] = (E_x)^(n)[i, j, k] + (Δt/ε₀) * ([(H_z)^(n+1/2)[i, j, k] - (H_z)^(n+1/2)[i, j - 1, z]]/Δy - [(H_y)^(n+1/2)[i, j, k] - (H_y)^(n+1/2)[i, j, k - 1]]/Δz)
    # (∂E/∂t)_y = [(∇×H)_y]/ε₀ -> (E_y)^(n+1)[i, j, k] = (E_y)^(n)[i, j, k] + (Δt/ε₀) * ([(H_x)^(n+1/2)[i, j, k] - (H_x)^(n+1/2)[i, j, z - 1]]/Δz - [(H_z)^(n+1/2)[i, j, k] - (H_z)^(n+1/2)[i - 1, j, k]]/Δx)
    # (∂E/∂t)_z = [(∇×H)_z]/ε₀ -> (E_z)^(n+1)[i, j, k] = (E_z)^(n)[i, j, k] + (Δt/ε₀) * ([(H_y)^(n+1/2)[i, j, k] - (H_y)^(n+1/2)[i - 1, j, z]]/Δx - [(H_x)^(n+1/2)[i, j, k] - (H_x)^(n+1/2)[i, j - 1, k]]/Δy)
    # the rnage of all 6 loops are based on the size of the nparray and the operation required (i.e the E calculations are from 1 to N b/c we subtract one ansd dont want a -1 index)

    for i in range(h_x.shape[0]):
        for j in range(h_x.shape[1]):
            for k in range(h_x.shape[2]):
                h_x[i, j, k] = old_h_x[i, j, k] - (delta_t/MU_0) * ((e_z[i, j + 1, k] - e_z[i, j, k])/delta_y - (e_y[i, j, k + 1] - e_y[i, j, k])/delta_z)
    for i in range(h_y.shape[0]):
        for j in range(h_y.shape[1]):
            for k in range(h_y.shape[2]):
                h_y[i, j, k] = old_h_y[i, j, k] - (delta_t/MU_0) * ((e_x[i, j, k + 1] - e_x[i, j, k])/delta_z - (e_z[i + 1, j, k] - e_z[i, j, k])/delta_x)
    for i in range(h_z.shape[0]):
        for j in range(h_z.shape[1]):
            for k in range(h_z.shape[2]):
                h_z[i, j, k] = old_h_z[i, j, k] - (delta_t/MU_0) * ((e_y[i + 1, j, k] - e_y[i, j, k])/delta_x - (e_x[i, j + 1, k] - e_x[i, j, k])/delta_y)
    for i in range(e_x.shape[0]):
        for j in range(1, N_y):
            for k in range(1, N_z):
                e_x[i, j, k] = old_e_x[i, j, k] + (delta_t/EPSILON_0) * ((h_z[i, j, k] - h_z[i, j - 1, k])/delta_y - (h_y[i, j, k] - h_y[i, j, k - 1])/delta_z)
    for i in range(1, N_x):
        for j in range(e_y.shape[1]):
            for k in range(1, N_z):
                e_y[i, j, k] = old_e_y[i, j, k] + (delta_t/EPSILON_0) * ((h_x[i, j, k] - h_x[i, j, k - 1])/delta_z - (h_z[i, j, k] - h_z[i - 1, j, k])/delta_x)
    for i in range(1, N_x):
        for j in range(1, N_y):
            for k in range(e_z.shape[2]):
                e_z[i, j, k] = old_e_z[i, j, k] + (delta_t/EPSILON_0) * ((h_y[i, j, k] - h_y[i - 1, j, k])/delta_x - (h_x[i, j, k] - h_x[i, j - 1, k])/delta_y)

    # copy the array's into old array for the next loop
    old_h_x = h_x.copy()
    old_h_y = h_y.copy()
    old_h_z = h_z.copy()
    old_e_x = e_x.copy()
    old_e_y = e_y.copy()
    old_e_z = e_z.copy()

    # save the value of each componenet to the array's that track the value of each at each time step
    e_x_hist[n] = e_x
    e_y_hist[n] = e_y
    e_z_hist[n] = e_z
    h_x_hist[n] = h_x
    h_y_hist[n] = h_y
    h_z_hist[n] = h_z


# while not broken
while True:
    # ask the user what x-index they would like to see the arry it
    x_index = input("What x index would you like to see the fields evolve at? (Type \"esc\" to end)")
    # if the x index entered equals esc end the program by breaking the loop
    if (x_index.lower() == "esc"):
        break
    # ask the user what y and z index they would like to see the fields evolve at
    y_index = input("What y index would you like to see the fields evolve at?")
    z_index = input("What z index would you like to see the fields evolve at?")

    # form a try exception to insure the indecies are inetgers
    try:
        x_index = int(x_index)
        y_index = int(y_index)
        z_index = int(z_index)
    except ValueError:
        print("Please enter an integer.")
        continue

    # if the indices are out of bounds, make the user re eneter 
    if not (0 <= x_index < N_x and 0 <= y_index < N_y and 0 <= z_index < N_z):
        print("Not in domain. Try again.")
        continue
    # print the values of the array's through time
    print("E_x: " + e_x_hist[:, x_index, y_index, z_index])
    print("E_y: " + e_y_hist[:, x_index, y_index, z_index])
    print("E_z: " + e_z_hist[:, x_index, y_index, z_index])
    print("H_x: " + h_x_hist[:, x_index, y_index, z_index])
    print("H_y: " + h_y_hist[:, x_index, y_index, z_index])
    print("H_z: " + h_z_hist[:, x_index, y_index, z_index])