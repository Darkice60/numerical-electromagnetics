# Darkice60
# main file for numerical-electromagnetics

# imports
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
from typing import Final

# user enters which calc they want to perform
choice = int(input("1. 1D Laplace\n2. 2D Laplace\n3. 2D Poisson\n4. 3D Poisson\n"))

# declare the electric constant for u
EPSILON_0: Final = 8.8541878128e-12

# functions builders to create the charge density function rho with 2 argumanets and 3 arguments
def create_rho_2d(rho_string):
    def rho(x, y):
        return eval(rho_string, {"x": x, "y": y, "math": math, "np": np})
    return rho

def create_rho_3d(rho_string):
    def rho(x, y, z):
        return eval(rho_string, {"x": x, "y": y, "z":z, "math": math, "np": np})
    return rho

# match case based on the users input
match choice:
    # 1D Lapcian Solver: ∇^2V = 0 -> ∂^2V/∂x^2 = 0 
    case 1:
        ## The boundary is a point, so the user must enter the initial and final electric potential
        initial = float(input("What is the initial electric potential\n"))
        final = float(input("What is the final electric potential\n"))
        # init np.array that represnts potential difference across the line
        v = np.array([initial, 0, 0, 0, 0, final], dtype=float)
        # init the counter accumulator to counter how many times the Finite-Difference Approximation is used (through the entire array, i.e. updated when all elements have been)
        c = 0
        # the loop runs until it is broken
        while True:
            # accumulate onto counter
            c += 1
            # make old_v equal to a copy of v before the FDA is applied so that an error bound can be checked
            old_v = v.copy()
            # use the FDA across the entire line. NOTE: EXCLUDING the boundary values
            for i in range(1, len(v) - 1):
                # apply the FDA at index i
                v[i] =(v[i-1] + v[i+1]) /2
            # find the difference between v vs. old_v for error bound
            diff = v - old_v
            # the highest value in diff is returned
            max_diff = np.max(np.abs(diff))
            # if the largest difference is less than 1e-10, break
            if (max_diff < 1e-10):
                break
        #print the np potential array and the loop runs required
        print(v)
        print(c)
    # 2D Lapcian Solver: ∇^2V = 0 -> ∂^2V/∂x^2 + ∂^2V/∂y^2 = 0
    case 2:
        ## The boundary is a square, so the user must enter the electric potential on the left, reight, top and bottom
        left = float(input("What is the minimum x electric potential\n"))
        right = float(input("What is the maximum x electric potential\n"))
        bottom = float(input("What is the minimum y in electric potential\n"))
        top = float(input("What is the maximum y electric potential\n"))
        # one must also ask the resolution of the calculation, for this solver it is the same as the domain
        # resolution is the number of of grid points in a row (so dim^2 is the total grid points)
        dim = int(input("How many rows and columns in the calculation grid?\n"))
        # init np.array that represnts potential difference across the square
        v = np.zeros((dim, dim), dtype=float)
        # set each side of the np.array equal to there given values, i.e. the left side of the array has the value of left
        v[:, 0] = left
        v[:, -1] = right
        v[-1, :] = bottom
        v[0, :] = top
        # init the counter accumulator to counter how many times the Finite-Difference Approximation is used (through the entire array, i.e. updated when all elements have been)
        c = 0
        # the loop runs until it is broken
        while True:
            # accumulate onto counter
            c += 1
            # make old_v equal to a copy of v before the FDA is applied so that an error bound can be checked
            old_v = v.copy()
            # use the FDA across the entire square. NOTE: EXCLUDING the boundary values
            for i in range(1, v.shape[0] - 1):
                for j in range(1, v.shape[1] - 1):
                    # apply the FDA at index i, j
                    v[i,j] = (v[i-1, j] + v[i+1, j] \
                            + v[i, j-1] + v[i, j+1]) / 4
            # find the difference between v vs. old_v for error bound
            diff = v - old_v
            # the highest value in diff is returned
            max_diff = np.max(np.abs(diff))
            # if the largest difference is less than 1e-10, break
            if (max_diff < 1e-10):
                break
        #print the np potential array and the loop runs required
        print(v)
        print(c)
        # flip the entire v np.array for proper rendering on contour graph
        v_plot = np.flipud(v)
        # find ∂V/∂y and ∂V/∂x numerically to calculate the electric field
        dV_dy, dV_dx = np.gradient(v_plot)
        # calculate the elctric field using E = -∇V
        E_x = -dV_dx
        E_y = -dV_dy
        # print the electric field values
        print(E_x)
        print(E_y)
        # plot the potential as a contour plot
        plt.contourf(v_plot)
        # plot the electric field as a quiver plot
        plt.quiver(E_x,E_y)
        # create a colourbar
        plt.colorbar()
        #show the graph
        plt.show()
    # 2D Poisson Solver: ∇^2V = -ρ(x, y)/ε₀ -> ∂^2V/∂x^2 + ∂^2V/∂y^2 = -ρ(x, y)/ε₀
    case 3:
        ## The boundary is a square, so the user must enter the electric potential on the left, reight, top and bottom
        left = float(input("What is the electric potential on the left boundary?\n"))
        right = float(input("What is the electric potential on the right boundary?\n"))
        bottom = float(input("What is the electric potential on the bottom boundary?\n"))
        top = float(input("What is the electric potential on the top boundary?\n"))
        # one must also enter the calculation resolution
        # resolution is the number of of grid points in a row (so dim^2 is the total grid points)
        dim = int(input("How many rows and columns in the calculation grid?\n"))
        # simultaineously we also need the actual size of the square domain (x, y)
        domain = float(input("What is the size of the square domain?\n"))
        # input string for ρ(x, y)
        rho_string = input("What is the charge density function?\n")
        # create ρ(x, y)
        rho = create_rho_2d(rho_string)
        # init np.array that represnts potential difference across the square
        v = np.zeros((dim, dim), dtype=float)
        # set each side of the np.array equal to there given values, i.e. the left side of the array has the value of left
        v[:, 0] = left
        v[:, -1] = right
        v[-1, :] = bottom
        v[0, :] = top
        # init the counter accumulator to counter how many times the Finite-Difference Approximation is used (through the entire array, i.e. updated when all elements have been)
        c = 0
        # h is the spatial grid spacing utilitized in the finite difference approximation
        # There are dim grid points, spanning a physical length domain, so there are dim - 1 intervals
        h = (domain)/ (dim-1)
        # the loop runs until it is broken
        while True:
            # accumulate onto counter
            c += 1
            # make old_v equal to a copy of v before the FDA is applied so that an error bound can be checked
            old_v = v.copy()
            # use the FDA across the entire square. NOTE: EXCLUDING the boundary values
            for i in range(1, v.shape[0] - 1):
                for j in range(1, v.shape[1] - 1):
                    # apply the FDA at index i, j
                    v[i,j] = (v[i-1, j] + v[i+1, j] \
                            + v[i, j-1] + v[i, j+1] + ((rho(i*h, j*h)/EPSILON_0) * h**2)) / 4
            # find the difference between v vs. old_v for error bound
            diff = v - old_v
            # the highest value in diff is returned
            max_diff = np.max(np.abs(diff))
            # if the largest difference is less than 1e-10, break
            if (max_diff < 1e-10):
                break
        #print the np potential array and the loop runs required
        print(v)
        print(c)
        # flip the entire v np.array for proper rendering on contour graph
        v_plot = np.flipud(v)
        # find ∂V/∂y and ∂V/∂x numerically to calculate the electric field
        dV_dy, dV_dx = np.gradient(v_plot, h, h)
        # calculate the elctric field using E = -∇V
        E_x = -dV_dx
        E_y = -dV_dy
        # print the electric field values
        print(E_x)
        print(E_y)

        # create array of even spaced numbers from 0 to domain with dim grid points
        x = np.linspace(0, domain, dim)
        y = np.linspace(0, domain, dim)

        # plot the potential as a contour plot
        contour = plt.contourf(x, y, v_plot)
        # plot the electric field as a quiver plot
        plt.quiver(x, y, E_x, E_y)
        # create a colourbar
        plt.colorbar(contour)
        # show the figure
        plt.show()
    # 3D Poisson Solver: ∇^2V = -ρ(x, y, z)/ε₀ -> ∂^2V/∂x^2 + ∂^2V/∂y^2 + ∂^2V/∂z^2 = -ρ(x, y, z)/ε₀
    # 3D Lapalace Solver is simply the special case where ρ(x, y, z) = 0
    case 4:
        ## The boundary is a cube, so the user must enter the electric potential on the left, right, top, bottom, above, and below
        left = float(input("What is the electric potential on the left boundary?\n"))
        right = float(input("What is the electric potential on the right boundary?\n"))
        bottom = float(input("What is the electric potential on the bottom boundary?\n"))
        top = float(input("What is the electric potential on the top boundary?\n"))
        down = float(input("What is the electric potential on the bottom of the cube boundary?\n"))
        up = float(input("What is the electric potential on the top of the cube boundary?\n"))
        # one must also enter the calculation resolution
        # resolution is the number of of grid points in a row (so dim^3 is the total grid points)
        dim = int(input("How many rows and columns in the calculation grid?\n"))
        # simultaineously we also need the actual size of the cubical domain (x, y, z)
        domain = float(input("What is the size of the cube domain?\n"))
        # input string for ρ(x, y, z)
        rho_string = input("What is the charge density function?\n")
        # create ρ(x, y, z)
        rho = create_rho_3d(rho_string)
        # init np.array that represnts potential difference across the cube
        v = np.zeros((dim, dim, dim), dtype=float)
        # set each face of the np.array cube equal to there given values, i.e. the left side of the array has the value of left
        v[0, :, :] = left
        v[-1, :, :] = right
        v[:, 0, :] = bottom
        v[:, -1, :] = top
        v[:, :, 0] = down
        v[:, :, -1] = up
        # init the counter accumulator to counter how many times the Finite-Difference Approximation is used (through the entire array, i.e. updated when all elements have been)
        c = 0
        # h is the spatial grid spacing utilitized in the finite difference approximation
        # There are dim grid points, spanning a physical length domain, so there are dim - 1 intervals
        h = (domain)/ (dim-1)
        # the loop runs until it is broken
        while True:
            # accumulate onto counter
            c += 1
            # make old_v equal to a copy of v before the FDA is applied so that an error bound can be checked
            old_v = v.copy()
            # use the FDA across the entire cube. NOTE: EXCLUDING the boundary values
            for i in range(1, v.shape[0] - 1):
                for j in range(1, v.shape[1] - 1):
                    for k in range(1, v.shape[2] - 1):
                        # apply the FDA at index i, j, k
                        v[i, j, k] = (v[i-1, j, k] + v[i+1, j, k] \
                                + v[i, j-1, k] + v[i, j+1, k] \
                                + v[i, j, k-1] + v[i, j, k+1] + ((rho(i*h, j*h, k*h)/EPSILON_0) * h**2)) / 6
            # find the difference between v vs. old_v for error bound
            diff = v - old_v
            # the highest value in diff is returned
            max_diff = np.max(np.abs(diff))
            # if the largest difference is less than 1e-10, break
            if (max_diff < 1e-10):
                break
        #print the np potential array and the loop runs required
        print(v)
        print(c)
        # find ∂V/∂x, ∂V/∂y, and ∂V/∂z numerically to calculate the electric field
        dV_dx, dV_dy, dV_dz = np.gradient(v, h, h, h)
        # calculate the elctric field using E = -∇V
        E_x = -dV_dx
        E_y = -dV_dy
        E_z = -dV_dz
        # print the electric field values
        print(E_x)
        print(E_y)
        print(E_z)

        # create array of even spaced numbers from 0 to domain with dim grid points
        x = np.linspace(0, domain, dim)
        y = np.linspace(0, domain, dim)
        # create the x, y coordinate matrix to plot the electric field and potential contour, indexed as a matrix
        X, Y= np.meshgrid(x, y, indexing = "ij")

        ## Since the potential and electric field inherently live in a 3D space,
        ## we have to display SLICES of the cubical surface where z is held at a constant value.
        ## We limit the amount of slices so the cube is actually able to be visualized
        # num_slices maakes it so that at max 6 slices are displayed
        num_slices = min(dim, 6)
        # choose the apprioate z values for the slices
        slices = np.linspace(0, dim - 1, num_slices, dtype=int)

        ## Potential Contour graphs
        # create the 2x3 figure, each fig is 12x12
        fig1, axes1 = plt.subplots(2, 3, figsize=(12,12))
        # loop through the 6 plots
        for i, ax in enumerate(axes1.flat):
            # if dim < 6, num_slices = 4, so to ensure the last 2 plots are empty, when i > the amount of slices
            if (i >= num_slices):
                # turn off that subplot and continue
                ax.axis("off")
                continue
            # k is the actula z-grid index, so slices[i] is the z value of the one plane we are graphing
            k = slices[i]
            # plot the potential of V(x, y, k)
            contour = ax.contourf(X, Y, v[:, :, k])
            # plot the x & y components of the Electric field through quiver on the contour grpah
            ax.quiver(X, Y, E_x[:, :, k], E_y[:, :, k], scale = 10)
            # create a colourbar for contour based on the value of the potential
            plt.colorbar(contour, ax=ax)
            # create a title for each subplot based one what z equals
            ax.set_title(f"z = {k * h}")
            # labels for x and y axis
            ax.set_xlabel("x")
            ax.set_ylabel("y")
        # making all 6 subplots have a tight layout
        fig1.tight_layout()

        ## E_z Contour graphs
        # create the 2x3 figure, each fig is 12x12
        fig2, axes2 = plt.subplots(2, 3, figsize=(12,12))
        # loop through the 6 plots
        for i, ax in enumerate(axes2.flat):
            # if dim < 6, num_slices = 4, so to ensure the last 2 plots are empty, when i > the amount of slices
            if (i >= num_slices):
                # turn off that subplot and continue
                ax.axis("off")
                continue
            # k is the actula z-grid index, so slices[i] is the z value of the one plane we are graphing
            k = slices[i]
            # plot the value of E_z at (x, y, k)
            contour = ax.contourf(X, Y, E_z[:, :, k])
            # plot the x & y components of the Electric field through quiver on the contour grpah
            ax.quiver(X, Y, E_x[:, :, k], E_y[:, :, k], scale = 10)
            # create a colourbar for contour based on the value of E_z
            plt.colorbar(contour, ax=ax)
            # create a title for each subplot based one what z equals
            ax.set_title(f"E_z, z = {k * h}")
            # labels for x and y axis
            ax.set_xlabel("x")
            ax.set_ylabel("y")
        # making all 6 subplots have a tight layout
        fig2.tight_layout()
        # Show the 2 figures
        plt.show()
    # if a valid output is not given
    case _:
        # end the program
        sys.exit()