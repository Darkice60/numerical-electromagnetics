# Darkice60
# main file for numerical-electromagnetics

import numpy as np
import math
import matplotlib.pyplot as plt
from typing import Final

choice = 3

EPSILON_0: Final = 8.8541878128e-12

def create_rho(rho_string):
    def rho(x, y):
        return eval(rho_string, {"x": x, "y": y, "math": math, "np": np})
    return rho

match choice:
    case 1:
        initial = float(input("What is the initial electric potential"))
        final = float(input("What is the final electric potential"))
        v = np.array([initial, 0, 0, 0, 0, final], dtype=float)
        c = 0
        while True:
            c += 1
            old_v = v.copy()
            for i in range(1, len(v) -1):
                v[i] =(v[i-1] + v[i+1]) /2
            diff = v - old_v
            max_diff = np.max(np.abs(diff))
            if (max_diff < 1e-10):
                break
        print(v)
        print(c)
    case 2:
        left = float(input("What is the minimum x electric potential"))
        right = float(input("What is the maximum x electric potential"))
        bottom = float(input("What is the minimum y in electric potential"))
        top = float(input("What is the maximum y electric potential"))
        dim = int(input("How many rows and columns in the calculation grid?"))
        v = np.zeros((dim, dim), dtype=float)
        v[:, 0] = left
        v[:, -1] = right
        v[-1, :] = bottom
        v[0, :] = top
        c = 0
        while True:
            old_v = v.copy()
            c += 1
            for i in range(1, v.shape[0] - 1):
                for j in range(1, v.shape[1] - 1):
                    v[i,j] = (v[i-1, j] + v[i+1, j] \
                            + v[i, j-1] + v[i, j+1]) / 4
            diff = v - old_v
            max_diff = np.max(np.abs(diff))
            if (max_diff < 1e-10):
                break
        print(v)
        print(c)
        v_plot = np.flipud(v)
        dV_dy, dV_dx = np.gradient(v_plot)
        E_x = -dV_dx
        E_y = -dV_dy
        print(E_x)
        print(E_y)
        plt.contourf(v_plot)
        plt.quiver(E_x,E_y)
        plt.colorbar()
        plt.show()
    case 3:
        left = float(input("What is the electric potential on the left boundary?"))
        right = float(input("What is the electric potential on the right boundary?"))
        bottom = float(input("What is the electric potential on the bottom boundary?"))
        top = float(input("What is the electric potential on the top boundary?"))
        dim = int(input("How many rows and columns in the calculation grid?"))
        domain = float(input("What is the size of the square domain?"))
        rho_string = input("What is the charge density function?")
        rho = create_rho(rho_string)
        v = np.zeros((dim, dim), dtype=float)
        v[:, 0] = left
        v[:, -1] = right
        v[-1, :] = bottom
        v[0, :] = top
        c = 0
        h = (domain)/ (dim-1)
        while True:
            old_v = v.copy()
            c += 1
            for i in range(1, v.shape[0] - 1):
                for j in range(1, v.shape[1] - 1):
                    v[i,j] = (v[i-1, j] + v[i+1, j] \
                            + v[i, j-1] + v[i, j+1] + ((rho(i*h, j*h)/EPSILON_0) * h**2)) / 4
            diff = v - old_v
            max_diff = np.max(np.abs(diff))
            if (max_diff < 1e-10):
                break
        print(v)
        print(c)
        v_plot = np.flipud(v)
        dV_dy, dV_dx = np.gradient(v_plot, h, h)
        E_x = -dV_dx
        E_y = -dV_dy
        print(E_x)
        print(E_y)

        x = np.linspace(0, domain, dim)
        y = np.linspace(0, domain, dim)

        contour = plt.contourf(x, y, v_plot)
        plt.quiver(x, y, E_x, E_y)
        plt.colorbar(contour)
        plt.show()
    case _:
        sys.exit()