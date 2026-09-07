# Darkice60
# main file for numerical-electromagnetics

import numpy as np
import math
import matplotlib.pyplot as plt
from typing import Final

choice = 4

EPSILON_0: Final = 8.8541878128e-12

def create_rho_2d(rho_string):
    def rho(x, y):
        return eval(rho_string, {"x": x, "y": y, "math": math, "np": np})
    return rho

def create_rho_3d(rho_string):
    def rho(x, y, z):
        return eval(rho_string, {"x": x, "y": y, "z":z, "math": math, "np": np})
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
            for k in range(1, len(v) -1):
                v[k] =(v[k-1] + v[k+1]) /2
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
            for k in range(1, v.shape[0] - 1):
                for j in range(1, v.shape[1] - 1):
                    v[k,j] = (v[k-1, j] + v[k+1, j] \
                            + v[k, j-1] + v[k, j+1]) / 4
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
        rho = create_rho_2d(rho_string)
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
            for k in range(1, v.shape[0] - 1):
                for j in range(1, v.shape[1] - 1):
                    v[k,j] = (v[k-1, j] + v[k+1, j] \
                            + v[k, j-1] + v[k, j+1] + ((rho(k*h, j*h)/EPSILON_0) * h**2)) / 4
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
    case 4:
        # left = float(input("What is the electric potential on the left boundary?"))
        # right = float(input("What is the electric potential on the right boundary?"))
        # bottom = float(input("What is the electric potential on the bottom boundary?"))
        # top = float(input("What is the electric potential on the top boundary?"))
        # down = float(input("What is the electric potential on the bottom of the cube boundary?"))
        # up = float(input("What is the electric potential on the top of the cube boundary?"))
        dim = int(input("How many rows and columns in the calculation grid?"))
        left, bottom, down = 0, 0, 0
        right, top, up = 10, 10, 10
        domain = float(input("What is the size of the cube domain?"))
        # rho_string = input("What is the charge density function?")
        rho_string = "0"
        rho = create_rho_3d(rho_string)
        v = np.zeros((dim, dim, dim), dtype=float)
        v[0, :, :] = left
        v[-1, :, :] = right
        v[:, 0, :] = bottom
        v[:, -1, :] = top
        v[:, :, 0] = down
        v[:, :, -1] = up
        c = 0
        h = (domain)/ (dim-1)
        while True:
            old_v = v.copy()
            c += 1
            for i in range(1, v.shape[0] - 1):
                for j in range(1, v.shape[1] - 1):
                    for k in range(1, v.shape[2] - 1):
                        v[i, j, k] = (v[i-1, j, k] + v[i+1, j, k] \
                                + v[i, j-1, k] + v[i, j+1, k] \
                                + v[i, j, k-1] + v[i, j, k+1] + ((rho(i*h, j*h, k*h)/EPSILON_0) * h**2)) / 6
            diff = v - old_v
            max_diff = np.max(np.abs(diff))
            if (max_diff < 1e-10):
                break
        print(v)
        print(c)
        dV_dx, dV_dy, dV_dz = np.gradient(v, h, h, h)
        E_x = -dV_dx
        E_y = -dV_dy
        E_z = -dV_dz
        # print(E_x)
        # print(E_y)
        # print(E_z)

        x = np.linspace(0, domain, dim)
        y = np.linspace(0, domain, dim)
        X, Y= np.meshgrid(x, y, indexing = "ij")

        num_slices = min(dim, 6)
        slices = np.linspace(0, dim - 1, num_slices, dtype=int)

        fig1, axes1 = plt.subplots(2, 3, figsize=(12,12))
        for i, ax in enumerate(axes1.flat):
            if (i >= num_slices):
                ax.axis("off")
                continue
            k = slices[i]
            contour = ax.contourf(X, Y, v[:, :, k])
            ax.quiver(X, Y, E_x[:, :, k], E_y[:, :, k], scale = 10)
            plt.colorbar(contour, ax=ax)
            ax.set_title(f"z = {k * h}")
            ax.set_xlabel("x")
            ax.set_ylabel("y")
        fig1.tight_layout()

        fig2, axes2 = plt.subplots(2, 3, figsize=(12,12))
        for i, ax in enumerate(axes2.flat):
            if (i >= num_slices):
                ax.axis("off")
                continue
            k = slices[i]
            contour = ax.contourf(X, Y, E_z[:, :, k])
            ax.quiver(X, Y, E_x[:, :, k], E_y[:, :, k], scale = 10)
            plt.colorbar(contour, ax=ax)
            ax.set_title(f"E_z, z = {k * h}")
            ax.set_xlabel("x")
            ax.set_ylabel("y")
        fig2.tight_layout()
        plt.show()
    case _:
        sys.exit()