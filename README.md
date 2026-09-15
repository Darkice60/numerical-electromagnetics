# numerical-electromagnetics
A collection of electromagnetism solvers and visualizers. These include 1D & 2D Laplace equation solvers, 2D & 3D Poisson equation solvers, and a Finite-Difference Time-Domain (FDTD) method solver for Maxwells equations.

# Features

The program has two main parts:

- The Laplace/Poisson Solver.
- The FDTD solver utilizing Yee's method.

The Laplace/Poisson solvers use finite-difference approximations to compute the potential and electric field within a defined domain. Its boundaries of the domain are an n-dimensional cube. 

The 2D Laplace/Poisson and the 3D Poisson calculator also include visualizations with matplotlib, showing the potential and electric field (in 3D there are multiple graphs showing potential, the direction of the electric field along the z-axis, and both have multiple layers at different z values).

The FDTD simulator utilizes Yee's method and its staggered grid. It is a numerical implementation of Faraday's Law of Induction and the Ampère–Maxwell Law.

# Dependencies
- Python 3
- Matplotlib
- Numpy