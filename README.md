# Heat Diffusion Simulation
Heat diffusion simulation is a pet project programmed in Python to simulate the cooking process of a piece of tofu in a pan by using DOLFINx to run Finite-Element solving methods on basis of the heat diffusion equation. The goal is to study how many times I have to turn a tofu in a cast iron pan, that the center cooks fastest, while not burning the outside.

# Usage
The project includes a Variational_Problem class which can be modified by passing different boundary and starting conditions and different weak forms. It should be capable to solve multiple time-dependent linear PDEs, not only the heat equation, although proper support for this is still to be implemented, as well as for non-linear PDEs. The Program further features the option to insert special operations on the variational problem (like flipping the tofu) at any time step. Further a clean way to define different constants on the different sub-domains of the mesh is included. This offers a clean and fast way to answer the tofu question or similar problems you might encounter.

#Features
- The script mesh_creation.py creates a valid geometric mesh and assigns physical regions to it, then usable in the simulation
- Work is done in the jupyter notebook where we can define the functional, material constants and boundary conditions of the problem
- Insert additional operations at chooseable discrete time steps (change constants, fields, etc.)
- Definable start and finish time as well as changeable time steps
- Creates a .xdmf file that logs the heat-field over time and saves it for later analysis
- Right now the siulation is only capable of 3d or 2d geometries. The option of interfacing different dimensional geometries might be added later


# Installation
There is no installation procedure planned right now. However, if you want to try. The script requires:
- Python 3.11
- ufl
- gmsh
- gmshio
- DOLFINx
- NumPy
- mpi4py
- json
- pathlib

An working virtual environment will maybe implemented in the future

# Example
You can find an example on how to use this script in the notebook already where the simulation is run for the piece of tofu and the core temperature is analized.

# Acknowledgement
I want to thank Jørgen S. Dokken for the very helpful tutorial on https://jsdokken.com/dolfinx-tutorial/.
