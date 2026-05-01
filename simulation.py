#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 08:53:54 2026

@author: LucasR

This script will take a time dependent PDE problem in variational form and
solve by doing discretized time steps and treating every step as a 
steady state problem.
"""

from dolfinx.fem.petsc import (
    assemble_vector,
    assemble_matrix,
    create_vector,
    apply_lifting,
    set_bc,
    
)
import sys
from petsc4py import PETSc

from dolfinx.fem import form, extract_function_spaces, Function

from ufl import lhs, rhs

import variational_class

def main():
    pass
  
def progress_bar(percent):
    bar_length = 40
    filled = int(bar_length * percent/100)
    bar = "█" * filled + "-" * (bar_length - filled)

    sys.stdout.write(f"\r|{bar}| {percent:6.2f}%")
    sys.stdout.flush()
    pass

def solution_writer(xdmf,grid,plotter,uh,t):
    # Write solution to file
    xdmf.write_function(uh, t)
    # Update plot
    grid.point_data["uh"][:] = uh.x.array
    plotter.render()
    plotter.write_frame()
    pass
    
def time_step_solution(func, solver, rhs,time_param,grid):
    t, T, DeltaT = time_param
    
    b, linear_form = rhs
    
    uh, u_n = func
    #step up time and print loading bar
    t += DeltaT
    percent = 100 * t / T
    
    progress_bar(percent)
    
    # Update the right hand side reusing the initial vector
    with b.localForm() as loc_b:
        loc_b.set(0)
        
    assemble_vector(b, linear_form)
    
    # Solve linear problem
    solver.solve(b, uh.x.petsc_vec)
    uh.x.scatter_forward()
        
    # Update solution at previous time step (u_n)
    u_n.x.array[:] = uh.x.array
    return uh, u_n, t

def run_simulation(VP,plots,xdmf):
    t = 0
    
    # unpack values
    grid, plotter = plots
    num_steps=int(VP.T/VP.DeltaT)
    
    # interpolate starting conditions
    u_n = Function(VP.FuncSpace)
    uh = Function(VP.FuncSpace)
    uh.name ="uh"
    u_n.name ="u_n"
    
    uh.interpolate(VP.u_n)
    u_n.interpolate(VP.u_n)
    
    # set up variational forms
    a = lhs(VP.Functional)
    L = rhs(VP.Functional)

    bilinear_form = form(a)
    linear_form = form(L)

    A = assemble_matrix(bilinear_form, bcs=VP.Bcs)
    A.assemble()
    b= create_vector(extract_function_spaces(linear_form))
    
    
    # Solve linear variational problem
    solver = PETSc.KSP().create(VP.Domain.comm)
    solver.setOperators(A)
    solver.setType(PETSc.KSP.Type.PREONLY)
    solver.getPC().setType(PETSc.PC.Type.LU)

    # run for every time step
    for i in range(num_steps):
        # rebuild lhs and rhs forms with new u_n
        bilinear_form, linear_form = VP.build_forms(u_n)
        
        # run single time step
        uh, u_n, t = time_step_solution([uh,u_n],solver,[b,linear_form],[t,VP.T,VP.DeltaT],grid)
        
        # save solutions in plot and xdmf
        solution_writer(xdmf,grid,plotter,uh,t)
    
    # print some info and destroy instances to prevent memory leaks
    print("\nSimulation Completed.\n\nWriting Files...")
    plotter.close()
    xdmf.close()
    print("Done!\n\nCleaning up...")
    A.destroy()
    b.destroy()
    solver.destroy()
    print("Done!")
    pass

if __name__ == "__main__":
    main()