#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 08:53:54 2026

@author: LucasR

This script will take a time dependent PDE problem in variational form and
solve by doing discretized time steps and treating every step as a 
steady state problem.
"""

import sys
from dolfinx.fem.petsc import (
    assemble_vector,
)

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

def run_simulation(t,T,DeltaT,uh,u_n,b,linear_form,solver,xdmf,A,plots):
    grid, plotter = plots
    num_steps=int(T/DeltaT)
    
    for i in range(num_steps):
        uh, u_n, t = time_step_solution([uh,u_n],solver,[b,linear_form],[t,T,DeltaT],grid)
        solution_writer(xdmf,grid,plotter,uh,t)
    
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