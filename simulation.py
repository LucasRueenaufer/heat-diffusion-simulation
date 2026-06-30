#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 09:13:23 2026

@author: LucasR

This script will take a time dependent PDE problem in variational form and
solve by doing discretized time steps and treating every step as a 
steady state problem.
"""
from tqdm.notebook import tqdm
from ufl import (
    TestFunction,
    TrialFunction,
    dx,
    inner,
    Measure,
    )
from dolfinx.fem.petsc import (
    assemble_vector,
    assemble_matrix,
    create_vector,
    apply_lifting,
    set_bc,
)

from dolfinx.fem import (
    locate_dofs_topological,
    dirichletbc,
    functionspace,
    extract_function_spaces,
    Function,
    form
    )

from dolfinx.mesh import locate_entities_boundary, meshtags

from dolfinx.io import XDMFFile

import numpy as np
from petsc4py import PETSc
import matplotlib.pyplot as plt

from dolfinx import geometry

import additional_operations as addop
   
from ufl import lhs, rhs    

class VariationalProblem:
    def __init__(self,domain,region,Functional,Inital_Condition,T,t,const,temp_dir,DeltaT=0.1):
        # Saev lots of variables
        self.temp_dir=temp_dir
        self.t = t
        self.DeltaT = DeltaT
        self.T= T
        self.Domain = domain
        self.FuncSpace = functionspace(self.Domain, ("Lagrange", 1))
        self.TrialFunc = TrialFunction(self.FuncSpace)
        self.TestFunc = TestFunction(self.FuncSpace)
        self.Measure = dx
        self.Const = const
        self.Bcs = []
        self.Fdim = domain.topology.dim - 1
        self.Tree = geometry.bb_tree(domain,self.Fdim+1)
        self.u_n = Function(self.FuncSpace)
        self.u_n.interpolate(Inital_Condition)
        self.uh =self.u_n.copy()
        self.uh.name = "uh"
        self.Region = region
        self.AvgT = []
        # declare functional
        self.Functional = Functional(self.TestFunc,self.TrialFunc,self.Measure,self.Const,self.DeltaT,self.u_n)
        print("Boundary Problem Initiated")
        
    def AddBC(self,boundary_conditions):
        # add boundary conditions to functional for wobin and neumann or
        # to list for dirichlet
        for condition in boundary_conditions:
            if condition.type == "Dirichlet":
                self.Bcs.append(condition.bc)
            else:
                self.Functional += condition.bc
        print("Boundary Conditions added to Variational Problem!")
                
    def TagFacet(self,boundaries):
        # tag all boundary facets that fulfill external locator
        facet_indices, facet_markers = [], []
        for marker, locator in boundaries:
            facets = locate_entities_boundary(self.Domain,self.Fdim, locator)
            facet_indices.append(facets)
            facet_markers.append(np.full_like(facets, marker))
        facet_indices = np.hstack(facet_indices).astype(np.int32)
        facet_markers = np.hstack(facet_markers).astype(np.int32)
        sorted_facets = np.argsort(facet_indices)
        
        # save in facet_tag format 
        self.Facet_tag = meshtags(
            self.Domain, self.Fdim, facet_indices[sorted_facets], facet_markers[sorted_facets]
        )
        
        # debug topology
        self.Domain.topology.create_connectivity(self.Fdim, self.Fdim+1)
        print("Tagged Boundary Facets!")
        
    def prepare_solver(self):
        # Build forms and assemble matrix and vector
        self.bilinear_form = form(lhs(self.Functional))
        self.linear_form = form(rhs(self.Functional))
        
        self.A = assemble_matrix(self.bilinear_form,bcs=self.Bcs)
        self.A.assemble()
        self.b = create_vector(extract_function_spaces(self.linear_form))
        
        # Initiate Solver
        self.solver = PETSc.KSP().create(self.Domain.comm)
        self.solver.setOperators(self.A)
        self.solver.setType(PETSc.KSP.Type.PREONLY)
        self.solver.getPC().setType(PETSc.PC.Type.LU)
        print("Solver Initiated!")
    
    def time_step(self):
        #step in time
        self.t += self.DeltaT
        
        # update vector by using initial form
        with self.b.localForm() as loc_b:
            loc_b.set(0)
        assemble_vector(self.b,self.linear_form)
        
        #prepare vector for dirichlet boundary conditions
        apply_lifting(self.b,[self.bilinear_form],[self.Bcs])
        self.b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        set_bc(self.b,self.Bcs)
        
        #solve linear problem
        self.solver.solve(self.b,self.uh.x.petsc_vec)
        self.uh.x.scatter_forward()
        
        # update last state
        self.u_n.x.array[:] = self.uh.x.array
    
    def destroy_solver(self):
        # print some info and destroy instances to prevent memory leaks
        print("Cleaning up...")
        self.A.destroy()
        self.b.destroy()
        self.solver.destroy()
        
    def prepare_writer(self):
        self.xdmf = XDMFFile(self.Domain.comm, self.temp_dir/"diffusion.xdmf", "w")
        self.xdmf.write_mesh(self.Domain)
        self.xdmf.write_function(self.uh, 0)

    def run_simulation(self):
        # initialize solver
        self.prepare_solver()
        self.prepare_writer()
        
        # use tqdm as iterator for nice progress bar
        for i in tqdm(range(int(self.T/self.DeltaT)),desc="Solving PDE"):
            # run single time step
            self.time_step()
            # write the solution
            self.xdmf.write_function(self.uh,self.t)
            
            # run any additional operation of the data
            addop.additional_operations(self)
        
        #plot the average core temp over time
        data=np.array(self.AvgT)
        plt.scatter(data[:, 0], data[:, 1])
        plt.xlabel("Time")
        plt.ylabel("Temperature")
        plt.grid(True)
        plt.show()
        
        # do not let any data leaks happen
        self.destroy_solver()
        self.xdmf.close()
        print("Simulation Finished!")
    
    
class BoundaryCondition:
    def __init__(self, type, marker, values, VP):
        self._type = type
        
        # define surface measure
        ds = Measure("ds", domain = VP.Domain, subdomain_data=VP.Facet_tag)
        
        # sort bcs into categories
        if type == "Dirichlet":
            u_D = VariationalProblem.FuncSpace
            u_D.interpolate(values)
            facets = VP.facet_tag.find(marker)
            dofs = locate_dofs_topological(VP.FuncSpace, VP.fdim, facets)
            self._bc = dirichletbc(u_D, dofs,VP.FuncSpace)
        elif type == "Neumann":
            self._bc = inner(values, VP.TestFunc) * ds(marker)
        elif type == "Robin":
            self._bc = values[0] * inner(VP.TrialFunc - values[1], VP.TestFunc) * ds(marker)
        else:
            raise TypeError("Unknown boundary condition: {0:s}".format(type))
    
    @property
    def bc(self):
    # return the bc
        return self._bc

    @property
    def type(self):
    # return the type of bc
        return self._type
    
def main():
    pass

if __name__ == "__main__":
    main()