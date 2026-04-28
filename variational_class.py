#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 09:13:23 2026

@author: LucasR

This script will take a time dependent PDE problem in variational form and
solve by doing discretized time steps and treating every step as a 
steady state problem.
"""

from ufl import (
    TestFunction,
    TrialFunction,
    dx,
    functionspace,
    inner,
    grad,
    div,
    )

from dolfinx.fem import (
    locate_dofs_topological,
    dirichletbc,
    )

from dolfinx.io import XDMFFile

from mpi4py import MPI

with XDMFFile(MPI.COMM_WORLD, "temp/meshing/sloped_skillet.xdmf", "r") as xdmf:
    domain = xdmf.read_mesh(name="Grid")
    cell_tags = xdmf.read_meshtags(domain, name="Grid")
    
    
    
class Variational_Problem:
    def __init__(self,domain,Functional,T):
        self.Domain = domain,
        self.FuncSpace = functionspace(domain, ("Lagrange", 1))
        self.TrialFunc = TrialFunction(self.FuncSpace)
        self.TestFunc = TestFunction(self.FuncSpace)
        self.Measure = dx
        self.Functional = Functional(self.TrialFunc,self.TestFunc,self.Measure)
        self.bcs = [None]
        self.fdim = domain.topology.dim - 1
        self.DeltaT = 0.1
    def AddBC(self,boundary_conditions):
        for condition in boundary_conditions:
            if condition.type == "Dirichlet":
                self.bcs.append(condition.bc)
            else:
                self.Functional += condition.bc
        
    
class BoundaryCondition:
    def __init__(self, type, marker, values, VP):
        self._type = type
        if type == "Dirichlet":
            u_D = Variational_Problem.FuncSpace
            u_D.interpolate(values)
            facets = facet_tag.find(marker)
            dofs = locate_dofs_topological(VP.FuncSpace, VP.fdim, facets)
            self._bc = dirichletbc(u_D, dofs)
        elif type == "Neumann":
            self._bc = inner(values, VP.TestFunc) * ds(marker)
        elif type == "Robin":
            self._bc = values[0] * inner(VP.TrialFunc - values[1], VP.TestFunc) * ds(marker)
        else:
            raise TypeError("Unknown boundary condition: {0:s}".format(type))

    @property
    def bc(self):
        return self._bc

    @property
    def type(self):
        return self._type

T=5

def Functional(u,v,Dx):
    return inner(grad(u),grad(v))*Dx

VP = Variational_Problem(domain, Functional, T)

# Define the Dirichlet condition
boundary_conditions = [
    BoundaryCondition("Dirichlet", 1,0,VP),
]