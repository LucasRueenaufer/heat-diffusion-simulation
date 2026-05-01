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
    inner,
    Measure,
    )

from dolfinx.fem import (
    locate_dofs_topological,
    dirichletbc,
    functionspace,
    Function,
    form
    )

from dolfinx.mesh import locate_entities, meshtags

import numpy as np
   
from ufl import lhs, rhs    
    
class VariationalProblem:
    def __init__(self,domain,Functional,Inital_Condition,T,const,DeltaT=0.1):
        # Saev lots of variables
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
        self.u_n = Function(self.FuncSpace)
        self.u_n.interpolate(Inital_Condition)
        
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
            facets = locate_entities(self.Domain,self.Fdim, locator)
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
    def build_forms(self,u):
        # update u_n with new value and build rhs and öhs forms
        self.u_n.x.array[:] = u.x.array
        a = lhs(self.Functional)
        L = rhs(self.Functional)
        return form(a), form(L)
    
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
            self._bc = dirichletbc(u_D, dofs)
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