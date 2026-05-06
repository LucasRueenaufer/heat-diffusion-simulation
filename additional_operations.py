#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 11:09:55 2026

@author: lucas
"""

import numpy as np
from dolfinx.fem import Function
from dolfinx import fem, geometry
import simulation as VC

def additional_operations(VP):
    """
    Is called during every simulation step. Write any additional
    operations that need to be fulfilled at the simulated timesteps
    here.
    """
    fliptimes=np.array([15,30,45])/VP.DeltaT
    fliptimes=fliptimes.astype(np.int32)
    
    # if we reach a flipping step, we flip the tofu
    if int(VP.t/VP.DeltaT) in fliptimes:
        VP.u_n.x.array[:]=flip_tofu_z(VP)
        print(f"Flipping Tofu at t={VP.t}")
        
def flip_tofu_z(VP):
    """
    Flip the tofu temperature in z by interpolation.

    Assumes:
      - rectangular tofu
      - same mesh before/after flip
      - pan remains unchanged
    """
    z_min = 0.8
    z_max = 0.8+2.5
    u_new = Function(VP.FuncSpace)
    
    u_new.x.array[:] = VP.uh.x.array

    tdim = VP.Fdim+1
    tree = geometry.bb_tree(VP.Domain, tdim)

    tofu_cells = VP.Region["tofu"]

    def flipped(x):
        xr = x.copy()
        xr[2] = z_min + z_max - xr[2]

        # find candidate cells
        xq = xr.T
        candidates = geometry.compute_collisions_points(tree, xq)
        colliding = geometry.compute_colliding_cells(VP.Domain, candidates, xq)

        cells = np.array(
            [colliding.links(i)[0] for i in range(len(xq))],
            dtype=np.int32,
        )

        # IMPORTANT: eval now has matching (x, cells)
        return VP.uh.eval(xq, cells).T

    u_new.interpolate(flipped, tofu_cells)
    u_new.x.scatter_forward()

    return u_new.x.array

# validate region map by comparing material map to all regions
def region_validation(region, region_map, material_const):
    # extract regions & materials
    region_cells = set(region.keys())
    region_materials = set(region_map.keys())

    # substract sets
    missing_region = region_cells - region_materials
    unknown_region = region_materials - region_cells

    # check if material actually exists in database
    unknown_materials = set()
    for reg, mat in region_map.items():
        if mat not in material_const:
            unknown_materials.add(reg)

    # check if all regions are valid and asigned a known material
    if unknown_materials:
        raise SystemExit(
            f"Undefined materials used in Regions: {sorted(unknown_materials)}"
    )

    if missing_region:
        raise SystemExit(
            f"No Materials assigned to regions: {list(missing_region)}\nAssignable materials are {list(region_materials)}"
        )

    if unknown_region:
        raise SystemExit(
            f"Materials assigned to unknown regions: {list(unknown_region)}\nAssignable regions are {list(region_cells)}"
        )

    if not missing_region and not unknown_region:
        print("All regions correctly mapped to materials")
    pass