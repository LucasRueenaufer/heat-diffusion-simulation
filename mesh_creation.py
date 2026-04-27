#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 09:41:11 2026

@author: LucasR

This script is providing meshes in a Dolfinx useable .xdmf file format, such
that different Finite-Element programs can be run on it.
"""

import gmsh
import glob
import meshio
from pathlib import Path

def main():
    temp_dir = Path("temp/meshing")
    temp_dir.mkdir(exist_ok=True,parents=True)
    skillet_mesher()
    mesh_to_dolfin()
    pass

def mesh_to_dolfin():
    mesh_list = glob.glob("temp/meshing/*.msh")
    
    for path in mesh_list:
        path = Path(path)
        
        mesh = meshio.read(path)
    
        # Extract tetrahedral cells
        tetra = mesh.get_cells_type("tetra")
        cell_data = mesh.get_cell_data("gmsh:physical", "tetra")
        
        mesh = meshio.Mesh(
            points=mesh.points,
            cells=[("tetra", tetra)],
            cell_data={"cell_tag": [cell_data]}
        )
    
        meshio.write(str(path.with_suffix(".xdmf")), mesh)
        
def skillet_mesher():
    # --- does the skillet need a handle? ---
    add_Handle = True
    
    # --- are we frying anything up? ---
    add_Tofu = True
    
    # --- initialize gmsh and name important variables ---
    gmsh.initialize()
    gmsh.model.add("skillet")

    model = gmsh.model
    occ = model.occ

    # --- Define Parameters of the skillet ---
    R_bottom_outer = 10     # bottom outer radius
    R_top_outer    = 12.5     # top outer radius 
    height         = 4      # inner height
    thickness      = 0.6    # wall thickness
    bottom_thickness = 0.8  # bottom thickness
    
    R_bottom_inner = R_bottom_outer - thickness
    R_top_inner    = R_top_outer - thickness
    
    total_height = height + bottom_thickness
    
    # --- define shape of skillet ---
    outer = occ.addCone(
        0, 0, 0,
        0, 0, total_height,
        R_bottom_outer,
        R_top_outer
        )
    
    # --- define geometry of cavisty ---
    inner = occ.addCone(
        0, 0, bottom_thickness,
        0, 0, height,
        R_bottom_inner,
        R_top_inner
        )
    
    # --- add optional Handle ---
    if add_Handle:
        # parameters of handle
        handle_length = 15
        handle_radius = 1
        
        # adding handle to skillet
        handle = occ.addCylinder(
            R_bottom_outer, 0 ,total_height*0.7, handle_length,0,0,handle_radius)
        
        skillet_body = occ.fuse([(3,outer)],[(3,handle)])[0]
    
    # --- subtract cavity from skillet ---
    skillet = occ.cut(skillet_body, [(3, inner)])[0]
    skillet_old = skillet
    # ---add optional tofu ---
    if add_Tofu:
        # parameters of tofu-block
        tofu_length = 8
        tofu_width = 6
        tofu_height = 3
        
        # adding tofu to skillet
        tofu = occ.addBox(-tofu_length/2,-tofu_width/2,bottom_thickness,tofu_length,tofu_width,tofu_height)
        
        # fragment domain into skillet and tofu, but keep tags
        frag_domains, frag_map = occ.fragment(skillet,[(3,tofu)])
    else:
        frag_domains, frag_map = occ.fragment(skillet)
    
    # --- synchonize ---
    occ.synchronize()
    
    # --- identify all material tags (volume) ---
    skillet_entities = [tag for (dim, tag) in frag_map[0]]
    
    if add_Tofu: 
        tofu_entities = [tag for (dim, tag) in frag_map[1]]
        
    # -- Add physical group to volume ---
    model.addPhysicalGroup(3, skillet_entities, tag=1)
    model.setPhysicalName(3, 1, "skillet")
    
    if add_Tofu:
        model.addPhysicalGroup(3, tofu_entities, tag=2)
        model.setPhysicalName(3, 2, "tofu")
    
    # --- mesh parameters ---
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.05)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.5)
    
    # --- generate mesh ---
    model.mesh.generate(3)
    
    # ---saving the result ---
    gmsh.write("temp/meshing/sloped_skillet.msh")
    
    gmsh.finalize()
    pass

if __name__ == "__main__":
    main()