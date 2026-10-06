# Crambin (PDB 1CRN) with a density map contoured as a mesh, plus crystallographic
# symmetry mates generated from the unit cell.
load 1crn.pdb, crambin
load crambin-density.ccp4, density
hide everything
show sticks, crambin
color gray70, crambin and elem C
isomesh mesh, density, 1.0
set mesh_color, palecyan
bg_color white
orient crambin
scene asu, store, message=Asymmetric unit with the density mesh
symexp sym, crambin, crambin, 6.0
show cartoon, sym*
color palegreen, sym*
hide sticks, sym*
orient
scene mates, store, message=Symmetry mates within 6 A of the asymmetric unit
scene asu
