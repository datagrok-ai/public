# Trp-cage miniprotein NMR ensemble (PDB 1L2Y, 38 models): every model is a state,
# the movie steps through the ensemble while the view rocks.
load 1l2y.pdb, trpcage
dss trpcage
hide everything
show cartoon, trpcage
color marine, trpcage and resi 1-4
color tv_blue, trpcage and resi 5-8
color palegreen, trpcage and resi 9-12
color yellow, trpcage and resi 13-16
color orange, trpcage and resi 17-20
show sticks, trpcage and resn TRP+PRO and not name N+C+O
color hotpink, trpcage and resn TRP and elem C
set cartoon_transparency, 0.15
bg_color white
orient trpcage
mset 1 -38
mplay
rock
