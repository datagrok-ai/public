# HIV-1 protease bound to indinavir (PDB 1HSG): representations, measurements,
# named scenes, and a keyframed camera movie that flies between them.
load 1hsg.pdb, protease
hide everything
show cartoon, protease
color gray80, protease and chain A
color lightblue, protease and chain B
select indinavir, protease and resn MK1
select pocket, byres (protease and polymer within 4.5 of indinavir)
show sticks, indinavir or (pocket and not name N+C+O)
color yellow, indinavir and elem C
color salmon, pocket and elem C
show spheres, protease and resi 308 and name O
set sphere_scale, 0.25

# Catalytic aspartates and the flap water that bridges the inhibitor to Ile50
distance catalytic, indinavir and name O2, protease and resi 25 and name OD1+OD2, cartesian
distance flapA, protease and resi 308 and name O, protease and chain A and resi 50 and name N
distance flapB, protease and resi 308 and name O, protease and chain B and resi 50 and name N
label pocket and name CA, resn
bg_color white

# Scenes: overview, binding pocket, pocket surface
orient protease
scene overview, store, message=HIV-1 protease homodimer with indinavir
orient indinavir
turn x, -25
scene pocket, store, message=Binding pocket: catalytic Asp25/Asp25' and the flap water
show surface, pocket
set surface_transparency, 0.55
set surface_color, wheat
scene surface, store, message=Pocket surface
scene overview
# Drop the named selections so their highlight overlay does not mask the atom colors
deselect

# Movie: 180 frames interpolated between the three scenes, then back to the overview
mset 1 x180
mview store, 1, scene=overview
mview store, 60, scene=pocket
mview store, 120, scene=surface
mview store, 180, scene=overview
mview interpolate
mplay
