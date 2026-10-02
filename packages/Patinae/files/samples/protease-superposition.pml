# Superposition of two HIV-1 protease inhibitor complexes: indinavir (1HSG) and
# nelfinavir (1OHR). The second structure is aligned onto the first (Kabsch fit on CA atoms)
# and the two inhibitors are compared in the shared pocket.
load 1hsg.pdb, indinavir_complex
load 1ohr.pdb, nelfinavir_complex
align nelfinavir_complex, indinavir_complex, method=sequence
hide everything
show cartoon
set cartoon_transparency, 0.5
color gray80, indinavir_complex
color slate, nelfinavir_complex
select indinavir, indinavir_complex and resn MK1
select nelfinavir, nelfinavir_complex and resn 1UN
show sticks, indinavir or nelfinavir
color yellow, indinavir and elem C
color hotpink, nelfinavir and elem C
show sticks, (indinavir_complex and resi 25+50 and not name N+C+O)
bg_color white
orient indinavir or nelfinavir
scene inhibitors, store, message=Indinavir (yellow) and nelfinavir (pink) after superposition
orient
scene dimers, store, message=Both homodimers superposed
scene inhibitors
deselect
