# Patinae demo script: loads the crambin structure next to this file and styles it
load 1crn.pdb, crambin
show cartoon, crambin
color marine, \
  crambin and ss h
color yellow, crambin and ss s
show sticks, crambin and resn CYS
zoom crambin
