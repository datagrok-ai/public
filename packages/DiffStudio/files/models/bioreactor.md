# Model
Controlled Fab-arm exchange assembly of a bispecific IgG1
from two parental antibodies (F405L and K409R), with a UF/DF
bioreactor wrapper ([source](https://doi.org/10.1074/jbc.RA117.000303)).
# Try
Interactive results based on input changes.
# Complexity
Each time you change inputs, this demo solves a system of 12 nonlinear ordinary differential equations.
# Quick reference

| Input | Meaning |
|---|---|
| Scenario | Load a saved input preset |
| Reduction | Length of 2-MEA reduction stage |
| Filtration | Length of UF/DF stage |
| Switch | UF → diafiltration handover |
| FFox | Parental F (F405L) antibody |
| KKox | Parental K (K409R) antibody |
| FFred | F homodimer, reduced hinges |
| KKred | K homodimer, reduced hinges |
| Ffree | Free F half-antibody |
| Kfree | Free K half-antibody |
| FKred | Bispecific, reduced hinges |
| FKox | Bispecific final product |
| MEAthiol | Cysteamine (2-MEA) dose |
| DO2 | Dissolved oxygen |
| yO2P | Headspace O₂ pressure |
| CYST | Cystamine (oxidized 2-MEA) |
| VL | Reactor liquid volume |
| Gas | Headspace gas flow rate |
| O₂ fraction | O₂ fraction in inlet gas |
| Temperature | Reactor temperature |
| Pressure | Headspace pressure |
| Initial | Simulation start time |
| Step | ODE solver step size |
