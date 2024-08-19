#!/bin/bash
set -e

config=${1}         # Config file name
ligand=${2}         # Ligand file name
# the autogrid paramters expected in the same file as receptor with pdbqt
config_basname="${config%.gpf}"
config_maps="${config_basname}.maps.fld"
ligand_basname="${ligand%.pdb}"

## Conda for automatic ligand and receptor creation
#mamba create -n autodock python=2.7
#mamba activate autodock
#mamba install -c insilichem autodocktools-prepare

# Activating conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate autodock

if [[ ! -f "${config_maps}" ]]; then
    autogrid4 -p "${config}" -l "${config_basname}.autogrid.log"
fi


prepare_ligand4.py -F -l "${ligand}"

# Autodock - GPU
adgpu-v1.5.3_linux_ocl_128wi --lfile "${ligand_basname}.pdbqt" --ffile "${config_maps}" --nrun 50 --resnam "${config_basname}.${ligand_basname}.autodock"





