#name: Generate Conformers
#description: Generates multiple conformers for a molecule using RDKit
#help-url: https://datagrok.ai/help/domains/chem/functions/conformers
#language: python
#meta.domain: chem
#top-menu: Chem | Calculate | Generate Conformers...
#input: string molecule = "CCCC" {semType: Molecule}  # Butane - shows anti/gauche conformations
#input: int num_conformers = 50 {caption: Num conformers} [Number of conformers to generate]
#input: bool optimize = true [Optimize with MMFF94 force field]
#input: double rms_threshold = 0.1 {caption: RMS threshold} [RMS threshold for conformer pruning - LOWER for butane]
#input: int max_attempts = 5000 {caption: Max attempts} [Maximum generation attempts]
#input: int random_seed = 42 {caption: Random seed} [Random seed for reproducibility]
#output: dataframe conformers

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom

def generate_conformers(mol, num_conformers=50, optimize=True,
                       rms_threshold=0.1, max_attempts=5000, random_seed=42):
    """
    Generate conformers for a molecule using RDKit's ETKDGv3 method.
    """

    # Add hydrogens for better geometry
    mol_with_h = Chem.AddHs(mol)

    # Set up conformer generation parameters
    params = rdDistGeom.ETKDGv3()
    params.randomSeed = random_seed
    params.maxAttempts = max_attempts
    params.pruneRmsThresh = rms_threshold  # CRITICAL: Lower threshold for butane
    params.useExpTorsionAnglePrefs = True
    params.useBasicKnowledge = True
    params.enforceChirality = False
    params.ignoreSmoothingFailures = True

    # Generate conformers
    conformer_ids = rdDistGeom.EmbedMultipleConfs(
        mol_with_h,
        numConfs=num_conformers,
        params=params
    )

    # If no conformers generated, try fallback
    if len(conformer_ids) == 0:
        AllChem.EmbedMolecule(mol_with_h, randomSeed=random_seed)
        conformer_ids = [0]

    # Optimize conformers if requested
    energies = []
    if optimize and len(conformer_ids) > 0:
        try:
            # Optimize with MMFF94
            results = AllChem.MMFFOptimizeMoleculeConfs(mol_with_h, maxIters=1000)

            # Calculate energies
            mp = AllChem.MMFFGetMoleculeProperties(mol_with_h, mmffVariant='MMFF94')
            for conf_id in conformer_ids:
                if mp:
                    ff = AllChem.MMFFGetMoleculeForceField(mol_with_h, mp, confId=conf_id)
                    energies.append(ff.CalcEnergy() if ff else np.nan)
                else:
                    energies.append(np.nan)
        except:
            # If optimization fails, just use the raw conformers
            energies = [np.nan] * len(conformer_ids)
    else:
        # Set energies to NaN if not optimized
        energies = [np.nan] * len(conformer_ids)

    # Create separate molecules for each conformer
    conformer_mols = []
    for conf_id in conformer_ids:
        # Create a copy of the molecule with just this conformer
        conf_mol = Chem.Mol(mol_with_h)
        conf_mol.RemoveAllConformers()
        conf_mol.AddConformer(mol_with_h.GetConformer(conf_id))
        conformer_mols.append(conf_mol)

    return conformer_mols, energies

# Parse input molecule
if "M  END" in molecule:
    mol = Chem.MolFromMolBlock(molecule)
else:
    mol = Chem.MolFromSmiles(molecule)

if mol is None:
    # Return empty dataframe if molecule is invalid
    conformers = pd.DataFrame(columns=['smiles', 'molblock', 'conformer', 'energy', 'rmsd'])
else:
    # Generate conformers
    conformer_mols, energies = generate_conformers(
        mol,
        num_conformers=num_conformers,
        optimize=optimize,
        rms_threshold=rms_threshold,
        max_attempts=max_attempts,
        random_seed=random_seed
    )

    # Calculate RMSD between conformers for diversity analysis
    rmsd_values = []
    if len(conformer_mols) > 1:
        # Use first conformer as reference
        ref_mol = conformer_mols[0]
        for i, conf_mol in enumerate(conformer_mols):
            if i == 0:
                rmsd_values.append(0.0)
            else:
                try:
                    rmsd = AllChem.GetBestRMS(ref_mol, conf_mol)
                    rmsd_values.append(rmsd)
                except:
                    rmsd_values.append(np.nan)
    else:
        rmsd_values = [0.0] * len(conformer_mols) if conformer_mols else []

    # Prepare data for DataFrame
    data = []
    for i, (conf_mol, energy, rmsd) in enumerate(zip(conformer_mols, energies, rmsd_values)):
        data.append({
            'smiles': Chem.MolToSmiles(mol),  # Original molecule SMILES
            'molblock': Chem.MolToMolBlock(conf_mol),
            'conformer': i + 1,
            'energy': energy,
            'rmsd': rmsd
        })

    # Create DataFrame
    conformers = pd.DataFrame(data)

# Ensure the output has the correct semType for molecules
conformers
