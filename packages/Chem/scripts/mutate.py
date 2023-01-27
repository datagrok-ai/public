#name: Mutate
#description: Mutate molecule
#help-url: https://datagrok.ai/help/domains/chem/functions/mutate
#language: python
#tags: demo, chem, rdkit
#top-menu: Chem | Mutate...
#input: string molecule = "CN1C(CC(O)C1=O)C1=CN=CC=C1" {semType: Molecule}
#input: int steps = 1 [Number of mutation steps]
#input: bool randomize = true [Randomize mutations]
#input: int maxRandomResults = 100 [Maximum random results to calculate]
#output: dataframe mutations {semType: Molecule}

import pickle
import random
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import EditableMol

# C, N, O, S, Cl
atomic_numbers = [6, 7, 8, 16, 17]


def add_atom(mol, atomic_number, atom):
    emol = EditableMol(mol)
    new_index = emol.AddAtom(Chem.Atom(atomic_number))
    emol.AddBond(atom.GetIdx(), new_index, Chem.BondType.SINGLE)
    mol_ = emol.GetMol()
    Chem.SanitizeMol(mol_)
    return mol_


def add_atoms(mol):
    mols = []
    for atomic_number in atomic_numbers:
        for atom in mol.GetAtoms():
            if atom.GetImplicitValence() == 0:
                continue
            mols.append(add_atom(mol, atomic_number, atom))
    return mols


def sanitize_emol(emol):
    try:
        m2 = emol.GetMol()
        Chem.SanitizeMol(m2)
        return m2
    except:
        return None


def add_bond_helper(mol, a1, a2, bond_type):
    emol = EditableMol(mol)
    emol.AddBond(a1.GetIdx(), a2.GetIdx(), bond_type)
    return sanitize_emol(emol)


def mutate_bond(mol, bond, bond_type):
    my_mol = deep_copy_mol(mol)
    my_bond = my_mol.GetBondWithIdx(bond.GetIdx())
    my_bond.SetBondType(bond_type)
    Chem.SanitizeMol(my_mol)
    return my_mol


def add_bond_for_atoms(mol, a1, a2):
    v1, v2 = a1.GetImplicitValence(), a2.GetImplicitValence()
    bond = mol.GetBondBetweenAtoms(a1.GetIdx(), a2.GetIdx())
    mols = []
    if bond is None:
        mols.append(add_bond_helper(mol, a1, a2, Chem.BondType.SINGLE))
        if v1 > 1 and v2 > 1:
            mols.append(add_bond_helper(mol, a1, a2, Chem.BondType.DOUBLE))
        if v1 > 2 and v2 > 2:
            mols.append(add_bond_helper(mol, a1, a2, Chem.BondType.TRIPLE))
        return mols

    if bond.GetBondType() == Chem.BondType.SINGLE:
        if v1 > 1 and v2 > 1:
            mols.append(mutate_bond(mol, bond, Chem.BondType.DOUBLE))
        if v1 > 2 and v2 > 2:
            mols.append(mutate_bond(mol, bond, Chem.BondType.TRIPLE))
        return mols

    if bond.GetBondType() == Chem.BondType.DOUBLE:
        if v1 > 2 and v2 > 2:
            mols.append(mutate_bond(mol, bond, Chem.BondType.TRIPLE))
        return mols
    return mols


def add_bonds(mol):
    mols = []
    num_atoms = mol.GetNumAtoms()
    for i1 in range(num_atoms):
        a1 = mol.GetAtomWithIdx(i1)
        if a1.GetImplicitValence() == 0:
            continue
        for i2 in range(i1 + 1, num_atoms):
            a2 = mol.GetAtomWithIdx(i2)
            if a2.GetImplicitValence() == 0:
                continue
            mols.extend(add_bond_for_atoms(mol, a1, a2))
    return mols


def get_sub_mols(mol):
    if mol is None:
        return []
    frags = Chem.rdmolops.GetMolFrags(mol, asMols=True)
    if len(frags) == 1:
        return frags
    is_one_atom = False
    for frag in frags:
        if frag.GetNumAtoms() == 1:
            is_one_atom = True
    if not is_one_atom:
        return []
    for frag in frags:
        if frag.GetNumAtoms() != 1:
            return [frag]


def create_emol(mol):
    return EditableMol(mol)


def deep_copy_mol(mol):
    pkl = pickle.dumps(mol)
    return pickle.loads(pkl)


def remove_bond(mol, bond):
    mols = []
    emol = create_emol(mol)
    emol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    new_mol = sanitize_emol(emol)
    mols.extend(get_sub_mols(new_mol))
    if bond.GetBondType() in [Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]:
        mols.append(mutate_bond(mol, bond, Chem.BondType.SINGLE))
    if bond.GetBondType() == Chem.BondType.TRIPLE:
        mols.append(mutate_bond(mol, bond, Chem.BondType.DOUBLE))
    return mols


def remove_bonds(mol):
    mols = []
    for bond in mol.GetBonds():
        mols.extend(remove_bond(mol, bond))
    return mols


def get_all_mutations(molecule, steps=1):
    mol = Chem.MolFromMolBlock(molecule, sanitize = True) if ("M  END" in molecule) else Chem.MolFromSmiles(molecule, sanitize = True)
    mutations = []
    if mol is not None:
        Chem.SanitizeMol(mol, sanitizeOps=Chem.rdmolops.SANITIZE_KEKULIZE)
        for step in range(0, steps):
            step_mutations = []
            for mol in [mol] if len(mutations) == 0 else mutations:
                step_mutations.extend(add_atoms(mol))
                step_mutations.extend(add_bonds(mol))
                step_mutations.extend(remove_bonds(mol))
                step_mutations = [x for x in step_mutations if x is not None]
                for mol in step_mutations:
                    Chem.SanitizeMol(mol, sanitizeOps=Chem.rdmolops.SANITIZE_KEKULIZE)
            mutations.extend(step_mutations)
        mutations = [Chem.MolToSmiles(molecule) for molecule in mutations]
    return list(set(mutations))


def get_random_mutations(molecule, steps=1, max_random_results=100):
    mol = Chem.MolFromMolBlock(molecule, sanitize = True) if ("M  END" in molecule) else Chem.MolFromSmiles(molecule, sanitize = True)
    mutations = []
    if mol is not None:
        Chem.SanitizeMol(mol, sanitizeOps=Chem.rdmolops.SANITIZE_KEKULIZE)
        while len(mutations) < max_random_results:
            _mol = deep_copy_mol(mol)
            for step in range(0, steps):
                mutation_type = random.randint(0, 2)
                if mutation_type == 0:
                    atoms = np.array(_mol.GetAtoms())
                    atoms = atoms[np.array([atom.GetImplicitValence() != 0 for atom in atoms], dtype=bool)]
                    if len(atoms) > 0:
                        _mol = add_atom(_mol, random.choice(atomic_numbers), random.choice(atoms))
                elif mutation_type == 1:
                    atoms = np.array(_mol.GetAtoms())
                    atoms = atoms[np.array([atom.GetImplicitValence() != 0 for atom in atoms], dtype=bool)]
                    if len(atoms) > 0:
                        a1 = random.choice(atoms)
                        atoms = np.array([_mol.GetAtomWithIdx(idx) for idx in range(a1.GetIdx() + 1, _mol.GetNumAtoms())])
                        atoms = atoms[np.array([atom.GetImplicitValence() != 0 for atom in atoms], dtype=bool)]
                        if len(atoms) > 0:
                            a2 = random.choice(atoms)
                            mols = add_bond_for_atoms(_mol, a1, a2)
                            if len(mols) > 0:
                                _mol = random.choice(mols)
                elif mutation_type == 2:
                    mols = remove_bond(_mol, random.choice(_mol.GetBonds()))
                    if len(mols) > 0:
                        _mol = random.choice(mols)
            Chem.SanitizeMol(_mol, sanitizeOps=Chem.rdmolops.SANITIZE_KEKULIZE)
            mutations.append(_mol)
    return [Chem.MolToSmiles(molecule) for molecule in mutations]

if len(molecule) < 3:
    molecule = "CCC"

if randomize:
    mutations = get_random_mutations(molecule, steps=steps, max_random_results=maxRandomResults)
else:
    mutations = get_all_mutations(molecule, steps=steps)

# Convert to Pandas DataFrame
mutations = pd.DataFrame(mutations, columns=['mutations'])
