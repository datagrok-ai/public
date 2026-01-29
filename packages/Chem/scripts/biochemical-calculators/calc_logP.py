#name: Calculate logP
#description: Calculates the octanol-water partition coefficient (logP) using RDKit's Crippen method.
#language: python
#meta.function_family: biochem-calculator
#input: dataframe table
#input: column molecules
#meta.role: transform
#meta.method_info.author: RDKit Team
#meta.method_info.year: 2024
#meta.method_info.package: rdkit-base
#meta.method_info.github: https://github.com/rdkit/rdkit
#output: dataframe result {action:join(table)}

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

def calculate_logp(mol_str: str, index: int):
    """Calculates logP for a single molecule."""
    if not isinstance(mol_str, str) or not mol_str.strip():
        return None

    try:
        mol = Chem.MolFromSmiles(mol_str) if 'M  END' not in mol_str else Chem.MolFromMolBlock(mol_str.replace('\\n', '\n'))
        if mol:
            return Descriptors.MolLogP(mol)
        return None
    except Exception as e:
        print(f"ERROR on logP for row {index}: {e}")
        return None

molecule_data = table[molecules]
logp_values = [calculate_logp(m, i) for i, m in enumerate(molecule_data)]
result = pd.DataFrame({'logP': logp_values}, index=table.index)
