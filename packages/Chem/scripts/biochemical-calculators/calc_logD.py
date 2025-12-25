#name: Calculate logD
#description: Calculates the distribution coefficient (logD) at a specified pH using pKa and logP.
#language: python
#meta.function_family: biochem-calculator
#environment: channels: [conda-forge, defaults], dependencies: [python=3.9, pip, rdkit, {pip: [pichemist]}]
#input: dataframe table
#input: column molecules {caption: Molecules column}
#input: double pH = 7.4 {caption: pH}
#meta.role: transform
#meta.method_info.author: Datagrok Hybrid Method
#meta.method_info.year: 2024
#meta.method_info.package: https://github.com/datagrok-ai/public
#output: dataframe result {action:join(table)}

import pandas as pd
import math
from rdkit import Chem
from rdkit.Chem import Descriptors
from pichemist.api import pichemist_from_dict
from pichemist.model import InputAttribute, OutputFragAttribute

def calculate_logd(mol_str: str, pH: float, index: int) -> float:
    """Calculates logD for a single molecule."""
    if not isinstance(mol_str, str) or not mol_str.strip():
        return None

    try:
        mol = Chem.MolFromSmiles(mol_str) if 'M  END' not in mol_str else Chem.MolFromMolBlock(mol_str.replace('\\n', '\n'))
        if not mol: return None

        logp = Descriptors.MolLogP(mol)

        input_dict = {0: {InputAttribute.MOL_OBJECT.value: mol, InputAttribute.MOL_NAME.value: f"mol_{index}"}}
        results = pichemist_from_dict(input_dict, method="pkamatcher", print_fragments=True)

        properties = results.get(0, {})
        acidic_pkas, basic_pkas = [], []
        for key in ['frag_acid_pkas_calc', 'frag_base_pkas_calc']:
            for frag_data in properties.get(key, {}).values():
                pka = frag_data.get(OutputFragAttribute.PKA)
                if pka is not None:
                    (acidic_pkas if 'acid' in key else basic_pkas).append(pka)

        if not acidic_pkas and not basic_pkas: return logp # Neutral species

        # Use the pH input parameter in the calculation
        sum_acids = sum(10**(pH - pka) for pka in acidic_pkas)
        sum_bases = sum(10**(pka - pH) for pka in basic_pkas)
        
        # logD equation for acids and bases
        denominator = 1.0 + sum_acids + sum_bases

        return logp - math.log10(denominator)
    except Exception as e:
        print(f"ERROR on logD for row {index}: {e}")
        return None


molecule_data = table[molecules]
logd_values = [calculate_logd(m, pH, i) for i, m in enumerate(molecule_data)]
result = pd.DataFrame({f'logD @ pH {pH}': logd_values}, index=table.index)