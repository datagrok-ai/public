#name: Calculate pKa
#description: Calculates acidic and basic pKa values, noting non-ionizable molecules.
#language: python
#meta.function_family: biochem-calculator
#environment: channels: [conda-forge, defaults], dependencies: [python=3.9, pip, rdkit, {pip: [pichemist]}]
#input: dataframe table
#input: column molecules
#input: bool pKa_acidic_list = true {caption: Acidic List}
#input: bool pKa_basic_list = true {caption: Basic List}
#input: bool pKa_strongest_acidic = false {caption: Strongest Acidic}
#input: bool pKa_strongest_basic = false {caption: Strongest Basic}
#meta.role: transform
#meta.method_info.author: Ghiandoni G.M, Frolov A.
#meta.method_info.year: 2024
#meta.method_info.package: bio-pichemist-env
#meta.method_info.github: https://github.com/AstraZeneca/peptide-tools
#output: dataframe result {action:join(table)}


import pandas as pd
from rdkit import Chem
from pichemist.api import pichemist_from_dict
from pichemist.model import InputAttribute, OutputFragAttribute
import traceback

def get_pka_values(mol_str: str, index: int) -> dict:
    """Calculates pKa values for a single molecule."""
    if not isinstance(mol_str, str) or not mol_str.strip():
        return {}
    try:
        mol = Chem.MolFromSmiles(mol_str) if 'M  END' not in mol_str else Chem.MolFromMolBlock(mol_str.replace('\\n', '\n'))
        if not mol: return {}
        input_dict = {0: {InputAttribute.MOL_OBJECT.value: mol, InputAttribute.MOL_NAME.value: f"mol_{index}"}}
        results = pichemist_from_dict(input_dict, method="pkamatcher", print_fragments=True)
        properties = results.get(0, {})
        acidic_pkas, basic_pkas = [], []
        for key in ['frag_acid_pkas_calc', 'frag_base_pkas_calc']:
            for frag_data in properties.get(key, {}).values():
                pka = frag_data.get(OutputFragAttribute.PKA)
                if pka is not None:
                    (acidic_pkas if 'acid' in key else basic_pkas).append(pka)
        return {'acidic': sorted(acidic_pkas), 'basic': sorted(basic_pkas)}
    except Exception as e:
        print(f"ERROR on pKa for row {index}: {e}")
        return {}

molecule_data = table[molecules]
pka_results = [get_pka_values(m, i) for i, m in enumerate(molecule_data)]
result_data = {}


if pKa_acidic_list:
    result_data['acidic_fragment(s)'] = [
        ", ".join(map(str, r.get('acidic', []))) if r.get('acidic') else "non-ionizable"
        for r in pka_results
    ]
if pKa_basic_list:
    result_data['basic_fragment(s)'] = [
        ", ".join(map(str, r.get('basic', []))) if r.get('basic') else "non-ionizable"
        for r in pka_results
    ]
if pKa_strongest_acidic:
    result_data['strongest_acidic_fragment'] = [
        min(r.get('acidic')) if r.get('acidic') else "non-ionizable"
        for r in pka_results
    ]
if pKa_strongest_basic:
    result_data['strongest_basic_fragment'] = [
        max(r.get('basic')) if r.get('basic') else "non-ionizable"
        for r in pka_results
    ]

result = pd.DataFrame(result_data, index=table.index)