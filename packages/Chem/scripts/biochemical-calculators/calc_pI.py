#name: Calculate pI
#description: Calculates Isoelectric Point (pI) using various pKa datasets.
#language: python
#meta.function_family: biochem-calculator
#environment: channels: [conda-forge, defaults], dependencies: [python=3.9, pip, rdkit, {pip: [pichemist]}]
#input: dataframe table
#input: column molecules
#input: bool pI_mean = true {caption: Mean}
#input: bool pI_IPC2_peptide = false {caption: IPC2 Peptide}
#input: bool pI_IPC_peptide = false {caption: IPC Peptide}
#input: bool pI_ProMoST = false {caption: ProMoST}
#input: bool pI_Gauci = false {caption: Gauci}
#input: bool pI_Grimsley = false {caption: Grimsley}
#input: bool pI_Thurlkill = false {caption: Thurlkill}
#input: bool pI_Lehninger = false {caption: Lehninger}
#input: bool pI_Toseland = false {caption: Toseland}
#meta.role: transform
#meta.method_info.author: Kozlova, L., Garbuzynskiy, S.O.
#meta.method_info.year: 2023
#meta.method_info.package: bio-pichemist-env
#meta.method_info.github: https://github.com/AstraZeneca/peptide-tools
#output: dataframe result {action:join(table)}

import pandas as pd
from rdkit import Chem
from pichemist.api import pichemist_from_dict
from pichemist.model import InputAttribute
import traceback

def calculate_pi_for_molecule(mol_str: str, index: int) -> dict:
    """Calculates all pI values for a single molecule string."""
    if not isinstance(mol_str, str) or not mol_str.strip():
        return {}

    try:
        mol = None
        if 'M  END' in mol_str:
            mol = Chem.MolFromMolBlock(mol_str.replace('\\n', '\n'), sanitize=True)
        else:
            mol = Chem.MolFromSmiles(mol_str, sanitize=True)

        if not mol:
            return {}

        input_dict = {
            0: {
                InputAttribute.MOL_OBJECT.value: mol,
                InputAttribute.MOL_NAME.value: f"molecule_{index}"
            }
        }
        results = pichemist_from_dict(input_dict, method="pkamatcher")

        properties = results.get(0, {})
        pi_data = properties.get('pI', {})
        return pi_data

    except Exception as e:
        print(f"ERROR processing molecule at row {index}: {e}")
        traceback.print_exc()
        return {}

molecule_data = table[molecules]
pi_results = [calculate_pi_for_molecule(m, i) for i, m in enumerate(molecule_data)]


result_data = {}
pi_calculation_map = {
    'pI_mean': 'pI mean', 'pI_IPC2_peptide': 'IPC2_peptide', 'pI_IPC_peptide': 'IPC_peptide',
    'pI_ProMoST': 'ProMoST', 'pI_Gauci': 'Gauci', 'pI_Grimsley': 'Grimsley',
    'pI_Thurlkill': 'Thurlkill', 'pI_Lehninger': 'Lehninger', 'pI_Toseland': 'Toseland'
}

requested_calcs = {key: name for key, name in pi_calculation_map.items() if globals().get(key)}
if not requested_calcs and pI_mean:
    requested_calcs['pI_mean'] = pi_calculation_map['pI_mean']

for key, name_in_lib in requested_calcs.items():
    col_name = f"pI ({name_in_lib.replace('_', ' ').replace('pI', '').strip()})"
    result_data[col_name] = [res.get(name_in_lib) for res in pi_results]

result = pd.DataFrame(result_data, index=table.index)
