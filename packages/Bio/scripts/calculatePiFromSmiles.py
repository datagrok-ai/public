# name: calculatePiFromSmiles
# description: Calculates pI from a column of SMILES strings.
# language: python
# tags: bio, chemistry, pichemist
# environment: bio-pichemist-env
# input: dataframe table
# input: string smiles_column_name "The name of the column with SMILES strings"
# output: dataframe result

import pandas as pd
from rdkit import Chem
from pichemist.api import pichemist_from_dict
from pichemist.model import InputAttribute

def _calculate_pi_single(smiles_string, mol_name='molecule'):
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            return {'pI_mean': None, 'pI_std': None, 'charge_at_pH7': None, 'error': 'Invalid SMILES'}
        input_dict = {0: {InputAttribute.MOL_NAME.value: mol_name, InputAttribute.MOL_OBJECT.value: mol}}
        results = pichemist_from_dict(input_dict, method="pkamatcher")
        result_dict = results[0]
        return {
            'pI_mean': result_dict['pI']['pI mean'],
            'pI_std': result_dict['pI']['std'],
            'charge_at_pH7': result_dict['QpH7']['Q at pH7.4 mean'],
            'error': None
        }
    except Exception as e:
        return {'pI_mean': None, 'pI_std': None, 'charge_at_pH7': None, 'error': str(e)}


results_list = []
smiles_column_data = table[smiles_column_name]
for i in table.index:
    smiles = smiles_column_data[i]
    results_list.append(_calculate_pi_single(smiles))
# ========================


result = pd.DataFrame(results_list)
result.index = table.index

result.rename(columns={
    'pI_mean': 'pI Mean',
    'pI_std': 'pI Std Dev',
    'charge_at_pH7': 'Charge @ pH 7.4',
    'error': 'pI Calculation Error'
}, inplace=True)