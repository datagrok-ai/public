#name: calculatePiFromMolblock
#description: Calculates pI from a column of MolBlocks.
#language: python
#tags: bio, chemistry, pichemist
#environment: bio-pichemist-env
#input: dataframe table
#input: string molblock_column_name "The name of the column with MolBlocks"
#output: dataframe result

import pandas as pd
from rdkit import Chem
from pichemist.api import pichemist_from_dict
from pichemist.model import InputAttribute

def _calculate_pi_single(smiles_string, mol_name='molecule'):
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            return {'pI_mean': None, 'pI_std': None, 'charge_at_pH7': None}

        input_dict = {0: {InputAttribute.MOL_NAME.value: mol_name, InputAttribute.MOL_OBJECT.value: mol}}
        results = pichemist_from_dict(input_dict, method="pkamatcher")
        result_dict = results[0]
        
        return {
            'pI_mean': result_dict['pI']['pI mean'],
            'pI_std': result_dict['pI']['std'],
            'charge_at_pH7': result_dict['QpH7']['Q at pH7.4 mean']
        }
    except Exception as e:
        print(f"Error processing molecule: {str(e)}")
        return {'pI_mean': None, 'pI_std': None, 'charge_at_pH7': None}

def _calculate_pi_from_molblock_single(molblock_string, mol_name='molecule'):
    try:
        mol = Chem.MolFromMolBlock(molblock_string, strictParsing=False, removeHs=True)
        if mol is None:
            return {'pI_mean': None, 'pI_std': None, 'charge_at_pH7': None}
        
        smiles = Chem.MolToSmiles(mol)
        return _calculate_pi_single(smiles, mol_name=mol_name)
    except Exception as e:
        print(f"Error processing MolBlock: {str(e)}")
        return {'pI_mean': None, 'pI_std': None, 'charge_at_pH7': None}

# Doing each row.. How the hell do i vectorize this? is it automatic?
mol_results_list = []
molblock_column_data = table[molblock_column_name]

for i in table.index:
    molblock = molblock_column_data[i]
    if pd.isna(molblock) or molblock == '':
        mol_results_list.append({'pI_mean': None, 'pI_std': None, 'charge_at_pH7': None})
    else:
        mol_results_list.append(_calculate_pi_from_molblock_single(str(molblock)))

result = pd.DataFrame(mol_results_list)
result.index = table.index

result.rename(columns={
    'pI_mean': 'pI Mean',
    'pI_std': 'pI Std Dev',
    'charge_at_pH7': 'Charge @ pH 7.4'
}, inplace=True)