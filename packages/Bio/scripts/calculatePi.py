#name: calculatePi
#description: Calculates pI from a column of SMILES strings or MolBlocks.
#language: python
#tags: bio, chemistry, pichemist
#environment: bio-pichemist-env
#input: dataframe table
#input: string molecule_column_name "The name of the column with molecules"
#output: dataframe result {action:join(table)}


import pandas as pd
from rdkit import Chem
from pichemist.api import pichemist_from_dict
from pichemist.model import InputAttribute

if molecule_column_name not in table.columns:
    print(f"ERROR: Column '{molecule_column_name}' not found!")
    raise ValueError(f"Column not found: {molecule_column_name}")


def _calculate_pi_single(smiles_string, mol_name='molecule'):
    """Calculate pI for a single SMILES string."""
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
        print(f"Error processing SMILES '{smiles_string}': {str(e)}")
        return {'pI_mean': None, 'pI_std': None, 'charge_at_pH7': None}

def _molblock_to_smiles(molblock_string):
    """Convert MOL block to SMILES with debugging."""
    try:
        mol = Chem.MolFromMolBlock(molblock_string, strictParsing=False, removeHs=False)
        if mol is not None:
            smiles = Chem.MolToSmiles(mol)
            print(f"Successfully converted to SMILES: {smiles}")
            return smiles
        mol = Chem.MolFromMolBlock(molblock_string, strictParsing=False, removeHs=True)
        if mol is not None:
            smiles = Chem.MolToSmiles(mol)
            print(f"Converted to SMILES (removeHs=True): {smiles}")
            return smiles
        return None
        
    except Exception as e:
        print(f"Exception converting molblock: {str(e)}")
        return None

def _is_molblock(molecule_string):
    """Check if string is a MOL block by looking for characteristic markers."""
    if not isinstance(molecule_string, str):
        return False
    return 'M  END' in molecule_string or 'M END' in molecule_string

def _process_molecule(molecule_string, row_index):
    """Process a single molecule string (either SMILES or MOL block)."""
    if pd.isna(molecule_string) or molecule_string == '':
        return {'pI_mean': None, 'pI_std': None, 'charge_at_pH7': None}
    
    molecule_str = str(molecule_string)
    
    # Check if it's a MOL block and convert to SMILES
    if _is_molblock(molecule_str):
        print(f"Row {row_index}: Detected MOL block, converting...")
        smiles = _molblock_to_smiles(molecule_str)
        if smiles is None:
            print(f"Row {row_index}: Failed to convert MOL block to SMILES")
            return {'pI_mean': None, 'pI_std': None, 'charge_at_pH7': None}
    else:
        # Assume it's already a SMILES string
        print(f"Row {row_index}: Treating as SMILES: {molecule_str[:50]}...")
        smiles = molecule_str
    
    # Calculate pI from SMILES
    result = _calculate_pi_single(smiles, mol_name=f'molecule_{row_index}')
    if result['pI_mean'] is not None:
        print(f"Row {row_index}: pI = {result['pI_mean']:.2f}")
    else:
        print(f"Row {row_index}: pI calculation failed")
    
    return result

molecule_column_data = table[molecule_column_name]
results_list = []

# this needs to be vectorized of course
for i in table.index:
    molecule = molecule_column_data[i]
    result = _process_molecule(molecule, i)
    results_list.append(result)

result = pd.DataFrame(results_list)
result.index = table.index

result.rename(columns={
    'pI_mean': 'pI Mean',
    'pI_std': 'pI Std Dev',
    'charge_at_pH7': 'Charge @ pH 7.4'
}, inplace=True)
