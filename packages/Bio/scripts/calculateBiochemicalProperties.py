#name: calculateBiochemicalProperties
#description: Calculates pI, pKa, logP, and logD for a molecule column.
#language: python
#tags: bio, chemistry, pichemist, properties
#environment: bio-pichemist-env
#input: dataframe TABLE
#input: string MOLECULES_COLUMN_NAME
#input: bool CALCULATE_PI = false
#input: bool CALCULATE_LOGP = false
#input: bool CALCULATE_LOGD = false
#input: bool IPC2_PEPTIDE = false
#input: bool IPC_PEPTIDE = false
#input: bool PROMOST = false
#input: bool GAUCI = false
#input: bool GRIMSLEY = false
#input: bool THURLKILL = false
#input: bool LEHNINGER = false
#input: bool TOSELAND = false
#output: dataframe result_df {action:join(TABLE)}

import math
import pandas as pd
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
from rdkit import Chem
from rdkit.Chem import Descriptors
from pichemist.api import pichemist_from_dict
from pichemist.model import InputAttribute, OutputFragAttribute

@dataclass
class MolecularAnalysisResult:
    mol_name: str = ""
    smiles: str = ""
    pi_mean: Optional[float] = None
    ipc2_peptide_pi: Optional[float] = None
    ipc_peptide_pi: Optional[float] = None
    promost_pi: Optional[float] = None
    gauci_pi: Optional[float] = None
    grimsley_pi: Optional[float] = None
    thurlkill_pi: Optional[float] = None
    lehninger_pi: Optional[float] = None
    toseland_pi: Optional[float] = None
    logp: Optional[float] = None
    logd: Optional[float] = None
    success: bool = True
    error_message: str = ""

def calculate_logd(logp: float, acidic_pkas: List[float], basic_pkas: List[float], ph: float = 7.4) -> float:
    denominator = 1.0
    for pka in acidic_pkas: denominator += 10**(ph - pka)
    for pka in basic_pkas: denominator += 10**(pka - ph)
    return logp - math.log10(denominator) if denominator > 0 else logp

def extract_pka_values(properties: Dict) -> tuple[List[float], List[float]]:
    acidic_pkas, basic_pkas = [], []
    for key in ['frag_acid_pkas_calc', 'frag_base_pkas_calc']:
        for frag_data in properties.get(key, {}).values():
            pka = frag_data.get(OutputFragAttribute.PKA)
            if pka is not None:
                (acidic_pkas if 'acid' in key else basic_pkas).append(pka)
    return acidic_pkas, basic_pkas

def analyze_molecule(
    molecule_string: str, mol_name: str, calculate_pi_flag: bool,
    pka_sets_to_use: List[str], calculate_logp_flag: bool, calculate_logd_flag: bool
) -> MolecularAnalysisResult:
    result = MolecularAnalysisResult(mol_name=mol_name)
    try:
        mol = Chem.MolFromSmiles(molecule_string) if 'M  END' not in molecule_string else Chem.MolFromMolBlock(molecule_string, strictParsing=False)
        if mol is None: raise ValueError("Failed to parse molecule")
        result.smiles = Chem.MolToSmiles(mol)
        if calculate_logp_flag or calculate_logd_flag:
            result.logp = Descriptors.MolLogP(mol)
    except Exception as e:
        result.success = False; result.error_message = str(e); return result

    try:
        pichemist_results = pichemist_from_dict(
            {0: {InputAttribute.MOL_NAME.value: mol_name, InputAttribute.MOL_OBJECT.value: mol}},
            method="pkamatcher", print_fragments=calculate_logd_flag
        )
        properties = pichemist_results.get(0, {})
    except Exception as e:
        result.success = False; result.error_message = f"Pichemist failed: {e}"; return result

    if calculate_pi_flag:
        pi_data = properties.get('pI', {})
        result.pi_mean = pi_data.get('pI mean')
        if pka_sets_to_use:
            pka_map = {
                'IPC2_peptide': 'ipc2_peptide_pi', 'IPC_peptide': 'ipc_peptide_pi', 'ProMoST': 'promost_pi',
                'Gauci': 'gauci_pi', 'Grimsley': 'grimsley_pi', 'Thurlkill': 'thurlkill_pi',
                'Lehninger': 'lehninger_pi', 'Toseland': 'toseland_pi'
            }
            for key, value in pi_data.items():
                if key in pka_sets_to_use:
                    setattr(result, pka_map[key], value)

    if calculate_logd_flag and result.logp is not None:
        acidic_pkas, basic_pkas = extract_pka_values(properties)
        result.logd = calculate_logd(result.logp, acidic_pkas, basic_pkas)

    return result

mols_list = TABLE[MOLECULES_COLUMN_NAME].tolist()

pka_methods = {
    'IPC2_peptide': IPC2_PEPTIDE, 'IPC_peptide': IPC_PEPTIDE, 'ProMoST': PROMOST,
    'Gauci': GAUCI, 'Grimsley': GRIMSLEY, 'Thurlkill': THURLKILL,
    'Lehninger': LEHNINGER, 'Toseland': TOSELAND
}
pka_sets_to_run = [key for key, value in pka_methods.items() if value]
run_pi_calculation = CALCULATE_PI or len(pka_sets_to_run) > 0

all_results = [analyze_molecule(
    mol_string, f"molecule_{i}", run_pi_calculation, pka_sets_to_run,
    CALCULATE_LOGP, CALCULATE_LOGD
) for i, mol_string in enumerate(mols_list)]

output_rows = [res.__dict__ for res in all_results]
result_df = pd.DataFrame(output_rows)

final_columns = []
if CALCULATE_PI: final_columns.append('pi_mean')
if CALCULATE_LOGP: final_columns.append('logp')
if CALCULATE_LOGD: final_columns.append('logd')
for pka_set in pka_sets_to_run:
    final_columns.append(f"{pka_set.lower()}_pi")

result_df = result_df[[col for col in final_columns if col in result_df.columns]]