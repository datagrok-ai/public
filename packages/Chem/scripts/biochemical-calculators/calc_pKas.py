#name: Calculate pKa
#description: Calculates acidic and basic pKa values, noting non-ionizable molecules.
#language: python
#meta.function_family: biochem-calculator
#environment: channels: [conda-forge, defaults], dependencies: [python=3.9, pip, rdkit, {pip: [pichemist]}]
#input: dataframe table
#input: column molecules {semType: Molecule}
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

import functools
import pandas as pd
from rdkit import Chem
from pichemist.api import pkas_and_charges_from_list
from pichemist.molecule import MolStandardiser, PeptideCutter
import pichemist.fasta.matcher as _fasta_matcher

# --- One-time speedup: cache the work pichemist's FASTA matcher throws away ---
# _pattern_match_rdkit re-parses each fragment (MolFromSmiles+AddHs) and recompiles
# each of the 186 amino-acid SMARTS on every call. Memoize both.
_compile_smarts = functools.lru_cache(maxsize=None)(Chem.MolFromSmarts)

@functools.lru_cache(maxsize=None)
def _mol_with_hs(smiles: str):
    m = Chem.MolFromSmiles(smiles)
    return Chem.AddHs(m) if m is not None else None

def _fast_pattern_match(smiles: str, smarts: str) -> int:
    mol = _mol_with_hs(smiles)
    if mol is None:
        return 0
    return len(mol.GetSubstructMatches(_compile_smarts(smarts)))

# matcher.py imported the name directly, so patch it there
_fasta_matcher.pattern_match = _fast_pattern_match

_standardiser = MolStandardiser()
_cutter = PeptideCutter()


def get_pka_values(mol_str: str, index: int) -> dict:
    """Calculates pKa values for a single molecule."""
    if not isinstance(mol_str, str) or not mol_str.strip():
        return {}
    try:
        mol = Chem.MolFromSmiles(mol_str) if 'M  END' not in mol_str else Chem.MolFromMolBlock(mol_str.replace('\\n', '\n'))
        if not mol:
            return {}
        # Same standardise + amide-bond fragmentation the full pipeline does,
        # then FASTA gate + pKa matcher only -- skips the unused pI/charge curve.
        mol = _standardiser.standardise_molecule(mol)
        frags = _cutter.break_amide_bonds_and_cap(mol)
        res = pkas_and_charges_from_list(frags, method="pkamatcher")
        base_pkas_calc, acid_pkas_calc = res[3], res[4]  # (pka, smiles) tuples
        return {
            'acidic': sorted(pka for pka, _ in acid_pkas_calc),
            'basic': sorted(pka for pka, _ in base_pkas_calc),
        }
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