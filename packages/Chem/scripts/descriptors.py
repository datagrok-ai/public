#name: Desc
#language: python
#input: string smiles
#input: dataframe df1
#input: string selected
#input: dataframe df2
#output: dataframe out

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
from rdkit.Chem import rdDepictor
from rdkit.Chem import Descriptors
from rdkit.Chem import Descriptors3D
from rdkit.Chem.Draw import rdMolDraw2D
import numpy as np
import pandas as pd
import multiprocessing

smiles = df1[smiles]
descriptors = df2[selected]


_descriptors_tree = {}

descriptors_modules = {
    'Crippen': {
        'name': 'Crippen',
        'description': "Atom-based calculation of LogP and MR using Crippen's approach",
    },

    'EState.EState_VSA': {
        'name': 'EState VSA',
        'description': 'Hybrid EState-VSA descriptors (like the MOE VSA descriptors) ',
    },

    'QED': {
        'description': 'QED stands for quantitative estimation of drug-likeness',
    },

    'GraphDescriptors': {
        'description': 'Topological/topochemical descriptors ',
    },

    'Descriptors': {
        'description': 'General descriptors ',
    },

    'MolSurf': {
        'description': 'MOE-like approximate molecular surface area descriptors',
    },

    'EState.EState': {
        'name': 'EState',
        'description': 'Basic EState descriptors',
    },

    'Fragments': {
        'description': 'Bunch of fragment descriptors from a file ',
    },

    'Lipinski': {
        'description': 'Lipinski parameters for molecules ',
    },

    'Descriptors3D': {
        'name': 'Descriptors 3D',
        'description': "Descriptors derived from a molecule's 3D structure"
    }
}


descriptors_params = {

    'EState_VSA1': {
        'type': 'double',
        'description': 'EState VSA Descriptor 1 (-inf < x < -0.39)',
        'tags': ['2D']
    },

    'EState_VSA10': {
        'type': 'double',
        'description': 'EState VSA Descriptor 10 ( 9.17 <= x < 15.00)',
        'tags': ['2D']
    },

    'EState_VSA11': {
        'type': 'double',
        'description': 'EState VSA Descriptor 11 ( 15.00 <= x < inf)',
        'tags': ['2D']
    },

    'EState_VSA2': {
        'type': 'double',
        'description': 'EState VSA Descriptor 2 ( -0.39 <= x < 0.29)',
        'tags': ['2D']
    },

    'EState_VSA3': {
        'type': 'double',
        'description': 'EState VSA Descriptor 3 ( 0.29 <= x < 0.72)',
        'tags': ['2D']
    },

    'EState_VSA4': {
        'type': 'double',
        'description': 'EState VSA Descriptor 4 ( 0.72 <= x < 1.17)',
        'tags': ['2D']
    },

    'EState_VSA5': {
        'type': 'double',
        'description': 'EState VSA Descriptor 5 ( 1.17 <= x < 1.54)',
        'tags': ['2D']
    },

    'EState_VSA6': {
        'type': 'double',
        'description': 'EState VSA Descriptor 6 ( 1.54 <= x < 1.81)',
        'tags': ['2D']
    },

    'EState_VSA7': {
        'type': 'double',
        'description': 'EState VSA Descriptor 7 ( 1.81 <= x < 2.05)',
        'tags': ['2D']
    },

    'EState_VSA8': {
        'type': 'double',
        'description': 'EState VSA Descriptor 8 ( 2.05 <= x < 4.69)',
        'tags': ['2D']
    },

    'EState_VSA9': {
        'type': 'double',
        'description': 'EState VSA Descriptor 9 ( 4.69 <= x < 9.17)',
        'tags': ['2D']
    },

    'VSA_EState1': {
        'type': 'double',
        'description': 'VSA EState Descriptor 1 (-inf < x < 4.78)',
        'tags': ['2D']
    },

    'VSA_EState10': {
        'type': 'double',
        'description': 'VSA EState Descriptor 10 ( 11.00 <= x < inf)',
        'tags': ['2D']
    },

    'VSA_EState2': {
        'type': 'double',
        'description': 'VSA EState Descriptor 2 ( 4.78 <= x < 5.00)',
        'tags': ['2D']
    },

    'VSA_EState3': {
        'type': 'double',
        'description': 'VSA EState Descriptor 3 ( 5.00 <= x < 5.41)',
        'tags': ['2D']
    },

    'VSA_EState4': {
        'type': 'double',
        'description': 'VSA EState Descriptor 4 ( 5.41 <= x < 5.74)',
        'tags': ['2D']
    },

    'VSA_EState5': {
        'type': 'double',
        'description': 'VSA EState Descriptor 5 ( 5.74 <= x < 6.00)',
        'tags': ['2D']
    },

    'VSA_EState6': {
        'type': 'double',
        'description': 'VSA EState Descriptor 6 ( 6.00 <= x < 6.07)',
        'tags': ['2D']
    },

    'VSA_EState7': {
        'type': 'double',
        'description': 'VSA EState Descriptor 7 ( 6.07 <= x < 6.45)',
        'tags': ['2D']
    },

    'VSA_EState8': {
        'type': 'double',
        'description': 'VSA EState Descriptor 8 ( 6.45 <= x < 7.00)',
        'tags': ['2D']
    },

    'VSA_EState9': {
        'type': 'double',
        'description': 'VSA EState Descriptor 9 ( 7.00 <= x < 11.00)',
        'tags': ['2D']
    },

    'qed': {
        'type': 'double',
        'description': 'Weighted sum of ADS mapped properties',
        'tags': ['2D']
    },

    'BalabanJ': {
        'type': 'double',
        'description': 'Balaban\'s J value',
        'tags': ['2D']
    },

    'BertzCT': {
        'type': 'double',
        'description': 'A topological index meant to quantify "complexity"',
        'tags': ['2D']
    },

    'Chi0': {
        'type': 'double',
        'description': 'From equations (1),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)',
        'tags': ['2D']
    },

    'Chi0n': {
        'type': 'double',
        'description': 'Similar to Hall Kier Chi0v, but uses nVal instead of valence. '
                       'This makes a big difference after we get out of the first row',
        'tags': ['2D']
    },

    'Chi0v': {
        'type': 'double',
        'description': 'From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)',
        'tags': ['2D']
    },

    'Chi1': {
        'type': 'double',
        'description': 'From equations (1),(11) and (12) of Rev. Comp. Chem. vol 2, 367-422, (1991)',
        'tags': ['2D']
    },

    'Chi1n': {
        'type': 'double',
        'description': 'Similar to Hall Kier Chi1v, but uses nVal instead of valence',
        'tags': ['2D']
    },

    'Chi1v': {
        'type': 'double',
        'description': 'From equations (5),(11) and (12) of Rev. Comp. Chem. vol 2, 367-422, (1991)',
        'tags': ['2D']
    },

    'Chi2n': {
        'type': 'double',
        'description': 'Similar to Hall Kier Chi2v, but uses nVal instead of valence. '
                       ' This makes a big difference after we get out of the first row',
        'tags': ['2D']
    },

    'Chi2v': {
        'type': 'double',
        'description': 'From equations (5),(15) and (16) of Rev. Comp. Chem. vol 2, 367-422, (1991)',
        'tags': ['2D']
    },

    'Chi3n': {
        'type': 'double',
        'description': 'Similar to Hall Kier Chi3v, but uses nVal instead of valence. '
                       'This makes a big difference after we get out of the first row',
        'tags': ['2D']
    },

    'Chi3v': {
        'type': 'double',
        'description': 'From equations (5),(15) and (16) of Rev. Comp. Chem. vol 2, 367-422, (1991)',
        'tags': ['2D']
    },

    'Chi4n': {
        'type': 'double',
        'description': 'Similar to Hall Kier Chi4v, but uses nVal instead of valence. '
                       'This makes a big difference after we get out of the first row.\n\n'
                       '**NOTE**: because the current path finding code does, by design, detect rings '
                       'as paths (e.g. in C1CC1 there is *1* atom path of length 3), values of Chi4v '
                       'may give results that differ from those provided by the old code in molecules '
                       'that have 3 rings',
        'tags': ['2D']
    },

    'Chi4v': {
        'type': 'double',
        'description': 'From equations (5),(15) and (16) of Rev. Comp. Chem. vol 2, 367-422, (1991).\n\n'
                       '**NOTE**: because the current path finding code does, by design, detect rings '
                       'as paths (e.g. in C1CC1 there is *1* atom path of length 3), values of Chi4v '
                       'may give results that differ from those provided by the old code in molecules '
                       'that have 3 rings',
        'tags': ['2D']
    },

    'HallKierAlpha': {
        'type': 'double',
        'description': 'The Hall-Kier alpha value for a molecule',
        'tags': ['2D']
    },

    'Ipc': {
        'type': 'double',
        'description': 'The information content of the coefficients of the characteristic polynomial of '
                       'the adjacency matrix of a hydrogen-suppressed graph of a molecule',
        'tags': ['2D']
    },

    'Kappa1': {
        'type': 'double',
        'description': 'Hall-Kier Kappa1 value',
        'tags': ['2D']
    },

    'Kappa2': {
        'type': 'double',
        'description': 'Hall-Kier Kappa2 value',
        'tags': ['2D']
    },

    'Kappa3': {
        'type': 'double',
        'description': 'Hall-Kier Kappa3 value',
        'tags': ['2D']
    },

    'MinAbsPartialCharge': {
        'type': 'double',
        'description': 'Minimal absolute partial charge',
        'tags': ['2D']
    },

    'NumRadicalElectrons': {
        'type': 'int',
        'description': 'The number of radical electrons the molecule has (says nothing about spin state)',
        'tags': ['2D']
    },

    'FpDensityMorgan2': {
        'type': 'double',
        'description': 'Morgan fingerprint, radius 2',
        'tags': ['2D']
    },

    'FpDensityMorgan3': {
        'type': 'double',
        'description': 'Morgan fingerprint, radius 3',
        'tags': ['2D']
    },

    'FpDensityMorgan1': {
        'type': 'double',
        'description': 'Morgan fingerprint, radius 1',
        'tags': ['2D']
    },

    'HeavyAtomMolWt': {
        'type': 'double',
        'description': 'The average molecular weight of the molecule ignoring hydrogens',
        'tags': ['2D']
    },

    'MaxAbsPartialCharge': {
        'type': 'double',
        'description': 'Maximum absolute partial charge',
        'tags': ['2D']
    },

    'MinPartialCharge': {
        'type': 'double',
        'description': 'Minimal partial charge',
        'tags': ['2D']
    },

    'ExactMolWt': {
        'type': 'double',
        'description': 'The exact molecular weight of the molecule',
        'tags': ['2D']
    },

    'MolWt': {
        'type': 'double',
        'description': 'The average molecular weight of the molecule',
        'tags': ['2D']
    },

    'NumValenceElectrons': {
        'type': 'int',
        'description': 'The number of valence electrons the molecule has',
        'tags': ['2D']
    },

    'MaxPartialCharge': {
        'type': 'double',
        'description': 'Maximum partial charge',
        'tags': ['2D']
    },

    'MolLogP': {
        'type': 'double',
        'description': 'Wildman-Crippen LogP value',
        'tags': ['2D']
    },

    'MolMR': {
        'type': 'double',
        'description': 'Wildman-Crippen MR value',
        'tags': ['2D']
    },

    'LabuteASA': {
        'type': 'double',
        'description': 'Labute\'s Approximate Surface Area (ASA from MOE)',
        'tags': ['2D']
    },

    'PEOE_VSA1': {
        'type': 'double',
        'description': 'MOE Charge VSA Descriptor 1 (-inf < x < -0.30)',
        'tags': ['2D']
    },

    'PEOE_VSA10': {
        'type': 'double',
        'description': 'MOE Charge VSA Descriptor 10 ( 0.10 <= x < 0.15)',
        'tags': ['2D']
    },

    'PEOE_VSA11': {
        'type': 'double',
        'description': 'MOE Charge VSA Descriptor 11 ( 0.15 <= x < 0.20)',
        'tags': ['2D']
    },

    'PEOE_VSA12': {
        'type': 'double',
        'description': 'MOE Charge VSA Descriptor 12 ( 0.20 <= x < 0.25)',
        'tags': ['2D']
    },

    'PEOE_VSA13': {
        'type': 'double',
        'description': 'MOE Charge VSA Descriptor 13 ( 0.25 <= x < 0.30)',
        'tags': ['2D']
    },

    'PEOE_VSA14': {
        'type': 'double',
        'description': 'MOE Charge VSA Descriptor 14 ( 0.30 <= x < inf)',
        'tags': ['2D']
    },

    'PEOE_VSA2': {
        'type': 'double',
        'description': 'MOE Charge VSA Descriptor 2 (-0.30 <= x < -0.25)',
        'tags': ['2D']
    },

    'PEOE_VSA3': {
        'type': 'double',
        'description': 'MOE Charge VSA Descriptor 3 (-0.25 <= x < -0.20)',
        'tags': ['2D']
    },

    'PEOE_VSA4': {
        'type': 'double',
        'description': 'MOE Charge VSA Descriptor 4 (-0.20 <= x < -0.15)',
        'tags': ['2D']
    },

    'PEOE_VSA5': {
        'type': 'double',
        'description': 'MOE Charge VSA Descriptor 5 (-0.15 <= x < -0.10)',
        'tags': ['2D']
    },

    'PEOE_VSA6': {
        'type': 'double',
        'description': 'MOE Charge VSA Descriptor 6 (-0.10 <= x < -0.05)',
        'tags': ['2D']
    },

    'PEOE_VSA7': {
        'type': 'double',
        'description': 'MOE Charge VSA Descriptor 7 (-0.05 <= x < 0.00)',
        'tags': ['2D']
    },

    'PEOE_VSA8': {
        'type': 'double',
        'description': 'MOE Charge VSA Descriptor 8 ( 0.00 <= x < 0.05)',
        'tags': ['2D']
    },

    'PEOE_VSA9': {
        'type': 'double',
        'description': 'MOE Charge VSA Descriptor 9 ( 0.05 <= x < 0.10)',
        'tags': ['2D']
    },

    'SMR_VSA1': {
        'type': 'double',
        'description': 'MOE MR VSA Descriptor 1 (-inf < x < 1.29)',
        'tags': ['2D']
    },

    'SMR_VSA10': {
        'type': 'double',
        'description': 'MOE MR VSA Descriptor 10 ( 4.00 <= x < inf)',
        'tags': ['2D']
    },

    'SMR_VSA2': {
        'type': 'double',
        'description': 'MOE MR VSA Descriptor 2 ( 1.29 <= x < 1.82)',
        'tags': ['2D']
    },

    'SMR_VSA3': {
        'type': 'double',
        'description': 'MOE MR VSA Descriptor 3 ( 1.82 <= x < 2.24)',
        'tags': ['2D']
    },

    'SMR_VSA4': {
        'type': 'double',
        'description': 'MOE MR VSA Descriptor 4 ( 2.24 <= x < 2.45)',
        'tags': ['2D']
    },

    'SMR_VSA5': {
        'type': 'double',
        'description': 'MOE MR VSA Descriptor 5 ( 2.45 <= x < 2.75)',
        'tags': ['2D']
    },

    'SMR_VSA6': {
        'type': 'double',
        'description': 'MOE MR VSA Descriptor 6 ( 2.75 <= x < 3.05)',
        'tags': ['2D']
    },

    'SMR_VSA7': {
        'type': 'double',
        'description': 'MOE MR VSA Descriptor 7 ( 3.05 <= x < 3.63)',
        'tags': ['2D']
    },

    'SMR_VSA8': {
        'type': 'double',
        'description': 'MOE MR VSA Descriptor 8 ( 3.63 <= x < 3.80)',
        'tags': ['2D']
    },

    'SMR_VSA9': {
        'type': 'double',
        'description': 'MOE MR VSA Descriptor 9 ( 3.80 <= x < 4.00)',
        'tags': ['2D']
    },

    'SlogP_VSA1': {
        'type': 'double',
        'description': 'MOE logP VSA Descriptor 1 (-inf < x < -0.40)',
        'tags': ['2D']
    },

    'SlogP_VSA10': {
        'type': 'double',
        'description': 'MOE logP VSA Descriptor 10 ( 0.40 <= x < 0.50)',
        'tags': ['2D']
    },

    'SlogP_VSA11': {
        'type': 'double',
        'description': 'MOE logP VSA Descriptor 11 ( 0.50 <= x < 0.60)',
        'tags': ['2D']
    },

    'SlogP_VSA12': {
        'type': 'double',
        'description': 'MOE logP VSA Descriptor 12 ( 0.60 <= x < inf)',
        'tags': ['2D']
    },

    'SlogP_VSA2': {
        'type': 'double',
        'description': 'MOE logP VSA Descriptor 2 (-0.40 <= x < -0.20)',
        'tags': ['2D']
    },

    'SlogP_VSA3': {
        'type': 'double',
        'description': 'MOE logP VSA Descriptor 3 (-0.20 <= x < 0.00)',
        'tags': ['2D']
    },

    'SlogP_VSA4': {
        'type': 'double',
        'description': 'MOE logP VSA Descriptor 4 ( 0.00 <= x < 0.10)',
        'tags': ['2D']
    },

    'SlogP_VSA5': {
        'type': 'double',
        'description': 'MOE logP VSA Descriptor 5 ( 0.10 <= x < 0.15)',
        'tags': ['2D']
    },

    'SlogP_VSA6': {
        'type': 'double',
        'description': 'MOE logP VSA Descriptor 6 ( 0.15 <= x < 0.20)',
        'tags': ['2D']
    },

    'SlogP_VSA7': {
        'type': 'double',
        'description': 'MOE logP VSA Descriptor 7 ( 0.20 <= x < 0.25)',
        'tags': ['2D']
    },

    'SlogP_VSA8': {
        'type': 'double',
        'description': 'MOE logP VSA Descriptor 8 ( 0.25 <= x < 0.30)',
        'tags': ['2D']
    },

    'SlogP_VSA9': {
        'type': 'double',
        'description': 'MOE logP VSA Descriptor 9 ( 0.30 <= x < 0.40)',
        'tags': ['2D']
    },

    'TPSA': {
        'type': 'double',
        'description': 'The polar surface area of a molecule based upon fragments.',
        'tags': ['2D']
    },

    'MaxAbsEStateIndex': {
        'type': 'double',
        'description': 'Maximum absolute EState index',
        'tags': ['2D']
    },

    'MaxEStateIndex': {
        'type': 'double',
        'description': 'Maximum EState index',
        'tags': ['2D']
    },

    'MinEStateIndex': {
        'type': 'double',
        'description': 'Minimum EState index',
        'tags': ['2D']
    },

    'MinAbsEStateIndex': {
        'type': 'double',
        'description': 'Minimum absolute EState index',
        'tags': ['2D']
    },

    'fr_Al_COO': {
        'type': 'int',
        'description': 'Number of aliphatic carboxylic acids',
        'tags': ['2D']
    },

    'fr_Al_OH': {
        'type': 'int',
        'description': 'Number of aliphatic hydroxyl groups',
        'tags': ['2D']
    },

    'fr_Al_OH_noTert': {
        'type': 'int',
        'description': 'Number of aliphatic hydroxyl groups excluding tert-OH',
        'tags': ['2D']
    },

    'fr_ArN': {
        'type': 'int',
        'description': 'Number of N functional groups attached to aromatics',
        'tags': ['2D']
    },

    'fr_Ar_COO': {
        'type': 'int',
        'description': 'Number of Aromatic carboxylic acide',
        'tags': ['2D']
    },

    'fr_Ar_N': {
        'type': 'int',
        'description': 'Number of aromatic nitrogens',
        'tags': ['2D']
    },

    'fr_Ar_NH': {
        'type': 'int',
        'description': 'Number of aromatic amines',
        'tags': ['2D']
    },

    'fr_Ar_OH': {
        'type': 'int',
        'description': 'Number of aromatic hydroxyl groups',
        'tags': ['2D']
    },

    'fr_COO': {
        'type': 'int',
        'description': 'Number of carboxylic acids',
        'tags': ['2D']
    },

    'fr_COO2': {
        'type': 'int',
        'description': 'Number of carboxylic acids',
        'tags': ['2D']
    },

    'fr_C_O': {
        'type': 'int',
        'description': 'Number of carbonyl O',
        'tags': ['2D']
    },

    'fr_C_O_noCOO': {
        'type': 'int',
        'description': 'Number of carbonyl O, excluding COOH',
        'tags': ['2D']
    },

    'fr_C_S': {
        'type': 'int',
        'description': 'Number of thiocarbonyl',
        'tags': ['2D']
    },

    'fr_HOCCN': {
        'type': 'int',
        'description': 'Number of C(OH)CCN-Ctert-alkyl or C(OH)CCNcyclic',
        'tags': ['2D']
    },

    'fr_Imine': {
        'type': 'int',
        'description': 'Number of Imines',
        'tags': ['2D']
    },

    'fr_NH0': {
        'type': 'int',
        'description': 'Number of Tertiary amines',
        'tags': ['2D']
    },

    'fr_NH1': {
        'type': 'int',
        'description': 'Number of Secondary amines',
        'tags': ['2D']
    },

    'fr_NH2': {
        'type': 'int',
        'description': 'Number of Primary amines',
        'tags': ['2D']
    },

    'fr_N_O': {
        'type': 'int',
        'description': 'Number of hydroxylamine groups',
        'tags': ['2D']
    },

    'fr_Ndealkylation1': {
        'type': 'int',
        'description': 'Number of XCCNR groups',
        'tags': ['2D']
    },

    'fr_Ndealkylation2': {
        'type': 'int',
        'description': 'Number of tert-alicyclic amines (no heteroatoms, not quinine-like bridged N)',
        'tags': ['2D']
    },

    'fr_Nhpyrrole': {
        'type': 'int',
        'description': 'Number of H-pyrrole nitrogens',
        'tags': ['2D']
    },

    'fr_SH': {
        'type': 'int',
        'description': 'Number of thiol groups',
        'tags': ['2D']
    },

    'fr_aldehyde': {
        'type': 'int',
        'description': 'Number of aldehydes',
        'tags': ['2D']
    },

    'fr_alkyl_carbamate': {
        'type': 'int',
        'description': 'Number of alkyl carbamates (subject to hydrolysis)',
        'tags': ['2D']
    },

    'fr_alkyl_halide': {
        'type': 'int',
        'description': 'Number of alkyl halides',
        'tags': ['2D']
    },

    'fr_allylic_oxid': {
        'type': 'int',
        'description': 'Number of allylic oxidation sites excluding steroid dienone',
        'tags': ['2D']
    },

    'fr_amide': {
        'type': 'int',
        'description': 'Number of amides',
        'tags': ['2D']
    },

    'fr_amidine': {
        'type': 'int',
        'description': 'Number of amidine groups',
        'tags': ['2D']
    },

    'fr_aniline': {
        'type': 'int',
        'description': 'Number of anilines',
        'tags': ['2D']
    },

    'fr_aryl_methyl': {
        'type': 'int',
        'description': 'Number of aryl methyl sites for hydroxylation',
        'tags': ['2D']
    },

    'fr_azide': {
        'type': 'int',
        'description': 'Number of azide groups',
        'tags': ['2D']
    },

    'fr_azo': {
        'type': 'int',
        'description': 'Number of azo groups',
        'tags': ['2D']
    },

    'fr_barbitur': {
        'type': 'int',
        'description': 'Number of barbiturate groups',
        'tags': ['2D']
    },

    'fr_benzene': {
        'type': 'int',
        'description': 'Number of benzene rings',
        'tags': ['2D']
    },

    'fr_benzodiazepine': {
        'type': 'int',
        'description': 'Number of benzodiazepines with no additional fused rings',
        'tags': ['2D']
    },

    'fr_bicyclic': {
        'type': 'int',
        'description': 'Bicyclic',
        'tags': ['2D']
    },

    'fr_diazo': {
        'type': 'int',
        'description': 'Number of diazo groups',
        'tags': ['2D']
    },

    'fr_dihydropyridine': {
        'type': 'int',
        'description': 'Number of dihydropyridines',
        'tags': ['2D']
    },

    'fr_epoxide': {
        'type': 'int',
        'description': 'Number of epoxide rings',
        'tags': ['2D']
    },

    'fr_ester': {
        'type': 'int',
        'description': 'Number of esters',
        'tags': ['2D']
    },

    'fr_ether': {
        'type': 'int',
        'description': 'Number of ether oxygens (including phenoxy)',
        'tags': ['2D']
    },

    'fr_furan': {
        'type': 'int',
        'description': 'Number of furan rings',
        'tags': ['2D']
    },

    'fr_guanido': {
        'type': 'int',
        'description': 'Number of guanidine groups',
        'tags': ['2D']
    },

    'fr_halogen': {
        'type': 'int',
        'description': 'Number of halogens',
        'tags': ['2D']
    },

    'fr_hdrzine': {
        'type': 'int',
        'description': 'Number of hydrazine groups',
        'tags': ['2D']
    },

    'fr_hdrzone': {
        'type': 'int',
        'description': 'Number of hydrazone groups',
        'tags': ['2D']
    },

    'fr_imidazole': {
        'type': 'int',
        'description': 'Number of imidazole rings',
        'tags': ['2D']
    },

    'fr_imide': {
        'type': 'int',
        'description': 'Number of imide groups',
        'tags': ['2D']
    },

    'fr_isocyan': {
        'type': 'int',
        'description': 'Number of isocyanates',
        'tags': ['2D']
    },

    'fr_isothiocyan': {
        'type': 'int',
        'description': 'Number of isothiocyanates',
        'tags': ['2D']
    },

    'fr_ketone': {
        'type': 'int',
        'description': 'Number of ketones',
        'tags': ['2D']
    },

    'fr_ketone_Topliss': {
        'type': 'int',
        'description': 'Number of ketones excluding diaryl, a,b-unsat',
        'tags': ['2D']
    },

    'fr_lactam': {
        'type': 'int',
        'description': 'Number of beta lactams',
        'tags': ['2D']
    },

    'fr_lactone': {
        'type': 'int',
        'description': 'Number of cyclic esters (lactones)',
        'tags': ['2D']
    },

    'fr_methoxy': {
        'type': 'int',
        'description': 'Number of methoxy groups -OCH3',
        'tags': ['2D']
    },

    'fr_morpholine': {
        'type': 'int',
        'description': 'Number of morpholine rings',
        'tags': ['2D']
    },

    'fr_nitrile': {
        'type': 'int',
        'description': 'Number of nitriles',
        'tags': ['2D']
    },

    'fr_nitro': {
        'type': 'int',
        'description': 'Number of nitro groups',
        'tags': ['2D']
    },

    'fr_nitro_arom': {
        'type': 'int',
        'description': 'Number of nitro benzene ring substituent',
        'tags': ['2D']
    },

    'fr_nitro_arom_nonortho': {
        'type': 'int',
        'description': 'Number of non-ortho nitro benzene ring substituents',
        'tags': ['2D']
    },

    'fr_nitroso': {
        'type': 'int',
        'description': 'Number of nitroso groups, excluding NO2',
        'tags': ['2D']
    },

    'fr_oxazole': {
        'type': 'int',
        'description': 'Number of oxazole rings',
        'tags': ['2D']
    },

    'fr_oxime': {
        'type': 'int',
        'description': 'Number of oxime groups',
        'tags': ['2D']
    },

    'fr_para_hydroxylation': {
        'type': 'int',
        'description': 'Number of para-hydroxylation sites',
        'tags': ['2D']
    },

    'fr_phenol': {
        'type': 'int',
        'description': 'Number of phenols',
        'tags': ['2D']
    },

    'fr_phenol_noOrthoHbond': {
        'type': 'int',
        'description': 'Number of phenolic OH excluding ortho intramolecular Hbond substituents',
        'tags': ['2D']
    },

    'fr_phos_acid': {
        'type': 'int',
        'description': 'Number of phosphoric acid groups',
        'tags': ['2D']
    },

    'fr_phos_ester': {
        'type': 'int',
        'description': 'Number of phosphoric ester groups',
        'tags': ['2D']
    },

    'fr_piperdine': {
        'type': 'int',
        'description': 'Number of piperdine rings',
        'tags': ['2D']
    },

    'fr_piperzine': {
        'type': 'int',
        'description': 'Number of piperzine rings',
        'tags': ['2D']
    },

    'fr_priamide': {
        'type': 'int',
        'description': 'Number of primary amides',
        'tags': ['2D']
    },

    'fr_prisulfonamd': {
        'type': 'int',
        'description': 'Number of primary sulfonamides',
        'tags': ['2D']
    },

    'fr_pyridine': {
        'type': 'int',
        'description': 'Number of pyridine rings',
        'tags': ['2D']
    },

    'fr_quatN': {
        'type': 'int',
        'description': 'Number of quarternary nitrogens',
        'tags': ['2D']
    },

    'fr_sulfide': {
        'type': 'int',
        'description': 'Number of thioether',
        'tags': ['2D']
    },

    'fr_sulfonamd': {
        'type': 'int',
        'description': 'Number of sulfonamides',
        'tags': ['2D']
    },

    'fr_sulfone': {
        'type': 'int',
        'description': 'Number of sulfone groups',
        'tags': ['2D']
    },

    'fr_term_acetylene': {
        'type': 'int',
        'description': 'Number of terminal acetylenes',
        'tags': ['2D']
    },

    'fr_tetrazole': {
        'type': 'int',
        'description': 'Number of tetrazole rings',
        'tags': ['2D']
    },

    'fr_thiazole': {
        'type': 'int',
        'description': 'Number of tetrazole rings',
        'tags': ['2D']
    },

    'fr_thiocyan': {
        'type': 'int',
        'description': 'Number of thiocyanates',
        'tags': ['2D']
    },

    'fr_thiophene': {
        'type': 'int',
        'description': 'Number of thiophene rings',
        'tags': ['2D']
    },

    'fr_unbrch_alkane': {
        'type': 'int',
        'description': 'Number of unbranched alkanes of at least 4 members (excludes halogenated alkanes)',
        'tags': ['2D']
    },

    'fr_urea': {
        'type': 'int',
        'description': 'Number of urea groups',
        'tags': ['2D']
    },

    'NumHDonors': {
        'type': 'int',
        'description': 'The number of Hydrogen Bond Donors',
        'tags': ['2D']
    },

    'NumHAcceptors': {
        'type': 'int',
        'description': 'The number of Hydrogen Bond Acceptors',
        'tags': ['2D']
    },

    'NumHeteroatoms': {
        'type': 'int',
        'description': 'The number of Heteroatoms',
        'tags': ['2D']
    },

    'NumRotatableBonds': {
        'type': 'int',
        'description': 'The number of Rotatable Bonds',
        'tags': ['2D']
    },

    'NOCount': {
        'type': 'int',
        'description': 'The number of Nitrogens and Oxygens',
        'tags': ['2D']
    },

    'NHOHCount': {
        'type': 'int',
        'description': 'The number of NHs or OHs',
        'tags': ['2D']
    },

    'HeavyAtomCount': {
        'type': 'int',
        'description': 'The number of heavy atoms a molecule',
        'tags': ['2D']
    },

    'FractionCSP3': {
        'type': 'double',
        'description': 'The fraction of C atoms that are SP3 hybridized',
        'tags': ['2D']
    },

    'NumAliphaticCarbocycles': {
        'type': 'int',
        'description': 'The number of aliphatic (containing at least one non-aromatic bond) carbocycles for a molecule',
        'tags': ['2D']
    },

    'NumAliphaticHeterocycles': {
        'type': 'int',
        'description': 'The number of aliphatic (containing at least one non-aromatic bond) heterocycles for a molecule',
        'tags': ['2D']
    },

    'NumAliphaticRings': {
        'type': 'int',
        'description': 'The number of aliphatic (containing at least one non-aromatic bond) rings for a molecule',
        'tags': ['2D']
    },

    'NumAromaticCarbocycles': {
        'type': 'int',
        'description': 'The number of aromatic carbocycles for a molecule',
        'tags': ['2D']
    },

    'NumAromaticHeterocycles': {
        'type': 'int',
        'description': 'The number of aromatic heterocycles for a molecule',
        'tags': ['2D']
    },

    'NumAromaticRings': {
        'type': 'int',
        'description': 'The number of aromatic rings for a molecule',
        'tags': ['2D']
    },

    'NumSaturatedCarbocycles': {
        'type': 'int',
        'description': 'The number of saturated carbocycles for a molecule',
        'tags': ['2D']
    },

    'NumSaturatedHeterocycles': {
        'type': 'int',
        'description': 'The number of saturated heterocycles for a molecule',
        'tags': ['2D']
    },

    'NumSaturatedRings': {
        'type': 'int',
        'description': 'The number of saturated rings for a molecule',
        'tags': ['2D']
    },

    'RingCount': {
        'type': 'int',
        'description': 'Ring count',
        'tags': ['2D']
    },

    'PMI1': {
        'type': 'double',
        'description': 'First (smallest) principal moment of inertia',
        'tags': ['3D']
    },

    'PMI2': {
        'type': 'double',
        'description': 'Second principal moment of inertia',
        'tags': ['3D']
    },

    'PMI3': {
        'type': 'double',
        'description': 'Third (largest) principal moment of inertia',
        'tags': ['3D']
    },

    'NPR1': {
        'type': 'double',
        'description': 'Normalized principal moments ratio 1 (=I1/I3)',
        'tags': ['3D']
    },

    'NPR2': {
        'type': 'double',
        'description': 'Normalized principal moments ratio 2 (=I2/I3)',
        'tags': ['3D']
    },

    'RadiusOfGyration': {
        'type': 'double',
        'description': 'Radius of gyration',
        'tags': ['3D']
    },

    'InertialShapeFactor': {
        'type': 'double',
        'description': 'Inertial shape factor',
        'tags': ['3D']
    },

    'Eccentricity': {
        'type': 'double',
        'description': 'Molecular eccentricity',
        'tags': ['3D']
    },

    'Asphericity': {
        'type': 'double',
        'description': 'Molecular asphericity',
        'tags': ['3D']
    },

    'SpherocityIndex': {
        'type': 'double',
        'description': 'Molecular spherocity index',
        'tags': ['3D']
    }
}

def np_none(shape):
  return np.full(shape, None, dtype=object)

def get_3d_descriptors_names():
  desc_3d = []
  for d in descriptors_params:
    if '3D' in descriptors_params[d]['tags']:
      desc_3d.append(d)
  return desc_3d

def _get_descriptors_list():
  descs = Descriptors._descList[:]
  descs_extra = [
    ("PMI1", Chem.Descriptors3D.PMI1),
    ("PMI2", Chem.Descriptors3D.PMI2),
    ("PMI3", Chem.Descriptors3D.PMI3),
    ("NPR1", Chem.Descriptors3D.NPR1),
    ("NPR2", Chem.Descriptors3D.NPR2),
    ("RadiusOfGyration", Chem.Descriptors3D.RadiusOfGyration),
    ("InertialShapeFactor", Chem.Descriptors3D.InertialShapeFactor),
    ("Eccentricity", Chem.Descriptors3D.Eccentricity),
    ("Asphericity", Chem.Descriptors3D.Asphericity),
    ("SpherocityIndex", Chem.Descriptors3D.SpherocityIndex),
   ]
  for de in descs_extra:
    descs.append(de)
  return descs

def _get_descriptors_funcs(descriptors):
  descriptors_funcs = {}
  descs = _get_descriptors_list()
  for desc in descs:
    if desc[0] in descriptors:
      descriptors_funcs[desc[0]] = desc[1]
  return descriptors_funcs

def _add_3d_coordinates(mol):
  mol = Chem.AddHs(mol)
  AllChem.EmbedMolecule(mol, AllChem.ETKDG())
  mol = Chem.RemoveHs(mol)
  return mol

def get_descriptor_type(name):
    return descriptors_params[name]['type'] if name in descriptors_params else 'string'
def get_descriptor_tags(name):
    return descriptors_params[name]['tags'] if name in descriptors_params else []

def set_type(value, value_type):
    return {'type': value_type, 'value': value.tolist()}

def get_descriptor_description(name):
    return descriptors_params[name]['description'] if name in descriptors_params else ''

def get_module_description(name):
    return descriptors_modules[name]['description'] if name in descriptors_modules else ''

def get_module_name(name):
    return descriptors_modules[name]['name'] \
        if ((name in descriptors_modules) and ('name' in descriptors_modules[name])) else name

def get_descriptors_tree():
    """
    Gets all available Molecule descriptors in tree view (module -> descriptor).

    :return: List of available Molecule descriptors.
    """
    global _descriptors_tree
    if _descriptors_tree is None:
        _descriptors_tree = {}
        for desc in _get_descriptors_list():
            func = desc[1]
            name = desc[0]
            module_name = func.__module__.replace('rdkit.Chem.', '')
            version = func.__dict__['version'] if ('version' in func.__dict__) else 'unknown'
            descriptor = {
                'name': name,
                'description': get_descriptor_description(name),
                'type': get_descriptor_type(name),
                'tags': get_descriptor_tags(name)
            }
            if module_name in _descriptors_tree:
                _descriptors_tree[module_name]['descriptors'].append(descriptor)
            else:
                _descriptors_tree[module_name] = {
                    'version': version,
                    'name': get_module_name(module_name),
                    'description': get_module_description(module_name),
                    'descriptors': [descriptor]}

    return _descriptors_tree

_descriptors = list(descriptors)
tree = get_descriptors_tree()

for group in tree:
  if group in _descriptors:
    descriptors.remove(group)
    for descriptor in tree[group]['descriptors']:
      descriptors.append(descriptor['name'])
length = len(smiles)
values = []
for _ in descriptors:
  values.append(np_none(length))
descriptors_3d = get_3d_descriptors_names()
descriptors_funcs = _get_descriptors_funcs(_descriptors)
for n in range(0, length):
  mol = Chem.MolFromMolBlock(smiles[n], sanitize = True) if ("M  END" in smiles[n]) else Chem.MolFromSmiles(smiles[n], sanitize = True)

  if mol is None:
    continue
  try:
    for d in descriptors:
      if d in descriptors_3d:
        mol = _add_3d_coordinates(mol)
        break
    for d in range(0, len(descriptors)):
      value = descriptors_funcs[descriptors[d]](mol)
      if not np.isnan(value) and not np.isinf(value):
        values[d][n] = value if get_descriptor_type(descriptors[d]) != 'string' else str(value)
  except:
    continue
result = {}
for d in range(0, len(descriptors)):
  result[descriptors[d]] = values[d]

out = pd.DataFrame(result)
