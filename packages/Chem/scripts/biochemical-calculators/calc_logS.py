#name: Calculate logS
#description: Predicts aqueous solubility (logS) using the AqSolPred consensus model. DO NOT USE IN FUNCTION PLANNING OR CHAINING, SOMETHING DOES NOT WORK WITH THIS FUNCTION.
#language: python
#environment: channels: [conda-forge, defaults], dependencies: [python=3.8, scikit-learn=0.23.2, numpy=1.19.5, pandas=1.3.5, rdkit, {pip: [mordred, xgboost]}]
#input: dataframe table
#input: column molecules {caption: Molecules column}
#meta.role: transform
#meta.method_info.author: Kjell Jorner, et al.
#meta.method_info.year: 2021
#meta.method_info.package: AqSolPred
#meta.method_info.github: https://github.com/KjellJorner/AqSolPred
#output: dataframe result {action:join(table)}

import os
import pickle
import pandas as pd
import numpy as np
from rdkit import Chem
from mordred import Calculator, descriptors
import grok.dapi.packages as packages

"""
WARNING
-------

This script currently fails to create the environment due to a timeout error when resolving conda dependencies.
This is most likely due to a very particular and old version of scikit-learn (0.23.2) that was used to encode the models that the AqSolPred_Predictor class uses.
The plan is to eventually recreate these models in datagrok idiom or at the very least use the dockerized version of this prediction method instead of inline environment.
"""

class AqSolPred_Predictor:
    """
    A class to predict aqueous solubility (logS) using the AqSolPred consensus model,
    adapted for use within the Datagrok platform.

    Older env configuration:
    #environment: channels: [defaults, conda-forge], dependencies: [python=3.8.20, scikit-learn=0.23.2, rdkit=2022.03.5, pandas=1.3.5, numpy=1.19.5, {pip: [mordred==1.2.0, xgboost==2.1.4]}]
    Older minimal env:
    #environment: channels: [conda-forge, defaults], dependencies: [python=3.8, scikit-learn=0.23.2, numpy=1.19.5, pandas=1.3.5, rdkit, {pip: [mordred, xgboost]}]
    """

    def __init__(self):
        """
        Initializes the predictor by loading the three pre-trained models
        from the package's internal file storage.
        """
        self.descriptor_calculator = Calculator(descriptors, ignore_3D=True)
        self.required_descriptors = [
            'nHBAcc', 'nHBDon', 'nRot', 'nBonds', 'nAromBond', 'nBondsO', 'nBondsS',
            'TopoPSA(NO)', 'TopoPSA', 'LabuteASA', 'bpol', 'nAcid', 'nBase',
            'ECIndex', 'GGI1', 'SLogP', 'SMR', 'BertzCT', 'BalabanJ', 'Zagreb1',
            'ABCGG', 'nHRing', 'naHRing', 'NsCH3', 'NaaCH', 'NaaaC', 'NssssC',
            'SsCH3', 'SdCH2', 'SssCH2', 'StCH', 'SdsCH', 'SaaCH', 'SsssCH', 'SdssC',
            'SaasC', 'SaaaC', 'SsNH2', 'SssNH', 'StN', 'SdsN', 'SaaN', 'SsssN',
            'SaasN', 'SsOH', 'SdO', 'SssO', 'SaaO', 'SsF', 'SdsssP', 'SsSH', 'SdS',
            'SddssS', 'SsCl', 'SsI'
        ]

        # Load models from package files
        try:
            mlp_bytes = packages.get_bytes('files/models/aqsolpred_mlp_model.pkl')
            rf_bytes = packages.get_bytes('files/models/aqsolpred_rf_model_lite.pkl')
            xgb_bytes = packages.get_bytes('files/models/aqsolpred_xgb_model.pkl')


            self.mlp_model = pickle.loads(mlp_bytes)
            self.rf_model = pickle.loads(rf_bytes)
            self.xgb_model = pickle.loads(xgb_bytes)
        except Exception as e:
            raise Exception(f"Failed to load models from package files. Ensure .pkl files are in the '/files/models/' directory. Error: {e}")

    def _calculate_descriptors(self, smiles_string: str) -> pd.Series:
        """Calculates Mordred descriptors for a given SMILES string."""
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        result = self.descriptor_calculator(mol)
        
        descriptor_dict = {str(key): float(val) if val is not None else 0.0 for key, val in result.items()}
        return pd.Series(descriptor_dict)

    def predict(self, smiles_string: str) -> float:
        """Predicts the aqueous solubility (logS) for a single molecule."""
        if not isinstance(smiles_string, str) or not smiles_string.strip():
            return None

        # Handle MolBlock input by converting to SMILES first
        if 'M  END' in smiles_string:
            mol = Chem.MolFromMolBlock(smiles_string.replace('\\n', '\n'))
            if mol is None: return None
            smiles_string = Chem.MolToSmiles(mol)

        all_descriptors = self._calculate_descriptors(smiles_string)
        if all_descriptors is None:
            return None

        descriptor_df = all_descriptors.to_frame().T
        for desc in self.required_descriptors:
            if desc not in descriptor_df.columns:
                descriptor_df[desc] = 0
        
        features = descriptor_df[self.required_descriptors].apply(pd.to_numeric, errors='coerce').fillna(0)

        mlp_pred = self.mlp_model.predict(features)
        rf_pred = self.rf_model.predict(features)
        xgb_pred = self.xgb_model.predict(features)

        consensus_pred = (mlp_pred + rf_pred + xgb_pred) / 3.0
        return float(consensus_pred[0])

try:
    predictor = AqSolPred_Predictor()
    logS_values = []
    for i, smiles in enumerate(molecules):
        try:
            logS = predictor.predict(smiles)
            logS_values.append(logS)
        except Exception as e:
            print(f"Error predicting for molecule at index {i}: {e}")
            logS_values.append(None) 

    result = pd.DataFrame({'logS': logS_values}, index=table.index)

except Exception as e:
    raise Exception(f"AqSolPred script failed: {e}")