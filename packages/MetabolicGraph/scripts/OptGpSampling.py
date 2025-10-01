#name: OptGpSampling
#description: Samples the metabolic map and returns binned values of samples
#language: python
#environment: channels: [Conda-forge], dependencies: [python=3.12, {pip: [cobra]}]
#input: string cobraModel
#input: int nSamples = 1000 {nullable: true}
#output: dataframe res
from cobra.core import Metabolite, Model, Reaction

from cobra.sampling import ACHRSampler, OptGPSampler, sample

from cobra.io.dict import model_from_dict

import json

# open local file

jsonMap = json.loads(cobraModel)

model = model_from_dict(jsonMap)

optgp = OptGPSampler(model, 1, seed=42, n_samples = nSamples)

optgp_samples = optgp.sample(nSamples)

import pandas as pd
import numpy as np

def bin_dataframe(df: pd.DataFrame, bins=20):
    binned_data = {}
    globMin = df.min().min()
    glomalMax = df.max().max()
    # store averages of original values
    averages = df.mean().to_frame().T
    for col in df.select_dtypes(include=[np.number]):  # Process only numerical columns
        bin_edges = np.linspace(globMin, glomalMax, bins + 1)
        bin_labels = [f"[{bin_edges[i]:.2f}, {bin_edges[i+1]:.2f})" for i in range(len(bin_edges)-1)]
        df[f"{col}_bin"] = pd.cut(df[col], bins=bin_edges, labels=bin_labels, include_lowest=True)
        
        binned_counts = df[f"{col}_bin"].value_counts().sort_index()
        binned_data[col] = binned_counts

    # Create a final dataframe
    binned_df = pd.DataFrame(binned_data)
    binned_df.index.name = "Bin Range"
    res_df = binned_df.reset_index()
    # Add another row at the end with average values of original values
    
    averages.index = ["Average"]
    binned_df = pd.concat([binned_df, averages], ignore_index=False)
    binned_df.index.name = "Bin Range"
    binned_df = binned_df.fillna(0)  # Fill NaN values with 0 for bins with no counts
    return binned_df.reset_index()

binned_df = bin_dataframe(optgp_samples, bins=30)
# add another column with ammount of actual samples, i.e. number of rows in the original dataframe
# binned_df['Samples'] = len(optgp_samples)

res = binned_df