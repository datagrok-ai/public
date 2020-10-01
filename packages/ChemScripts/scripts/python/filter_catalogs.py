#name: Filter by Catalogs
#description: Finds undesireable molecules based on various criteria
#help-url: https://datagrok.ai/help/domains/chem/functions/filter-catalogs
#language: python
#sample: chem/smiles.csv
#tags: demo, chem, rdkit
#input: dataframe data [Input data table]
#input: column smiles {type:categorical, semType: Molecule} [Molecules, in SMILES format]
#input: string catalog = BRENK {choices: ["BRENK", "NIH", "PAINS_A", "PAINS_B", "PAINS_C", "ZINC"]}
#output: dataframe filter {action:join(data)} [Column with filter matches]

import numpy as np
from rdkit import Chem
from rdkit.Chem import FilterCatalog

smiles = data[smiles]

params = FilterCatalog.FilterCatalogParams()
catalogs = {
    "BRENK": FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK,
    "NIH": FilterCatalog.FilterCatalogParams.FilterCatalogs.NIH,
    "PAINS_A": FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A,
    "PAINS_B": FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B,
    "PAINS_C": FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C,
    "ZINC": FilterCatalog.FilterCatalogParams.FilterCatalogs.ZINC
}
params.AddCatalog(catalogs[catalog])
filterCatalog = FilterCatalog.FilterCatalog(params)

length = len(smiles)
filter = np.zeros(length, dtype=bool)
for n in range(0, length):
    mol = Chem.MolFromSmiles(smiles[n])
    if mol is None:
        continue
    filter[n] = filterCatalog.HasMatch(mol)

# Convert to Pandas DataFrame
filter = pd.DataFrame(filter, columns=[catalog])
