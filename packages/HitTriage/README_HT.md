# Hit Triage

Hit Triage is a dynamic application designed for chemists and biologists to streamline the process of selecting molecular hits. It facilitates collaboration by allowing users to load datasets, calculate molecular properties, filter based on these properties or substructure search, select hits, save campaigns, and share their work with others. The application is built around templates and campaigns, providing a structured yet flexible approach to hit triage process.

## Templates

Templates in Hit Triage contain essential configurations for conducting a campaign. template configuration includes:

- **Name** : Identifies the template.

- **Campaign Prefix** : A code used as a prefix for campaign names (e.g., TMP-1, TMP-2).

- **Data Ingestion Settings** : Controls how molecular data is ingested into the template. Users can choose between file upload or a query option. For the query, Hit Triage searches for functions in any Datagrok package tagged with `HitTriageDataSource` and queries with same tag. For example, a package can export the following function:

```//name: Demo File Ingestion
//input: int numberOfMolecules [Molecules count]
//tags: HitTriageDataSource
//output: dataframe result
export async function demoFileIngest(numberOfMolecules: number): Promise<DG.DataFrame> {
const df = grok.data.demo.molecules(numberOfMolecules);
df.name = 'Variable Molecules number';
return df;
}
```
Or users can write a query in the query editor, save it and share with others. Example for loading molecules from Chembl database:

```--name: _someChemblStructure
    --friendlyName: Load Some Chembl structures
    --input: int numberOfMolecules = 1000
    --tags: HitTriageDataSource
    --connection: Chembl
    select
    canonical_smiles, molregno
    from
    compound_structures
    limit @numberOfMolecules
```

The application will detect that the function/query requeires an input parameter and will prompt the user to provide it in campaigns form. The function\query must return a dataframe with a column containing molecules.

- **Additional fields** : Users can configure additional fields for the template, which will be prompted for input during campaign creation. These fields include name, type, and whether they are required or not. For example, additional field for a campaign can be a target protein name, Head scientist name, deadlile, etc.

- **Compute functions**

Compute functions are used to calculate molecular properties. For example, mass, solubility, mutagenicity, partial charges, toxicity risks, etc. By default, Hit design will include compute functions from `Chem` package, which are molecular descriptors, Structural alerts, Toxicity risks and Chemical properties. Users can add additional compute functions by tagging them with `HitDesignFunction` tag and writing them in normal datagrok style. The First two inputs of these functions should be `Dataframe` `table` and `Column` `molecule`, and rest can be any other input. Function should perform a certain task, modify the dataframe in desired way and return the modified dataframe. For example, we can create a function that retrieves the `Chembl` mol registration number by smiles string:

```typescript
//name: Chembl molregno
//tags: HitTriageFunction
//input: dataframe table [Input data table] {caption: Table}
//input: column molecules {caption: Molecules; semType: Molecule}
//output: dataframe result
export async function chemblMolregno(table: DG.DataFrame, molecules: DG.Column): Promise<DG.DataFrame> {
  const name = table.columns.getUnusedName('CHEMBL molregno');
  table.columns.addNewInt(name);
  for (let i = 0; i < molecules.length; i++) {
    const smile = molecules.get(i);
    if (!smile) {
      table.set(name, i, null);
      continue;
    }
    const canonical = grok.chem.convert(smile, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles);
    const resDf: DG.DataFrame = await grok.data.query('Chembl:ChemblMolregNoBySmiles', {smiles: canonical});
    const res: number = resDf.getCol('molregno').toList()[0];
    table.set(name, i, res);
  }
  return table;
}
```

This function will go through every molecule in the dataframe, convert them to canonical smiles and call the query from Chembl database, that will retrieve the molregno number. The result will be added as a new column to the dataframe. If this function is defined in the `Chembl` package, after building and deploying it to stand, it will be automatically added to the compute functions list in Hit Design.

Datagrok scripts can also be used as compute functions. For example, you can create a js script that adds a new column to the dataframe. This script also needs to have `HitTriageFunction` tag and should accept `Dataframe` `table` and `Column` `molecules` as first two inputs:

```javascript
//name: Demo script HT
//description: Hello world script
//language: javascript
//input: dataframe df
//input: column col
//input: int a
//tags: HitTriageFunction
//output: dataframe res

df.columns.addNewInt('Some number col').init(() => a)
res = df

```

Or a python script that calculates the number of atoms in the molecule and multiplies it by a specified value. In case of python, you need to return the dataframe containing columns that you want to append:

```python
#name: HTPythonDemo
#description: Calculates number of atoms in mulecule in python and also multiplies it by specified value 'multiplier'
#language: python
#tags: HitTriageFunction
#input: dataframe table [Data table]
#input: column col {semType: Molecule}
#input: int multiplier
#output: dataframe result

from rdkit import Chem
import numpy as np
# in python, column is passed as column name and dataframes are in pandas format.
# first, get the column.
molecules = table[col]
length = len(molecules)
# create array of same length
resCol = np.full(length, None, dtype=object)
for n in range(0, length):
	if molecules[n] == "":
		continue
	try:
		mol = Chem.MolFromMolBlock(molecules[n], sanitize = True) if ("M  END" in molecules[n]) else Chem.MolFromSmiles(molecules[n], sanitize = True)
		if mol is None or mol.GetNumAtoms() == 0:
			continue
		resCol[n] = mol.GetNumAtoms() * multiplier
	except:
		continue
result = pd.DataFrame({'Number of Atoms * mult': resCol})
```

Similarly, queries with same `HitTriageFunction` tag will be added to the compute functions list. The query needs to have at least one input, first of which must be `list<string>`, representing the list of molecules. The query must return a dataframe, which should contain column `molecules` in order to join result with initial dataframe. `molecules` column will be used as key for joining tables. For example, we can create a query that looks for the molecule in Chembl database and returns the molregno number:

```sql
--name: ChemblMolregNoBySmilesDirect
--friendlyName: Chembl Molregno by smiles direct
--input: list<string> molecules
--tags: HitTriageFunction
--connection: Chembl
select molregno, molecules from compound_structures c
	INNER JOIN unnest(@molecules) molecules
    ON molecules.molecules
 = c.canonical_smiles
```

Or a query that calculates fraction of sp3 hybridized carbons in the molecule using RDKit SQL cartridge:

```sql
--name: SP3Fraction
--friendlyName: SP3 fraction of carbons
--input: list<string> molecules
--tags: HitTriageFunction
--connection: Chembl
select molecules, mol_fractioncsp3(Cast(molecules as mol))
from unnest(@molecules) as molecules
where is_valid_smiles(Cast(molecules as cstring))
```

- **Submit function** : Users can define custom submit functions (tagged with `HitTriageSubmitFunction`) to further process or save the filtered and computed dataset. This could include saving to a private database or additional calculations.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/template.png?raw=true)

## Campaigns

Campaigns are built based on templates and encompass the actual hit triage process.

- **Creation** : Users select a template, fill out additional information and data source, and start the campaign.

- **Collaboration** : Campaigns are automatically shared. Users can copy and share URLs for seamless collaboration.

- **Functionality**: Once a campaign starts, you can add extra calculated columns, apply changes, fileter, save or submit the campaign.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/campaign.png?raw=true)

## Getting started

Users can continue ongoing campaigns either directly by a link or by selecting it from the campaigns table.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/HT_Continue_campaign.gif?raw=true)

Users can create a new template by clicking on the `New Template` button in the `Templates` dropdown.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/HT_create_template.gif?raw=true)

Users can start a new campaign by choosing a template and filling out the required information. 

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/HT_create_campaign.gif?raw=true)

After the campaign starts, new calculated columns will be added. Users can filter, modify or add viewers to the campaign and then save them. Once saved, reloading the campaign will restore the saved state.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/HT_save_campaign.gif?raw=true)