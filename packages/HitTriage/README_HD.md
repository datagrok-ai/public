# HitDesign

HitDesign streamlines the process of designing new molecules and collaborating on molecular ideas.
You can sketch molecules,
calculate or predict chemical and biological properties, enrich data with information
from the proprietary databases, filter based on these properties or substructure search, select hits, 
move hits to different stages, save campaigns, and share your work with others.

The application is built around templates and campaigns, providing a structured yet flexible approach to
hit design process. In addition to the built-in functions for data ingestion, property calculation, and
data submission, you can define functions specific to your company or use case, allowing for the
seamless integration. 

## Templates

Templates in HitDesign contain essential configurations for conducting a campaign. Template configuration includes:

- **Name** : Identifies the template.

- **Campaign Prefix** : A code used as a prefix for campaign names (e.g., TMP-1, TMP-2).

- **Additional fields** : Configure additional fields for the template, which will be prompted for input during campaign creation. These fields include name, type, and whether they are required or not. For example, additional field for a campaign can be a target protein name, Head scientist name, deadlile, etc.

- **Stages** : Define stages for the campaign. Tiles view provides a versatile way to organize molecules in the campaign. Users can drag and drop molecules between stages. For example, stages can be used to organize molecules by their readiness for synthesis.

- **Compute functions** : HitDesign aggregates compute functions tagged with `HitTriageFunction` from Datagrok packages. Users can select from these functions to perform calculations (e.g., mass, solubility, mutagenicity, partial charges, toxicity risks, etc.) on the dataset.
Every time user changes given molecule or adds new molecule to a dataframe, compute functions are executed automatically for that row.

- **Submit function** : Define custom submit functions (tagged with `HitTriageSubmitFunction`) to further process or save the filtered and computed dataset. This could include saving to a private database or additional calculations.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/template-HD.png?raw=true)

## Campaigns

Campaigns are built based on templates and encompass the actual hit design process.

- **Creation** : Select a template, fill out additional information and data source, and start the campaign.

- **Collaboration** : Campaigns are automatically shared. Users can copy and share URLs for seamless collaboration.

- **Functionality**: Once a campaign starts, you can add extra calculated columns, add new rows to dataframe and sketch molecules, apply changes, fileter, save or submit the campaign.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/campaign-HD.png?raw=true)

Campaigns can be groupped either by template key, status or be ungroupped. The groupping setting can be changed by clicking 'Group By' icon next to 'Comntinue Campaign' header.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/HD-modify-groupping?raw=true)

## Getting started

Continue ongoing campaigns either directly by a link or by selecting it from the campaigns table.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/continue_campaign_HD.gif?raw=true)

Create a new template by clicking on the `New Template` button in the `Templates` dropdown.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/create-template-HD.gif?raw=true)

Start a new campaign by choosing a template and filling out the required information. Once a new campaign starts, an empty dataframe will be added, where you can add new molecules. Upon adding or changing molecules, all compute functions will be executed on that row and the results in coresponding columns will be updated.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/HD-start-campaign.gif?raw=true)

After the campaign starts, users can sketch new molecules, filter, modify or add viewers to the campaign and then save them. Once saved, reloading the campaign will restore the saved state.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/HD-save-campaign.gif?raw=true)

Hit design campaign consists of two views, a main design view and a tiles view. You can access the tiles view from by clicking the 'Progrss Tracker' button on the ribbon pannel. Tiles view provides a versatile way to organize molecules in the campaign. Users can drag and drop molecules between stages, which were defined in the template.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/HD-tiles.gif?raw=true)

As mentioned previously, "Stages" or tiles are defined in the template, but users can also modify them after starting the campaign. To do so, open the Progress Tracker and click the "Modify Stages" button. You can add, remove or rename stages. After clicking OK, the changes will be saved and reflected in the Progress Tracker view. If you remove a stage which had some molecules in it, those molecules will be automatically moved to the first stage.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/HD-tiles-modify-stages.gif?raw=true)

Campaign status can be changed through 'Submit' menu. You can set status to anything. When you start typing, you will be suggested other statuses that are already in use.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/HD-modify-status.gif?raw=true)

## Adding custom compute and submit functions

HitDesign allows users to define custom compute and submit functions, and these functions can be written in any Datagrok package that is installed in the environment. 

### Compute Functions

Compute functions are used to calculate molecular properties. For example, mass, solubility, mutagenicity, partial charges, toxicity risks, etc. By default, Hit design will include compute functions from `Chem` package, which are molecular descriptors, Structural alerts, Toxicity risks and Chemical properties. Users can add additional compute functions by tagging them with `HitDesignFunction` tag and writing them in normal datagrok style. The First two inputs of these functions should be `Dataframe` `table` and `Column` `molecule`, and rest can be any other input. Function should perform a certain task, modify the dataframe in desired way and return the modified dataframe. For example, we can create a function that retrieves the `Chembl` mol registration number by smiles string:

```
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

This function will go through every molecule in the dataframe, convert them to canonical smiles and call the query from Chembl database, that will retrieve the molregno number. The result will be added as a new column to the dataframe. If this function is defined in the `Chembl` package, after building and deploying it to stand, it will be automatically added to the compute functions list in HitDesign.

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

Similarly, queries with same `HitTriageFunction` tag will be added to the compute functions list. The query needs to have at least one input, first of which must be `list<string>`, representing the list of molecules. The query must return a dataframe, which should contain column `molecules` in order to join result with initial dataframe. `molecules` column will be used as key for joining tables. For example, you can create a query that looks for the molecule in Chembl database and returns the molregno number:

```sql
--name: ChemblMolregNoBySmiles
--friendlyName: Chembl Molregno by smiles
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

### Submit Functions

Submit functions are used to save or submit the filtered and computed dataset. This could include saving to a private database or additional calculations. Submit functions are defined in the same way as compute functions, but they are tagged with `HitTriageSubmitFunction` tag. The function should accept only two inputs, `Dataframe` `df` and `String` `molecules`, which are the resulting dataframe and name of molecules column respectively. For example, we can create a function that saves the filtered and computed dataset to a database:

```typescript
//name: Sample File Submit
//tags: HitTriageSubmitFunction
//input: dataframe df [dataframe]
//input: string molecules [molecules column name]
export async function demoFileSubmit(df: DG.DataFrame, molecules: string): Promise<void> {
    const smiles = df.getCol(molecules).toList();
    myCustomDatabase.add(smiles) // template function which could for example save the smiles to a database
}
```

## Permission Management

Hit Design allows users to manage permissions for campaigns, granting `view` or `edit` access to specific `users` or `groups`.

By default, any new campaign will be shared with all users. To manage who can view or edit the campaign, click on the `share` icon located in the campaigns table, or open the campaign and click on the `share` button in the ribbon menu. If all users are removed from `edit` permissions, only the creator of the campaign will be able to edit it. Users who do not have edit permission but can view the campaign, will be able to see the campaign but not make any changes or delete it.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/HD-permissions.gif?raw=true)

Similarly, users can manage view permissions, which will restrict access to viewing the campaign or accessing it directly via link. If all users are removed from `view` permissions, only the creator of the campaign will be able to view it.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/HD-view-permissions.gif?raw=true)