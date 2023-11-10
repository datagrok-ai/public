# HitDesign

HitDesign is a dynamic application designed for chemists and biologists to streamline the process of designing and selecting molecular hits. It facilitates collaboration by allowing users to sketch and add molecules, calculate molecular properties, filter based on these properties or substructure search, select hits, move hits to different stages, save campaigns, and share their work with others. The application is built around templates and campaigns, providing a structured yet flexible approach to hit design process.

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

## Getting started

Continue ongoing campaigns either directly by a link or by selecting it from the campaigns table.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/continue_campaign_HD.gif?raw=true)

Create a new template by clicking on the `New Template` button in the `Templates` dropdown.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/create-template-HD.gif?raw=true)

Start a new campaign by choosing a template and filling out the required information. Once a new campaign starts, an empty dataframe will be added, where you can add new molecules. Upon adding or changing molecules, all compute functions will be executed on that row and the results in coresponding columns will be updated.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/HD-start-campaign.gif?raw=true)

After the campaign starts, users can sketch new molecules, filter, modify or add viewers to the campaign and then save them. Once saved, reloading the campaign will restore the saved state.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/HD-save-campaign.gif?raw=true)

Hit design campaign consists of two views, a main design view and a tiles view. You can access the tiles view from the views list. Tiles view provides a versatile way to organize molecules in the campaign. Users can drag and drop molecules between stages, which were defined in the template.

![hitDesignReadmeImg](https://github.com/datagrok-ai/public/blob/master/help/uploads/hittriage/HD-tiles.gif?raw=true)

## Adding custom compute and submit functions

HitDesign allows users to define custom compute and submit functions, and these functions can be written in any Datagrok package that is installed in the environment. 

** Compute functions **

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

** Submit functions **

Submit functions are used to save or submit the filtered and computed dataset. This could include saving to a private database or additional calculations. Submit functions are defined in the same way as compute functions, but they are tagged with `HitTriageSubmitFunction` tag. The function should accept only two inputs, `Dataframe` `df` and `String` `molecules`, which are the resulting dataframe and name of molecules column respectively. For example, we can create a function that saves the filtered and computed dataset to a database:

```
//name: Sample File Submit
//tags: HitTriageSubmitFunction
//input: dataframe df [dataframe]
//input: string molecules [molecules column name]
export async function demoFileSubmit(df: DG.DataFrame, molecules: string): Promise<void> {
    const smiles = df.getCol(molecules).toList();
    myCustomDatabase.add(smiles) // template function which could for example save the smiles to a database
}
```
