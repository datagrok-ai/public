# HitTriage

HitTriage is a dynamic application designed for chemists and biologists to streamline the process of selecting molecular hits. It facilitates collaboration by allowing users to load datasets, calculate molecular properties, filter based on these properties or substructure search, select hits, save campaigns, and share their work with others. The application is built around templates and campaigns, providing a structured yet flexible approach to hit triage process.

## Templates

Templates in HitTriage contain essential configurations for conducting a campaign. template configuration includes:

- **Name** : Identifies the template.

- **Campaign Prefix** : A code used as a prefix for campaign names (e.g., TMP-1, TMP-2).

- **Data Ingestion Settings** : Controls how molecular data is ingested into the template. Users can choose between file upload or a query option. For the query, HitTriage searches for functions in any Datagrok package tagged with `HitTriageDataSource`. For example:

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
    The application will detect that the function requeires an input parameter and will prompt the user to provide it in campaigns form. The function must return a dataframe with a column containing molecules

- **Additional fields** : Users can configure additional fields for the template, which will be prompted for input during campaign creation. These fields include name, type, and whether they are required or not. For example, additional field for a campaign can be a target protein name, Head scientist name, deadlile, etc.

- **Compute functions** : HitTriage aggregates compute functions tagged with `HitTriageFunction` from Datagrok packages. Users can select from these functions to perform calculations (e.g., mass, solubility, mutagenicity, partial charges, toxicity risks, etc.) on the dataset.

- **Submit function** : Users can define custom submit functions (tagged with `HitTriageSubmitFunction`) to further process or save the filtered and computed dataset. This could include saving to a private database or additional calculations.

![template](https://datagrok.ai/help/uploads/hittriage/files/images/template.png)

## Campaigns

Campaigns are built based on templates and encompass the actual hit triage process.

- **Creation** : Users select a template, fill out additional information and data source, and start the campaign.

- **Collaboration** : Campaigns are automatically shared. Users can copy and share URLs for seamless collaboration.

- **Functionality**: Once a campaign starts, you can add extra calculated columns, apply changes, fileter, save or submit the campaign.

![template](https://datagrok.ai/help/uploads/hittriage/files/images/campaign.png)

## Getting started

Users can continue ongoing campaigns either directly by a link or by selecting it from the campaigns table.

![template](https://datagrok.ai/help/uploads/hittriage/files/images/HT_Continue_campaign.gif)

Users can create a new template by clicking on the `New Template` button in the `Templates` dropdown.

![template](https://datagrok.ai/help/uploads/hittriage/files/images/HT_create_template.gif)

Users can start a new campaign by choosing a template and filling out the required information. 

![template](https://datagrok.ai/help/uploads/hittriage/files/images/HT_create_campaign.gif)

After the campaign starts, users can filter, modify or add viewers to the campaign and then save them. once saved, reloading the campaign will restore the saved state.

![template](https://datagrok.ai/help/uploads/hittriage/files/images/HT_save_campaign.gif)