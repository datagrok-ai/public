---
title: Info panes for chemical data
---

The [info panes](../../../../navigation/panels/info-panels.md) are displayed
automatically in the **Context Panel** when you select a molecule or mixture in a dataset.

| Panel <div style={{ width:130 }}></div> | Description <div style={{ width:700 }}></div> |
|-----|-----------|
| Gasteiger partial charges | [RDKit-based script](https://dev.datagrok.ai/script/7acf813d-4f65-51f2-bc3f-503cde26c460) |  
| Descriptors | The following descriptors are shown by default:<br/><li>FractionCSP3</li><li>HeavyAtomCount</li><li>NHOHCount</li><br/> To calculate other descriptors, click the **SELECT** button and choose from the list. The info pane updates accordingly |
| Properties | Shows molecular properties, such as formula, molecular weight, LogP, and others |
| Retrosynthesis | Shows the most efficient synthetic pathways and commercially available starting materials for your target molecules (based on [AiZynthFinder](https://github.com/MolecularAI/aizynthfinder)) |
| Databases | Depending on the plugins installed, allows to search by substructure or similarity in databases like ChEMBL, Chemspace, DrugBank, PubChem, as well as the Enamine store catalog |
| CDD Vault | Shows vault data for the current molecule |
| SureChEMBL | Shows patent information for target molecule based on similarity or substructure search |
| AutoDock | Displays docking results for molecules that have undergone AutoDock analysis |
| Admetica | [Calculates ADMET](https://github.com/datagrok-ai/public/tree/master/packages/Admetica). In addition, the **Summary** info pane visualizes ADMET in a pie chart |
| DiffDock | Provides an interactive interface for running molecular docking using NVIDIA's DiffDock model |
| Drug likeness | [Calculates drug likeness](drug-likeness.md) and displays the results |
| Structural alerts | [Calculates and displays structural alerts](structural-alerts.md) |
| Toxicity | [Calculates toxicity](toxicity-risks.md) and displays the results |
| Identifiers | Calculates identifiers like SMILES, InChi, ChEMBL ID, etc. |
| 3D Structure | Shows an interactive 3D view of the molecule |
| 2D Structure | Shows a 2D view of the molecule |
| Mixture | Available for mixtures. Shows all mixture components and their properties (structure, name, relation, etc.) in a tabular view. The table can be added to the workspace |
| Mixture Tree | Available for mixtures. Shows the mixture structure as a hierarchical tree with the mixfile version and all nested components |

<br/>

:::note developers

You can [create custom info panes](../../../../../develop/how-to/ui/add-info-panel.md).

:::
