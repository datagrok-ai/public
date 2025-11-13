---
title: Info panes for chemical data
---

The [info panes](../../../../navigation/panels/info-panels.md) are displayed
automatically in the **Context Panel** when you select a molecule in a dataset.

|Panel|Description|
|-----|-----------|
|Gasteiger partial charges | [RDKit-based script](https://dev.datagrok.ai/script/7acf813d-4f65-51f2-bc3f-503cde26c460) |  
|Descriptors |The following descriptors are shown by default: FractionCSP3, HeavyAtomCount, NHOHCount.<br/> To calculate other descriptors, click the **SELECT** button and choose from the list. The info pane updates accordingly|
|Properties | Shows molecular properties, such as formula, molecular weight, LogP, and others|
| Databases | Depending on the plugins installed, allows to search by substructure or similarity in databases like ChEMBL, Chemspace, DrugBank, PubChem, as well as the Enamine store catalog |
|Drug likeness | [Calculates drug likeness](drug-likeness.md) and displays the results  |
|Structural alerts | [Calculates and displays structural alerts](structural-alerts.md)  |
|Toxicity | [Calculates toxicity](toxicity-risks.md) and displays the results  |
|ADME/Tox | [Calculates ADMET](https://github.com/datagrok-ai/public/tree/master/packages/Admetica). In addition, the **Summary** info pane visualizes ADMET in a pie chart |
|Identifiers | Calculates identifiers like SMILES, InChi, ChEMBL ID, etc.|
|3D Structure| Shows an interactive 3D view of the molecule|
|2D Structure| Shows a 2D view of the molecule|
| Retrosynthesis|Shows the most efficient synthetic pathways and commercially available starting materials for your target molecules (based on [AiZynthFinder](https://github.com/MolecularAI/aizynthfinder))|
|CDD Vault|Shows vault data for the current molecule|
|SureChEMBL|Shows patent information for target molecule based on similarity or substructure search|

<br/>

:::note developers

You can [create custom info panes](../../../../../develop/how-to/ui/add-info-panel.md).

:::
