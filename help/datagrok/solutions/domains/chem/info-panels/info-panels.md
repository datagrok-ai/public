---
title: Info panes for chemical data
---

The [info panes](../../../../navigation/panels/info-panels.md) are displayed
automatically in the **Context Panel** when you select a molecule in a dataset.

|Panel|Description|
|-----|-----------|
|[Gasteiger partial charges](link) | [RDKit-based script](https://dev.datagrok.ai/script/7acf813d-4f65-51f2-bc3f-503cde26c460) |  
|Descriptors |The following descriptors are shown by default:<br/><li>FractionCSP3</li><li>HeavyAtomCount</li><li>NHOHCount</li><br/> To calculate other descriptors, click the **SELECT** button and choose from the list. The info pane updates accordingly|
|Properties | Shows molecular properties, such as formula, molecular weight, LogP, and others|
| Databases | Depending on the plugins installed, allows to search by substructure or similarity in databases like ChEMBL, Chemspace, DrugBank, PubChem, as well as the Enamine store catalog |
|Drug likeness | [Calculates drug likeness](drug-likeness.md) and displays the results  |
|Structural alerts | [Calculates and displays structural alerts](structural-alerts.md)  |
|Toxicity | [Calculates toxicity](toxicity-risks.md) and displays the results  |
|ADME/Tox | [Calculates ADMET](https://github.com/datagrok-ai/public/tree/master/packages/ADMETox). In addition, the **Summary** info pane visualizes ADMET in a pie chart |
|Identifiers | Calculates identifiers like SMILES, InChi, ChEMBL ID, etc.|
|3D Structure| Shows an interactive 3D view of the molecule|
|2D Structure| Shows a 2D view of the molecule|

<br/>

:::note developers

You can [create custom info panes](../../../../../develop/how-to/add-info-panel.md).

:::
