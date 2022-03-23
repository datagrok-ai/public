<!-- TITLE: Tests: Chem info panels -->
<!-- SUBTITLE: -->

# Tests: Chem info panels

Platform supports special info-panels for chemical data that are displayed
on [Property Panel](../../overview/navigation.md#properties) for specific smiles value.

## Testing scenario

1. Open *smiles.csv* file

1. Click on first cell in "smiles" column

* [Property Panel](../../overview/navigation.md#properties) switches to display the info panels as tabs

1. Expand "Structure" tab on [Property Panel](../../overview/navigation.md#properties)

* Tab shows structure of the molecule

1. Expand "Properties" tab on [Property Panel](../../overview/navigation.md#properties)

* Tab shows parameters "Formula", "MW", "acceptorCount", "donorCount", "c", "logS", "
  polarSurfaceArea", "rotatableBondCount", "stereoCenterCount" and "name"

1. Expand "SDF" tab on [Property Panel](../../overview/navigation.md#properties)

* Tab shows .sdf file contents

1. Expand "3D" tab on [Property Panel](../../overview/navigation.md#properties)

* Tab shows molecular structure in 3D space
* Available rotation of the molecule in all axes
* Zoom in and out are available

1. Expand "Toxicity" tab on [Property Panel](../../overview/navigation.md#properties)

* Tab shows molecule toxicity by four parameters: "Mutagenicity", "Tumorigenicity", "Irritating effects" and "
  Reproductive effects"
* Degree of toxicity is color coded (none - green, high - dark red)

1. Expand "Drug Likeness" tab on [Property Panel](../../overview/navigation.md#properties)

* Tab shows similar molecules with their similarity score

1. Expand "Identifiers" tab on [Property Panel](../../overview/navigation.md#properties)

* Tab shows list of ChEMBL identifier for molecule

1. Expand "ChEMBL similarity" tab on [Property Panel](../../overview/navigation.md#properties)

* Tab shows list of similar molecules by ChEMBL

1. Expand "ChEMBL search" tab on [Property Panel](../../overview/navigation.md#properties)

* Tab shows list of molecules containing selected structure in ChEMBL

1. Expand "Structural Alerts" tab on [Property Panel](../../overview/navigation.md#properties)

* Tab shows list of structural alerts of current molecule
* Substructure alerts is highlighted in red

1. Expand "PubChem" tab on [Property Panel](../../overview/navigation.md#properties)

* Panel has tree structure of information about molecule from PubChema
* All sub-tabs open and contain relevant information

1. Expand "Gasteiger Partial Charges" tab on [Property Panel](../../overview/navigation.md#properties)

* Tab shows Gasteiger partial charges visualization based on RDKit
* Info panel is result of execution of corresponding Python script "Gasteiger Partial Charges"

1. Click on "Edit parameters" icon near "Gasteiger Partial Charges" tab

* "Contours" field appeared inside the tab

1. Change "Counters value to 50 and click "Apply" button

* Displayed Gasteiger partial charges visualization inside tab corresponds entered value of counters

1. Expand "Solubility prediction" tab on [Property Panel](../../overview/navigation.md#properties)

* Tab shows shows calculated value of molecule solubility prediction
* Info panel is result of execution of corresponding Grok script "Solubility predictions"
* Prediction calculated based on chem descriptors
* Predictive model is used for prediction

1. Go to other cells in "smiles" column. (you can use the "↓" key)

* Info panels display calculated data depending on cell values that is in focus
* If calculation takes time, corresponding animation is displayed in panel body

1. Click on "Settings" (gear icon) near "Properties" info-panel and select "Do not show"

* Properties" info-panel no longer displayed on [Property Panel](../../overview/navigation.md#properties)
* Corresponding message is shown on balun

1. Open **Tools | Settengs | Info Panels**

* Dialog tab shows disabled info panels ("Properties" info-panel)

1. Click on "UNBLOCK" near "Properties" info-panel

* Now the "Properties" info-panel displayed again
* Corresponding message is shown on balun

1. Click on "Help" (icon "?") near "Structural Alerts"

* Context help switched to page "Structural Alerts" which shows wiki for this panel
