# Chem changelog

## 1.7.2 (2023-09-05)

### Bug Fixes

* Chem: Scaffold Tree: Checkbox shouldn't be set, when group is expanded
* Chem: Scaffold Tree: Fix the behaviour of allowGenerate property

## 1.7.1 (2023-08-31)

### Features

* GROK-13571: Chem | Ability to terminate substructure search if substructure has been changed

### Bug Fixes

* GROK-13327: Chem | Substructure Search: two identical panels open on the Filter Panel
* GROK-13791: Chem | Chemical space (using t-SNE) fails on smiles dataset
* [#2135](https://github.com/datagrok-ai/public/issues/2135):
  * The structure rendering is too small.
* [#2322](https://github.com/datagrok-ai/public/issues/2322): Properties panel is unexpectedly reset on changing viewer properties if there is a scaffold tree filter in filters panel 
* GROK-13848: Chem: Substructure search results flickering

## 1.7.0 (2023-08-09)

### Features

* GROK-13172: Chem | implement substructure search using JSSubstructLibrary

## 1.6.22 (2023-08-07)

### Bug Fixes

* GROK-13713: Chem | Incorrect molecule rendering

## 1.6.21 (2023-08-02)

*Dependency: datagarok-api >= 1.16.0*

### Features

* Calculate drug likeness, toxicity and alerts for whole table from widgets
* Color Coding for toxicity
* Scaffold Tree improvements:
  * [#2154](https://github.com/datagrok-ai/public/issues/2154): Scaffold Tree: harmonization.

### Bug Fixes

* GROK-13586: _chemFindSimilar fails with 'Cannot read properties of null (reading 'rows')'
* [#2135](https://github.com/datagrok-ai/public/issues/2135):
  * The counts and controls are partially hidden
  * When change the drawing of a scaffold, then selection (checkbox) of all other scaffolds gets reset (to unchecked)
* [#2139](https://github.com/datagrok-ai/public/issues/2139): Scaffold tree stops working after adding an invalid structure

## 1.6.20 (2023-07-21)

This release focuses on improvements and bug fixes.

*Dependency: datagarok-api >= 1.14.0*

### Features

* Set default values in all dialogs where appropriate.
* Unified the layout of the search panel results.
* Work with DBs APIs: Added links to sources, similarity search viewer style, reproduce for ChemicalSpace, Enamine.
* **Elemental Analysis**: the resulting column name now includes the specific column for which calculations were conducted.
* **R-Groups Analysis**: a new option to choose between searching for MCS exact atoms or MCS exact bonds.
* **Similarity Search**: a new **Follow Current Row** setting prevent recalculation when the current row is changed. 
* For proper handling of properties and rendering, we now check for smarts and molecular fragments separately.
* Ability to copy data from the **Descriptors** and **Properties** tabs on the **Context Pane**.
* Moved **Descriptors** and **Fingerprints** from the  **Context Pane**  to the **Top Menu** ( **Chem** > **Calculate**).
* Added **Substructure Search** to the **Top Menu** ( **Chem** > **Search**)
* Modified the tooltip for dialog and drag-n-drop to prevent it from overlapping with the data.
* Implemented RDKit rendering for Chembl, ChemblAPI, PubChem, and DrugBank databases if OCL is used currently.
* UI polishing: harmonized input field names, repositioned elements in dialogs, added tooltips, and organized the top menu items into groups: Calculate, ADME/Tox, Search, Analyze, and Transform.
* Scaffold Tree improvements:
  * [#1730](https://github.com/datagrok-ai/public/issues/1730): Implemented Scaffold Tree integration to the **Filters Panel**.
  * [#1998](https://github.com/datagrok-ai/public/issues/1998): Now a scaffold tree shows a confirmation dialog before dropping all trees on the **Clear** icon click.
  * Added the **Allow Generate** property for the viewer in order to control autogeneration.

### Bug Fixes

* GROK-13105: Substructure search doesn't work after similarity search.
* GROK-13118: Activity cliffs selects non-numeric column as activity.
* GROK-13123: Substructure Search: error when the Filter Panel is opened with substructure filtering.
* [#1492](https://github.com/datagrok-ai/public/issues/1492): Elemental analysis: malformed data handling.
* GROK-11898: Orientation for smiles in Structural Alerts.
* GROK-12115: Hamburger menu closing while switching a sketcher.
* GROK-12905: Sketcher is not opening from the **Filter Panel**.
* GROK-12929: Scripts don't work if called from the package.
* GROK-12933: Drug likeness: set score precision.
* GROK-12961: Elemental Analysis: `Unsupported operation: NaN.round()` error on some data.
* GROK-12962: Similarity Search: doesn't work on malformed data.