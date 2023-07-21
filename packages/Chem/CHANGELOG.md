# Cheminformatics changelog

## 1.6.20 (2023-07-21)

This release focuses on improving analysis stability and usability.

*Dependency: datagarok-api >= 1.14.0*

### Features

* Set default values in all dialogs where appropriate.
* Unified the layout of the search panel results.
* Work with DBs APIs: Added links to sources, similarity search viewer style, reproduce for ChemicalSpace, Enamine.
* Performed **Elemental Analysis** refinement. Now the resulting column name includes the name of the column for which Elemental Analysis is performed.
* Added the ability to choose whether to search for MCS exact atoms or MCS exact bonds when performing R-Groups Analysis.
* Added the ability to avoid recalculating the similarity search when the current row is changed. To achieve this, we introduced the **Follow Current Row** checkbox in the settings of the similarity search viewer. 
* For proper handling of properties and rendering, we now check for smarts and molecular fragments separately.
* Provided the ability to copy data from the **Descriptors** and **Properties** tabs on the **Context Pane**.
* Moved **Descriptors** and **Fingerprints** from the  **Context Pane**  to the **Top Menu** ( **Chem** > **Calculate**).
* Added **Substructure Search** to the **Top Menu** ( **Chem** > **Search**)
* Modified the tooltip for dialog and drag-n-drop to prevent it from overlapping with the data.
* Implemented RDKit rendering for Chembl, ChemblAPI, PubChem, and DrugBank databases if OCL is used currently.
* Performed UI polishing. We’ve adjusted input field names, repositioned elements in certain dialogs, added tooltips, and organized the top menu items into groups: Calculate, ADME/Tox, Search, Analyze, and Transform.
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