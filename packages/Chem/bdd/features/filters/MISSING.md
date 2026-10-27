# Chem filters — what the features cannot say yet

Written during the round of 26 Oct 2026 that turned closed issues into regression scenarios
(`substructure-filter-regressions.feature`), with the existing steps only: no library, binding or code change.

## Still open

1. **The Scaffold Tree filter card's readings in Chem's suite** (public#2135: editing one scaffold reset the
   checkboxes of the others, reported on the filter card). The viewer variant is in
   `scaffold-tree/scaffold-tree-regressions.feature`; the card hosts a Scaffold Tree viewer that no view owns,
   so `the "…" reading of {widget}` cannot reach it. UsageAnalysis has the steps (`bindings/nx.ts`: `the {string}
   reading of the scaffold tree filter should be {int}`, its node areas), read through the panel's filter
   object. Wanted: those steps in the library's filter-panel tier (or Chem's bindings), so that:

   ```gherkin
   Scenario: Editing one scaffold of the filter card keeps another one checked (public#2135)
     When user clicks on filter icon in toolbar
     And user picks "Add Filter | Scaffold Tree Filter..." from the viewer menu of filter panel
     And user clicks on the "cell 1 of x" area of grid viewer in "Select columns..." dialog
     And user clicks on OK button in "Select columns..." dialog
     And user hovers over filter panel
     And user clicks on "Sketch scaffolds manually" icon inside filter panel
     And user types "c1ccncc1" into molecule input of sketcher dialog
     And user presses Enter in molecule input of sketcher dialog
     And user clicks on "Add" button in sketcher dialog
     And user hovers over filter panel
     And user clicks on "Sketch scaffolds manually" icon inside filter panel
     And user types "c1ccccc1" into molecule input of sketcher dialog
     And user presses Enter in molecule input of sketcher dialog
     And user clicks on "Add" button in sketcher dialog
     # NEW STEP (library): the scaffold tree filter's readings and areas
     Then the "hits of node 2" reading of the scaffold tree filter should be 32
     When user clicks on the "checkbox of node 2" area of the scaffold tree filter
     Then 32 rows should pass the filter
     When user hovers over the "node 1" area of the scaffold tree filter
     And user clicks on the "edit icon of node 1" area of the scaffold tree filter
     And user clears molecule input of sketcher dialog
     And user types "c1ccc2ncccc2c1" into molecule input of sketcher dialog
     And user presses Enter in molecule input of sketcher dialog
     And user clicks on "Save" button in sketcher dialog
     Then the "checked of node 2" reading of the scaffold tree filter should be "true"
     And 32 rows should pass the filter
     And no errors should have been logged
   ```

   The signal is there: `ScaffoldTreeViewer.getWidgetStatus()` (`checked of node <i>`, `checked nodes`, the
   node areas) on `filter.viewer` of the panel's Scaffold Tree filter.

2. **A file-backed table refreshed from its file** (public#2647: Tables | File | Refresh cleared the substructure
   card's criterion). Needs two steps: spgi-100 opened from a file in the user's files (the workspace's "a file
   put into the user's files" fixture, then opened from Browse, so the table keeps its file), and the Refresh of
   the table's File pane in the context panel, ended on the table's reload.

   ```gherkin
   Scenario: Refreshing a file-backed table keeps the substructure card's criterion (public#2647)
     # NEW STEP: Given spgi-100 is opened from a file in the user's files
     When user clicks on filter icon in toolbar
     And user clicks on "Sketch" text in "Structure" filter card
     And user types "c1ccncc1" into molecule input of sketcher dialog
     And user presses Enter in molecule input of sketcher dialog
     And user clicks on OK button in sketcher dialog
     Then 17 rows should pass the filter
     # NEW STEP: When user refreshes the table from its file (the File pane's Refresh), waiting for the reload
     Then the "structure of Structure" reading of filter panel should be "c1ccncc1"
     And 17 rows should pass the filter
     And no errors should have been logged
   ```

   The issue's second case (a cloned view whose Structure card was closed, then Refresh) needs the same two steps.

## Dropped

- **public#2753** (a structure drawn on a switched-off card did not filter when the card was switched on): the
  premise is gone. A switched-off card's body is `aria-disabled` and takes no clicks (the card's Sketch link is
  under the disabled body, measured 26 Oct 2026), so nothing can be drawn on it — by design (Olesia,
  27 Oct 2026). A card switched off with a
  structure and on again filters again: `card-interplay.feature` ("The substructure and category cards
  intersect and switch off on their own"). The panel's master switch, the other half of public#3328, is a
  scenario of `substructure-filter-regressions.feature`.
