# Scaffold Tree — what the features cannot say yet

Written during the round of 26 Oct 2026 that turned closed issues into regression scenarios
(`scaffold-tree-regressions.feature`), with the existing steps only: no library, binding or code change.
Every scaffold there is sketched by hand, so none of them runs the Python tree generation.

## Still open

1. **The sketcher dialog's invalid state** (public#2139: a scaffold made invalid in the sketcher and fixed
   again did not filter). Needs a step that claims the molecule input of the sketcher dialog holds an invalid
   structure (the error state the sketcher shows for a SMILES RDKit cannot read), so that the fix is claimed
   after a real invalid state, not after a typo the sketcher ignored.

   ```gherkin
   Scenario: A scaffold made invalid and fixed again still filters (public#2139)
     When user picks "Chem > Analyze > Scaffold Tree" from the top menu
     And user hovers over Scaffold Tree viewer
     And user clicks on "Sketch scaffolds manually" icon inside Scaffold Tree viewer
     And user types "c1ccncc" into molecule input of sketcher dialog
     And user presses Enter in molecule input of sketcher dialog
     # NEW STEP: Then the molecule input of sketcher dialog should show an invalid structure
     When user clears molecule input of sketcher dialog
     And user types "c1ccncc1" into molecule input of sketcher dialog
     And user presses Enter in molecule input of sketcher dialog
     And user clicks on "Add" button in sketcher dialog
     Then the "hits of node 1" reading of Scaffold Tree viewer should be 17
     When user clicks on the "checkbox of node 1" area of Scaffold Tree viewer
     Then 17 rows should pass the filter
     And no errors should have been logged
   ```

2. **The count a node shows** (public#2550: the node counters did not follow the table's filter). `hits of node
   <i>` is the node's bitset count over the whole table, not the label the node draws. Wanted in
   `ScaffoldTreeViewer.getWidgetStatus()` (Chem `src/widgets/scaffold-tree.ts`): a `shown of node <i>` reading
   of the counter label.

   ```gherkin
   Scenario: Node counts follow the table's filter (public#2550)
     # a pyridine node sketched by hand, as in scaffold-tree-regressions.feature
     Then the "hits of node 1" reading of Scaffold Tree viewer should be 17
     When user remembers the "shown of node 1" reading of Scaffold Tree viewer
     And user filters rows where "Stereo Category" is "S_UNKN"
     # NEW SIGNAL (Chem): the counter label the node draws
     Then the "shown of node 1" reading of Scaffold Tree viewer should be lower than remembered
   ```

   The exact count is to be measured once the reading exists.

3. **The filter panel's hover text** (public#2504: the panel's ? did not describe the Scaffold Tree filter as it
   does other filters). Needs a step that hovers the filter indicator of the filter panel and reads its tooltip
   text; the panel's `summary of <caption>` reading is what the counter's tooltip lists, but nothing yet claims the
   tooltip itself. The Scaffold Tree filter card's readings are also needed (see `filters/MISSING.md`, item 1).

   ```gherkin
   Scenario: The panel's hover describes the Scaffold Tree filter (public#2504)
     # a Scaffold Tree filter card with a checked pyridine node
     # NEW STEP: When user hovers over the filter indicator of filter panel
     Then tooltip should contain the text "Structure"
   ```

4. **The filter card variant of public#2135**: `filters/MISSING.md`, item 1.

## Seen while writing (not a closed issue, proposed for a check by hand)

- **A node checked right after it is added filters out every row, and stays so.** On spgi-100 with OpenChemLib:
  Chem | Analyze | Scaffold Tree, Sketch scaffolds manually, `c1ccncc1`, Enter, Add, and the node's checkbox
  clicked at once. The viewer then reads `hits of node 1` = 17 and `checked nodes` = 1, but 0 rows pass the
  filter for the whole 15 s poll (measured twice in a row on localhost, 26 Oct 2026, and in three scenarios of
  the first run). With the check made after `hits of node 1` has come to 17, 17 rows pass. Editing another
  node afterwards recomputes the filter (32 rows for a checked benzene). The regression scenarios wait for the
  hits before checking, as `scaffold-tree.feature` does. A person needs to click within the node's search time
  to meet it; the suspect is the filter applied from the node's bitset before its search has set it.
