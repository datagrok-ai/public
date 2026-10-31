# Chem filters — what the features cannot say yet

Written during the round of 26 Oct 2026 that turned closed issues into regression scenarios
(folded into `card-interplay.feature`, `card-state.feature` and `column-popup.feature`), with the existing steps only: no library, binding or code change.

## Still open

1. **A file-backed table refreshed from its file** (public#2647: Tables | File | Refresh cleared the substructure
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
  intersect and switch off on their own"). The panel's master switch, the other half of public#3328, is claimed in
  `card-interplay.feature`.
