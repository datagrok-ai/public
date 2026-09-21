---
feature: chem
realizes_atlas: []   # untyped: no matching atlas type
realizes: [chem.cp.panels-synthon-search]
priority: p2
target_layer: manual-only
coverage_type: regression
manual_only_reason: |
  The Synthon Search context panels are fully present and correctly wired on dev
  (Databases > Synthon Search > Substructure Search / Similarity Search accordion
  panes; inputs Space / Max hits / Include synthons / Cutoff), but the underlying
  RDKit SynthonSpaceSearch server compute (Chem:synthonSearchFunc -> synthon_search.py)
  does not return within any Playwright-viable budget on dev: the panel search spun
  for >9 min without reaching a grid or "No matches", and direct synthonSearchFunc
  calls timed out at 90s and 175s (live MCP recon 2026-08-20, dev.datagrok.ai, no
  console errors). Every Expected-results bullet in all three scenarios — including
  Scenario 3's "No matches" empty branch — is gated behind that same round-trip, so
  no assertable end-state is reachable and an automated spec could only false-RED.
  Execute manually against a warmed environment where synthon compute completes.
  Recon note: the molecule column in System:DemoFiles/chem/smiles.csv is
  `canonical_smiles`, not `smiles` (correct the step prose when automating).
  Operator decision needed: is dev's synthon compute stall a bug (file a ticket) or
  a cold-cache artifact to warm before restoring automation (set target_layer back
  to playwright)?
related_bugs: []
---

# Chem — Synthon Search context panels (Substructure and Similarity)

## Setup

1. Open the table `System:DemoFiles/chem/smiles.csv` via the Files browser
   (Browse > Files > Demo Files > chem > smiles.csv). The table contains a
   `smiles` column detected as semType Molecule.
2. Click the first cell of the `smiles` column to select it and open the
   Context Panel on the right side of the screen.
3. In the Context Panel, navigate to the Databases tab.
4. Confirm that the Databases tab contains the Synthon Search sub-section with
   at least two entries: Substructure Search and Similarity Search. If the
   Synthon Search sub-section is absent, the Chem package has no synthon-data
   files installed and the scenario cannot continue — record this as a
   prerequisite failure.

## Scenarios

### Scenario 1: Substructure Search panel queries the synthon space and respects control inputs

Steps:
1. In the `smiles.csv` table view, click the first cell of the `smiles` column
   to select it. The Context Panel shows the Databases tab.
2. In the Context Panel > Databases > Synthon Search, locate the Substructure
   Search sub-section and expand it if collapsed.
3. Verify that the Space selector in the panel lists "Syntons_5567.csv" (the
   only available synthon space) and that the Max hits field shows the default
   value 100. The search runs automatically after the widget loads.
4. Wait for the search to complete (the loading indicator disappears). Observe
   the results grid that appears inside the panel: count the visible rows and
   confirm the count is greater than zero and less than 100 (if the query
   returns 100 rows it may have hit the cap without error, which is also
   acceptable — the critical assertion is that the count is strictly greater than
   zero, confirming real hits were returned rather than an empty or error state).
5. Confirm that an "Open compounds as table" icon (arrow pointing into a
   square) is visible above the results grid.
6. Change the Max hits input from 100 to 10 (click the input field, clear it,
   type 10, then press Tab or Enter to trigger re-search). Wait for the search to
   complete. The results grid row count must now be at most 10.
7. Enable the "Include synthons" toggle (click the checkbox or toggle to turn it
   on). Wait for re-search. Inspect the results grid column headers: at least one
   column matching the pattern "synthon_N" (e.g. "synthon_1") must be present
   among the added columns, each rendered as Molecule cells.

Expected:
- The Substructure Search panel returns a non-empty results grid (row count > 0)
  for a typical SMILES molecule against the Syntons_5567 synthon space.
- The "Open compounds as table" icon is visible when results are present.
- Reducing Max hits to 10 causes the results grid row count to drop to at most 10.
- Enabling "Include synthons" adds synthon_N columns (semType Molecule) to the
  results grid.

### Scenario 2: Similarity Search panel applies the cutoff slider and opens results as a table

Steps:
1. With the `smiles.csv` table still open, click the first cell of the `smiles`
   column again if the selection was lost.
2. In the Context Panel > Databases > Synthon Search, locate the Similarity
   Search sub-section and expand it.
3. Wait for the initial search to complete (default cutoff 0.50, Max hits 100,
   Space "Syntons_5567.csv"). Observe the results grid row count and confirm it
   is greater than zero.
4. Note the row count at cutoff 0.50 as the baseline.
5. Locate the Cutoff slider in the Similarity Search panel. Drag the slider
   toward the right end to set Cutoff to approximately 0.90 (the slider label
   or description text next to it should read 0.90). Wait for re-search.
   The results grid row count after re-search must be less than or equal to the
   baseline row count observed at cutoff 0.50.
6. Click the "Open compounds as table" icon above the results grid. A new
   table view named "Synthon Similarity Search Results" must open in the
   Datagrok workspace (visible as a new tab in the main view area). Confirm
   that the new table's row count matches the row count shown in the panel
   results grid before the icon was clicked.

Expected:
- The Similarity Search panel returns a non-empty results grid at the default
  cutoff 0.50 for a typical SMILES molecule.
- Raising the Cutoff to 0.90 does not increase the row count relative to the
  baseline at 0.50 (tighter cutoff yields fewer or equal hits).
- Clicking "Open compounds as table" opens a new table named "Synthon Similarity
  Search Results" whose row count equals the panel grid row count.

### Scenario 3: Empty-result and error states do not leave the panel in an indefinite spinner

Steps:
1. In the `smiles.csv` table view, select a cell containing a complex or
   unusual molecule that is unlikely to match any entry in the Syntons_5567
   synthon space. A long peptide-like SMILES or a large macrocycle works; if
   no such cell is present in smiles.csv, type
   `C1CC2CCCCC2C(=O)N1C3=CC=C(C=C3)Cl` directly into the Sketcher of the
   Substructure Search panel's Space query field, or open a new table with
   that structure as the only row and click its cell.
2. In the Context Panel > Databases > Synthon Search > Substructure Search
   panel, wait for the search to complete. The panel must display the text
   "No matches" rather than a results grid or an indefinite spinning indicator.
   The panel must not remain in the loading state after the server response.
3. Verify the browser console (DevTools > Console) contains no uncaught error
   originating from the synthon-search widget code during or after the
   empty-result render.

Expected:
- When the substructure query matches no structures in the synthon space, the
  panel displays the text "No matches" explicitly and exits the loading state.
- No uncaught error appears in the browser console during the empty-result
  render path.

## Notes

- target_layer rationale: all three scenarios require clicking molecule grid
  cells to open the Context Panel, interacting with widget controls (Space
  selector, Max hits input, Include synthons toggle, Cutoff slider, Open-as-
  table icon), and then asserting on DOM grid row counts, column headers, and
  table view creation. These are DOM-driven interactions and state assertions
  that playwright handles natively; apitest cannot drive context-panel widget
  controls.
- Synthon space available on dev: Syntons_5567.csv is the only file under
  `public/packages/Chem/files/synthon-data/` as of 2026-08-20. If additional
  spaces are deployed, the Space selector will list more entries — the scenarios
  remain valid with any space that returns hits for typical SMILES structures.
- The Cutoff slider interaction in Scenario 2 Step 5 fires on the `mouseup`
  event (`cutoffInput.root.querySelector('.ui-input-editor')?.addEventListener('mouseup', …)`).
  The Automator must trigger a `mouseup` after dragging the slider handle to
  the target position to ensure the re-search fires.
- The "Open compounds as table" icon uses `ui.iconFA('arrow-square-down', …)`.
  The Automator should locate this icon by its aria-label attribute
  ("Open compounds as table") or by its Font Awesome class `fa-arrow-square-down`
  within the panel container.

---
{
  "order": 6,
  "datasets": ["System:DemoFiles/chem/smiles.csv"]
}
