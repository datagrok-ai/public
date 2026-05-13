# Migration Report — r-group-analysis.md

## Step mapping

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| Block 1 Step 1 — "Open linked datasets" (smiles + sar-small) | Setup steps 1–2 (provisioning) + Block A step 1 / Block B step 1 (each block opens its own dataset explicitly) | preserved (split for clarity) |
| Block 1 Step 2 — "Go to Chem > Analyze > R-groups analysis. A dialog opens." | Block A step 2 | preserved |
| Block 1 Step 3 — "In the dialog, click the MCS button" | Block A step 3 | preserved |
| Block 1 Step 4 — "Select the Visual analysis checkbox" | Block A step 4 | preserved (rephrased — checkbox defaults to checked per `r-group-analysis.ts:97`, so step verifies default state) |
| Block 1 Step 5 — "Click the OK button." | Block A step 5 | preserved |
| Block 1 Expected — "balloon 'No R Groups were found' for 'smiles' dataset, trellis plot for 'sar-small' dataset" | Block A step 6 verification (smiles branch) + Block B step 4 verification (sar-small branch) | preserved as verification (split into two blocks — the original Block 1 conflated two datasets into a single dual-expectation; migrated body separates them per dataset since each is a distinct invariant) |
| Block 2 Step 1 — "Open linked datasets (use 'sar-small' dataset)" | Setup step 1 (already provisions both) + Block B step 1 (opens sar-small) | preserved |
| Block 2 Step 2 — "Go to Chem > Analyze > R-groups analysis. A dialog opens." | Block B step 2 | preserved |
| Block 2 Step 3 — "In the dialog, click the MCS button." | Block B step 3 | preserved |
| Block 2 Step 4 — "Click the OK button. Expected: a trellis plot with the results." | Block B step 4 (action + verification) | preserved |
| Block 2 Step 5 — "Run R-groups analysis once more." | Block B step 5 | preserved |
| Block 2 Step 6 — "Click MCS." | Block B step 6 | preserved |
| Block 2 Step 7 — "In the dialog uncheck the Replace latest checkbox. Expected: the second trellis plot should be displayed. In the grid, there should be two sets of columns resulting from the R-Group analysis." | Block B step 7 (action + verification) | preserved |
| Block 2 Step 8 — "Run R-groups analysis once more." | Block B step 8 | preserved |
| Block 2 Step 9 — "Click MCS." | Block B step 9 | preserved |
| Block 2 Step 10 — "In the dialog check the Replace latest checkbox. Expected: the latest results (columns and trellis plot) are replaced." | Block B step 10 (action + verification) | preserved |
| Block 2 Step 11 — "Run R-groups analysis without clicking the MCS. Expected: balloon 'No core was provided'" | Block B step 11 (action + verification) | preserved |

## Decisions

- **Why this `target_layer`:** chose `playwright` per chain YAML `output_plan`
  (`scenario-chains/chem.yaml` line ~880, `target_layer: playwright`, rationale
  "the 'No R Groups were found' balloon invariant requires DOM
  balloon-notification assertion that playwright handles natively"). The
  scenario drives a dialog widget (R-Groups Analysis), interacts with the
  sketcher canvas via the MCS button, toggles checkboxes, and asserts
  balloon-notification DOM — all classic playwright-mode UI driving per
  Migrator heuristic 1.
- **Why this `strategy`:** `simple` per chain `output_plan` (single migration
  step, no decomposition, no fixture sharing across files).
- **Why this `coverage_type`:** chose `edge` because this scenario IS the bug
  reproduction path for `GROK-16329` (chain rev 2 explicit:
  `pyramid_layer: bug-focused`; chain notes: "Block 1 IS the bug repro; the
  scenario was authored to lock in the no-rgroups invariant"). The MCS
  empty-result branch on `smiles.csv` is a textbook edge case — empty result
  is a legitimate outcome of the analysis function and the test asserts the
  graceful-message contract instead of the happy-path "decomposition
  succeeded" outcome. Atlas critical path `chem.cp.r-groups-analysis` (p1)
  explicitly names this as "the sharp edge". `coverage_type: edge` organically
  satisfies the section-level A-STRUCT-02 requirement (first edge|perf
  scenario in the Chem section); NO SR-01 carry-forward needed.
- **Sibling tests consulted:**
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/r-group-analysis-spec.ts`
    (existing playwright spec per `existing-test-index.yaml`,
    `test_name: "Chem: R-Groups Analysis"`, `features_covered:
    [chem.r-group-analysis]`) — Automator will regenerate / extend from
    this migrated scenario.
  - `public/packages/Chem/src/tests/projects-tests.ts` (Chem category,
    `test_name: r-group-analysis` + `r-group-analysis-sync`, package layer) —
    these are JS API-side R-Group tests via `uses-grok.functions` /
    `uses-grok.dapi`; consulted to confirm a Playwright UI surface is the
    distinct delivery target (not duplicate of package-layer JS API tests).
  - Sibling migrated scenarios in same Chem section:
    `info-panels.md` (coverage_type: regression, playwright);
    `Advanced/scaffold-tree-functions.md`; `Advanced/scaffold-tree.md`
    (coverage_type: regression, playwright) — borrowed frontmatter shape and
    house-style "Setup / Scenarios > ### Block / Notes" structure.
- **Helpers reused:** (none) — `helpers-registry.yaml` has no helpers covering
  R-Groups dialog driving. No candidate helper proposed: the R-Groups dialog
  is invoked at most once or twice per spec and inline driving is more
  readable than a single-use helper.
- **Bug library consulted:** `yes — bugs [GROK-16329]`. Chem bug library is
  curated; this scenario IS the reproduction path for `GROK-16329` per chain
  rev 2 (Block A smiles + MCS + empty-result expectation). Atlas
  `known_issues` (lines 1646–1657) names `GROK-16329` with sub_features
  `[chem.analyze.r-groups, .top-menu, .decomposition]` — full intersection
  with scenario coverage. NOT emitted under
  `bug_focused_candidates[]` in the chain (chain rev 2 notes: "GROK-16329
  affects [all of which this scenario covers] — so bug is covered by this
  scenario; not emitted under bug_focused_candidates[]").
- **Decision log queried:** `yes — empty for r-group-analysis`. No
  `failed_attempts` entries match `feature == chem AND failure_key ==
  r-group-analysis` (decision-log.yaml has chem-section entries for
  scaffold-tree wave and atlas-curation infrastructure but none for
  R-Groups). No prior approaches off the table for this scenario.

### Source-text fix per chain rev 1 → 2 directive (Migrator authoritative cross-reference)

The chain rev 1 flagged the "Visual analysis checkbox" label as unverified
against the atlas. Per Olena's "Migrator fixes by cross-referencing src"
directive, I read
`public/packages/Chem/src/analysis/r-group-analysis.ts` and confirmed:

- **"Visual analysis"** label is correct (verbatim).
  `r-group-analysis.ts:97`: `const visualAnalysisCheck = ui.input.bool('Visual
  analysis', {value: true});` — checkbox label is `"Visual analysis"`,
  defaults to checked. Tooltip: "Add trellis plot after analysis is
  completed". Silent fix applied to migrated body: clarified Block A step 4
  to "Confirm the Visual analysis checkbox is checked (it defaults to true)"
  rather than the original "Select the Visual analysis checkbox", since the
  source confirms the default is `true` — selecting it would be a no-op on a
  freshly-opened dialog.
- **"Replace latest"** label is correct (verbatim).
  `r-group-analysis.ts:99`: `const replaceLatest = ui.input.bool('Replace
  latest', {value: true});`. Conditional rendering at `r-group-analysis.ts:162`
  — only shown when prior analysis exists for the dataframe. Captured in
  Notes section.
- **"No core was provided"** balloon text is correct (verbatim).
  `r-group-analysis.ts:241`: `grok.shell.error('No core was provided');`.
- **"No R-Groups were found"** balloon text — original scenario said "No R
  Groups were found" (with a space, no hyphen, Groups capitalized).
  `r-group-analysis.ts:192`: `grok.shell.error('No R-Groups were found');` —
  the canonical product wording is **"No R-Groups were found"** (hyphenated,
  capitalized "G"). Silent fix applied — migrated body uses the canonical
  hyphenated form. Captured in Notes.

### Bug-focused / GROK-16329 reproduction mapping

Block A (smiles.csv → MCS → Visual analysis → OK → expect `No R-Groups were
found` balloon) IS the literal reproduction path for `GROK-16329`. The bug
manifested as "Cannot set properties of null" — downstream code attempted to
set properties on a null decomposition result. The fix (1.20.0) introduced
graceful empty-result handling via
`r-group-analysis.ts:191–192`:

```ts
if (!res.yAxisColName && !res.xAxisColName)
  grok.shell.error('No R-Groups were found');
```

The migrated Block A locks this exact branch: structurally diverse SMILES
input → MCS returns a substructure but `rGroupDecomp` returns
`xAxisColName` / `yAxisColName` empty → balloon shown. Regression of this
fix would either silently swallow the empty case (no balloon) or re-throw the
null reference. `related_bugs: [GROK-16329]` set in frontmatter.

### Negative-path invariant — "No core was provided" balloon

Block B step 11 (run R-Groups → don't click MCS → click OK → expect "No core
was provided" balloon) is a distinct negative-path invariant guarding
`rGroupDecomp` lines 240–243:

```ts
if (!core) {
  grok.shell.error('No core was provided');
  return;
}
```

This is NOT GROK-16329 (which is about empty decomposition result, not empty
core input) but a sibling defensive-coding invariant in the same function.
Preserved as Block B step 11 verification.

## Opt-outs (SCOPE_REDUCTION proposals)

(none)

## Deferred items (NOT opt-outs)

(none)

## Edge cases

- **MCS returns a substructure but decomposition yields no R-groups (smiles.csv
  on diverse SMILES):** preserved as Block A step 6 verification. This IS the
  GROK-16329 invariant — the central edge case the scenario locks. No atlas
  movement needed (atlas `known_issues` line 1646–1657 already captures the
  pattern via the source_bug entry).
- **MCS-without-input (sketcher empty, OK clicked):** preserved as Block B
  step 11 verification (negative-path invariant — guards
  `rGroupDecomp` empty-core check at `r-group-analysis.ts:240–243`).
- **Sequential R-Group runs (Replace Latest matrix):** preserved as Block B
  step 5–10 verifications. Three permutations exercised (default OFF after
  first run is N/A since checkbox doesn't appear yet; explicit OFF on second
  run; explicit ON on third run). This is multi-pass-stateful UI driving —
  appropriate as in-scenario steps rather than a separate scenario.
- **Visual analysis checkbox unchecked (skip trellis plot):** NOT covered by
  this scenario. The original scenario does not exercise the
  Visual-analysis=OFF branch (in `r-group-analysis.ts:190` the
  `visualAnalysisCheck.value!` gate controls trellis-plot creation when
  results exist). Flagged for atlas curator as a `manual_only` candidate or
  follow-up scenario — not currently in the chain rev 2 dependency_graph for
  this section.
- **Column prefix override:** NOT exercised. The dialog has a "Column prefix"
  input (`r-group-analysis.ts:95`, default `"R"`). Not in original scenario;
  not flagged for new coverage — orthogonal axis that the existing
  package-layer JS API tests already exercise.

## Unresolved ambiguities

- **Block B step 4 "trellis plot" location.** The original Block 2 Step 4
  expected "a trellis plot with the results" but does not specify whether the
  plot replaces the table view, opens as a separate viewer dock, or appears
  inline. From the source code (`r-group-analysis.ts:196–200`,
  `view.trellisPlot({xColumnNames: [res.xAxisColName], yColumnNames:
  [res.yAxisColName]})`) the trellis plot is added to the table view via the
  `view.trellisPlot` method — likely as a dock-pane viewer. Migrated body
  uses "a trellis plot with R-group results is added to the view"
  (non-specific). Flag for QA pair review or spec-time UI inspection — the
  spec selector and dock target may need confirmation against the running
  build.
- **Block B step 7 column count.** Original says "there should be two sets of
  columns resulting from the R-Group analysis". Migrated body preserves "two
  sets of R-Group columns". The exact suffix pattern (e.g. `R1`, `R1_1` vs
  `R1`, `R2_1`) is determined by `r-group-analysis.ts:214–232`
  (`getPrefixIdx` logic + Core/R column naming at lines 300–310). The spec
  may need to assert column-count delta (e.g. `+N` columns) rather than
  exact names — flag for Automator at spec time.
