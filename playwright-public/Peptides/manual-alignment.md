---
feature: peptides
sub_features_covered:
  - peptides.panels.manual-alignment
  - peptides.widgets.manual-alignment
  - peptides.model.fire-bitset-changed
  - peptides.workflow.start-analysis
  - peptides.workflow.sar-dialog
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - manual-alignment-spec.ts
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-05-29-peptides-migrate-02
    timestamp: 2026-05-29T00:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-05-30-peptides-automate-02
    timestamp: 2026-05-31T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-05-30-peptides-automate-02
    timestamp: 2026-05-30T22:01:01Z
    spec_runs:
      - spec: manual-alignment-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 131
        failure_keys: []
---

## Setup

1. Open the linked Peptides demo dataset (`System:DemoFiles/bio/peptides.csv`) — a TableView with the peptides grid and a `Macromolecule`-semtype `AlignedSequence` column should be the active view.
2. Launch SAR analysis via top-menu `Bio | Analyze | SAR...`, accept the default config in the **Analyze Peptides** dialog (Activity column auto-detected, Scaling `-lg`), click **OK**. The default SAR layout attaches: per-position columns split out on the peptides grid with WebLogo column-headers rendered, `Sequence Variability Map` (`MonomerPosition`) viewer, `Most Potent Residues` viewer, `Logo Summary Table` viewer. Confirm the `PeptidesModel` singleton is present on the active DataFrame (`DataFrame.temp['peptidesModel']` non-null) — the `manualAlignment` panel depends on `PeptidesModel.getInstance(currentDf)` returning a non-null controller in its Apply handler.
3. Navigate the peptides grid to a deterministic row anchor for the alignment edit — set the current row to a populated peptide sequence (e.g. row 0 of the peptides grid via `DataFrame.currentRowIdx = 0`). Record the row's original `AlignedSequence` value and the per-position values from the split position columns (`col('0')`, `col('1')`, ... per the position-column split that `startAnalysis` performed). These values are the round-trip-baseline that the Reset path must restore and that the Apply path must overwrite.

## Scenarios

### Scenario 1 — Manual Alignment panel surfaces for a Monomer cell, Apply rewrites the sequence column + per-position columns, and downstream SAR viewers refresh against the updated sequence

Atlas anchor: `interactions[manual-alignment-edits-applied-to-sequence]` (`coverage_type: regression`, `related_bugs: []`, atlas description: "Select a Monomer cell, open the Manual Alignment context panel, edit the aligned sequence in the text area, click Apply — the underlying sequence column reflects the change, downstream SAR viewers refresh."). Exercises the full `peptides.panels.manual-alignment` → `peptides.widgets.manual-alignment` → `peptideUtils.getSeqHelper().getSeqHandler().splitter()` → `AlignedSequence.set(idx, newSequence)` + per-position-column `.set(i.toString(), idx, part)` write-back → `PeptidesModel.getInstance(currentDf).updateGrid()` chain.

1. With the SAR analysis active and the current row anchored per Setup step 3, click directly on a Monomer-semtype cell on the peptides grid — pick a populated cell in one of the per-position columns (e.g. column `0` of the current row). The cell's semtype is `Monomer` (set by `setMonomerRenderer` on the per-position columns during `startAnalysis`); selecting it pins the cell as the Datagrok context-panel input (`grok.shell.o`).
2. Open the Datagrok context panel (right-side panel) and confirm the **Manual Alignment** panel section is present — this is the `manualAlignment` `@panel` registered for `Monomer`-semtype cells (`peptides.panels.manual-alignment`). The panel was conditionally rendered because `PeptidesModel.getInstance(currentDf)` returned a non-null model (Setup step 2 guarantees this); without an active SAR model the panel section is absent.
3. Expand the **Manual Alignment** panel — the `manualAlignmentWidget` mounts: a text area (with `pep-textinput` class) pre-populated with the current row's `AlignedSequence` value (per `alignedSequenceCol.get(currentDf.currentRowIdx)`), a `Reset` icon-button (`pep-snippet-editor-icon pep-reset-icon`), and an `Apply` button labeled `Apply changes`. Confirm the text-area value matches the Setup-step-3 recorded original sequence value.
4. Edit the aligned sequence in the text area — modify a single monomer at a known position (e.g. replace the monomer at split-position index 2 with a different alphabet-valid monomer, keeping the dash-separator format intact so the `SeqHandler.splitter()` still produces a same-length split). The text-area value reflects the edited string.
5. Click **Apply** — the click handler fires the documented widget logic: `PeptideUtils.getSeqHelper().getSeqHandler(alignedSequenceCol)` resolves the sequence handler, `sh.splitter(newSequence)` produces the per-position split, `alignedSequenceCol.set(affectedRowIndex, newSequence)` rewrites the row's sequence, the per-position columns receive `currentDf.set(i.toString(), affectedRowIndex, part)` for each split index whose column exists, and `PeptidesModel.getInstance(currentDf).updateGrid()` refreshes the SAR grid against the new state.
6. Confirm the underlying `AlignedSequence` column reflects the edit — read back `alignedSequenceCol.get(affectedRowIndex)` via the JS API; the value must equal the edited string from step 4 (not the original baseline from Setup step 3).
7. Confirm the per-position columns reflect the edit — read back `currentDf.col('2').get(affectedRowIndex)` (the position-column for the edited split index); the value must equal the new monomer at split-position 2 (not the original baseline). For all other split-position indices, the per-position-column values must remain unchanged (the edit was a single-monomer substitution at index 2).
8. Confirm the SAR layout refreshed against the new sequence — the peptides grid re-renders the edited row's per-position cells with the new monomer color (each `Monomer`-semtype cell is rendered by `setMonomerRenderer`'s per-cell renderer keyed on the cell value); the `Sequence Variability Map` viewer's underlying `MonomerPositionStats` cache is invalidated and recomputed against the updated per-position columns on the next render pass (the `updateGrid()` call schedules a redraw that picks up the new column values).

Expected (assertion summary):
- The **Manual Alignment** panel section is present in the context panel when a `Monomer`-semtype cell is selected and an active `PeptidesModel` exists; absent when no active SAR model is on the table.
- The text-area pre-populates with the current row's original `AlignedSequence` value.
- After Apply, `alignedSequenceCol.get(affectedRowIndex)` equals the edited string.
- After Apply, the per-position column for the edited split index equals the new monomer; all other per-position columns at the same row are unchanged.
- The SAR grid re-renders the edited row's per-position cells reflecting the new monomer colors and the WebLogo column-headers regenerate against the updated split-column values.
- No null-receiver console error from the Apply handler (the `PeptidesModel.getInstance(currentDf)` call inside the Apply path resolves to a non-null controller per Setup step 2).

### Scenario 2 — Reset button restores the text-area to the original sequence without rewriting the column; selection broadcast across viewers stays consistent across the alignment edit round-trip

Atlas anchor: same `interactions[manual-alignment-edits-applied-to-sequence]` entry, extended to assert (a) the `Reset` path (text-area baseline restoration without column write-back) and (b) that selection state broadcast via `peptides.model.fire-bitset-changed` remains consistent on the edited row. The atlas description names the widget's Apply + Reset surface; the regression-class assertion is that the post-edit broadcast path stays functional and that the Reset path is non-destructive (does NOT rewrite the column on click).

1. Starting from the post-Apply state established by Scenario 1 (the edited row now carries the edited sequence in the `AlignedSequence` column and the corresponding per-position-column update), re-open the **Manual Alignment** panel for the same Monomer-semtype cell (same row, same per-position column) — the text-area now pre-populates with the EDITED sequence (per `alignedSequenceCol.get(currentDf.currentRowIdx)`), confirming the panel reads the live column value on each mount rather than caching a stale value.
2. Click the **Reset** icon-button — per the documented widget logic, the click handler runs `sequenceInput.value = alignedSequenceCol.get(currentDf.currentRowIdx)!`. Since `alignedSequenceCol.get(currentDf.currentRowIdx)` is the EDITED sequence (Scenario 1 wrote it back), the Reset action re-binds the text-area to the same edited value — it does NOT restore the Setup-step-3 baseline (Reset's contract is "re-bind text-area to current column value", not "restore previous-edit baseline"). Confirm the text-area value matches the post-Apply edited string and the `AlignedSequence` column value is unchanged by the Reset click (Reset does NOT call `.set(...)` on either the sequence column or the per-position columns).
3. Edit the text-area again — type a fresh substitution (different from the Scenario 1 edit) at a different split-position index (e.g. index 4). Do NOT click Apply yet.
4. Click **Reset** while the text-area carries the unsaved fresh edit — confirm the text-area re-binds to the post-Scenario-1 column value (the unsaved fresh edit is discarded), and confirm both the `AlignedSequence` column and the per-position columns are unchanged from the post-Scenario-1 state (Reset on an unsaved edit is the documented baseline-restoration shape).
5. Exercise the selection-broadcast surface on the edited row — click a populated monomer glyph in a WebLogo column-header for the edited position; this fires the `setWebLogoRenderer` → `modifySelection` → `getSelectionBitset` → `PeptidesModel.fireBitsetChanged('WebLogo')` chain. Confirm the broadcast does not crash on the post-edit model state (no `Cannot read ... on null` / `Cannot read 'fire' on undefined` console error on the listener path) — the `peptides.model.fire-bitset-changed` listener fan-out must remain functional after the column write-back from Scenario 1's Apply.
6. Confirm `DataFrame.selection.trueCount > 0` after the WebLogo click (the broadcast wired the selection into the DataFrame's BitSet), and confirm the `Sequence Variability Map` viewer reflects the new selection (cross-surface mirror consistency holds across the alignment edit round-trip).

Expected:
- On panel re-open after Scenario 1's Apply, the text-area pre-populates with the EDITED sequence (panel reads live column value, no stale cache).
- Reset re-binds the text-area to the current column value without writing back to any column; both `AlignedSequence` and per-position columns retain their post-Scenario-1 state.
- Reset on an unsaved edit discards the in-progress edit without mutating any column.
- The WebLogo column-header click after the alignment edit fires `fireBitsetChanged` without throwing; the broadcast updates `DataFrame.selection` and the `Sequence Variability Map` viewer mirrors the selection.

## Notes

- **Atlas provenance.** Both scenarios derive from `interactions[manual-alignment-edits-applied-to-sequence]` (atlas `peptides.yaml`, `coverage_type: regression`, `related_bugs: []`, atlas description quoted in each Scenario heading). The interaction entry's `sub_features:` list is the canonical seed (`peptides.panels.manual-alignment`, `peptides.widgets.manual-alignment`, `peptides.model.fire-bitset-changed`); this scenario extends that seed with `peptides.workflow.start-analysis` and `peptides.workflow.sar-dialog` because the `manualAlignment` panel's Apply handler depends on `PeptidesModel.getInstance(currentDf)` returning a non-null controller — without an active SAR model, the panel section either does not render or the Apply path null-receiver-crashes (the prerequisite ordering is part of the regression-class assertion shape). Both extension sub_features are already in the section's coverage union via `peptides.md`, `sar.md`, `sar-save-reopen.md`, etc.; re-coverage here exercises them as the launch-leg prerequisite for the manual-alignment surface.
- **target_layer rationale.** Multi-step UI flow with context-panel section discovery, panel expansion + widget mount assertion, text-area interaction (typing, value readback), button-click chain with column write-back side-effects, post-edit broadcast-path assertion against a mutated model state, and cross-viewer DOM rendering observation. `apitest` is not viable because (a) the principal regression-class assertion is that the `manualAlignment` panel section surfaces in the right-side context panel when a `Monomer`-semtype cell is selected (a UI-surface assertion with no JS-API surrogate — the panel registration is `@panel` keyed on semtype, not directly callable as a function), (b) the Apply handler's side-effects include `PeptidesModel.updateGrid()` which schedules a grid re-render that must be observed at the DOM layer to confirm the new monomer colors render, and (c) the post-edit broadcast surface in Scenario 2 step 5 requires WebLogo column-header click — a canvas-rendered UI surface with no JS-API surrogate. Same layer choice as sibling regression scenarios (`peptides.md`, `sar.md`, `collaborative-selection.md`, `sar-save-reopen.md`).
- **coverage_type rationale.** Atlas `interactions[manual-alignment-edits-applied-to-sequence].coverage_type: regression` is canonical per the rule "When picking an `interactions[]` entry as seed... `coverage_type:` = the entry's `coverage_type:` (canonical)." This is the standard regression-class shape for cross-widget state-mutation round-trips on the panel surface (analogous to the settings-dialog → column-add round-trip class shared with GROK-14357 silent-propagation patterns, though no curated bug currently anchors this specific surface — `related_bugs: []` per atlas).
- **Coverage contribution.** Per the Gate F coverage audit at `scenario-chains/peptides.yaml :: gate_f_verdict.gaps[1].proposal` (cycle 03 revision 14), this scenario realizes proposal sub-clause **(e)** (`manual-alignment-edits` scenario — atlas interaction `manual-alignment-edits-applied-to-sequence` — "covers panels.manual-alignment + widgets.manual-alignment"). New sub_features added to the section's coverage union: `peptides.panels.manual-alignment` (previously uncovered, on the `panels.manual-alignment` entry enumerated in `gaps[1].uncovered_sub_feature_families[]`) and `peptides.widgets.manual-alignment` (previously uncovered, on the `widgets.{activity-distribution, manual-alignment}` family). `peptides.model.fire-bitset-changed` was already in the union via `collaborative-selection.md`; re-coverage here is the post-edit broadcast-path assertion that the listener fan-out remains functional after the column write-back. `peptides.workflow.start-analysis` and `peptides.workflow.sar-dialog` were already in the union via `peptides.md`, `sar.md`, `peptide-space.md`, `export-mutation-cliffs.md`, `sar-similarity-threshold-matrix.md`, `export-invariant-map.md`, `sar-save-reopen.md`; re-coverage here is the launch-leg prerequisite for the manual-alignment panel surface (without an active SAR model, the panel section either does not render or the Apply path null-receiver-crashes). Net addition to section's coverage union: two new sub_features (`peptides.panels.manual-alignment`, `peptides.widgets.manual-alignment`), bringing the union toward 30/72 ≈ 41.7% — incremental progress toward the 70% threshold; cycle-03-revision-14 gaps[1].proposal sub-clause (f) remains outstanding for subsequent Test Designer dispatches (`peptide-sar-demo-dashboard` scenario covering `app` + `demos.macromolecule-sar-fasta`).
- **Related bugs.** None (`related_bugs: []`). The atlas interaction entry carries `related_bugs: []` and no curated bug in `bug-library/peptides.yaml` anchors the `peptides.panels.manual-alignment` / `peptides.widgets.manual-alignment` surfaces specifically. The closest sibling bug class is GROK-14357 (silent UI-state propagation failure — Logo Summary Table settings dialog) which shares the "action accepted, no error, no effect" failure mode but anchors a different viewer surface. The scenario's regression-class assertion shape (column write-back observable; post-edit broadcast path stays functional; Reset is non-destructive) is the structural defense against the silent-propagation class for this widget surface, independent of a current curated bug.
- **Setup composition.** Setup step 2 launches SAR via the top-menu entry path (matching the sister atlas-driven scenarios `export-mutation-cliffs.md`, `export-invariant-map.md`, `sar-similarity-threshold-matrix.md`, `collaborative-selection.md`, `sar-save-reopen.md`); the `manualAlignment` panel's prerequisite is "active `PeptidesModel` on the table" — both entry paths (top-menu and context-panel `Launch SAR`) satisfy this prerequisite identically, so the choice does not affect the manual-alignment assertion surface. Setup step 3 anchors the current row to a deterministic baseline (row 0 with a populated peptide sequence) so the Apply / Reset write-back assertions can read deterministic values back via the JS API.
- **Cross-scenario coverage shape.** Scenario 1 exercises the panel-discovery + Apply + column-write-back path (the principal `peptides.widgets.manual-alignment` Apply assertion; the principal `peptides.panels.manual-alignment` panel-discovery assertion). Scenario 2 exercises the Reset path (non-destructive baseline restoration) and the post-edit broadcast surface (`peptides.model.fire-bitset-changed` stays functional after the column write-back). Together they cover the atlas-stated surfaces ("the underlying sequence column reflects the change, downstream SAR viewers refresh") plus the implicit "Reset must not rewrite the column" and "post-edit broadcast must not crash" follow-ons that the widget's two-button interface and the SAR model's listener fan-out architecture jointly anchor.
- **Deferrals.** None. Scenario steps map to existing Playwright-driveable surfaces: per-position-column cell click (semtype-keyed context-panel input pinning), right-side context panel section discovery (`Manual Alignment` panel name match), panel expansion (collapsible-panel toggle), text-area value readback (`pep-textinput` DOM element via `.value`), button click (`Apply changes` button label match; `pep-reset-icon` class match for Reset), JS-API readback of `alignedSequenceCol.get(...)` and `currentDf.col(...).get(...)` via the cell renderer's column accessors, WebLogo column-header monomer click (canvas-rendered with mouse event), `DataFrame.selection.trueCount` readback via the JS API, console / error-balloon absence assertion. The `MonomerPositionStats` cache invalidation in Scenario 1 step 8 is an implementation-detail follow-on (the visible end-user surface is the SAR grid re-render of the edited row's per-position cells; cache invalidation is implied by that re-render and does not require a separate dedicated step).
- **See:** `public/help/datagrok/solutions/domains/bio/peptides-sar.md#Interactive Features` (the atlas `help_docs` entry for `peptides-sar.md` includes the `Interactive Features` heading with `covers_sub_features: []` — feature-level navigation pointer; no rich-object `covers_sub_features[]` mapping exists for `peptides.panels.manual-alignment` or `peptides.widgets.manual-alignment` specifically, so this citation is feature-level rather than sub_feature-specific).

## Original trailing metadata

```json
{
  "order": 10,
  "datasets": [
    "System:DemoFiles/bio/peptides.csv"
  ]
}
```
