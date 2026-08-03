---
feature: peptides
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [manual-alignment-edits-applied-to-sequence]
realizes: []
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

# Peptides — Manual Alignment panel: edit a sequence, Apply and Reset

Selecting a monomer cell on the SAR grid surfaces a **Manual Alignment** panel in the context panel, where the aligned sequence can be edited by hand. Clicking **Apply** should rewrite the sequence column and the per-position columns for that row and refresh the SAR viewers; clicking **Reset** should discard any unsaved edit without touching the underlying data. This scenario also checks that the collaborative-selection broadcast still works after an alignment edit.

## Setup

1. Open the linked Peptides demo dataset (`System:DemoFiles/bio/peptides.csv`) — a TableView with the peptides grid and a `Macromolecule`-semtype `AlignedSequence` column should be the active view.
2. Launch SAR analysis via top-menu `Bio | Analyze | SAR...`, accept the default config in the **Analyze Peptides** dialog (Activity column auto-detected, Scaling `-lg`), click **OK**. The default SAR layout attaches: per-position columns split out on the peptides grid with WebLogo column-headers rendered, `Sequence Variability Map` (`MonomerPosition`) viewer, `Most Potent Residues` viewer, `Logo Summary Table` viewer. Confirm the `PeptidesModel` singleton is present on the active DataFrame (`DataFrame.temp['peptidesModel']` non-null) — the `manualAlignment` panel depends on `PeptidesModel.getInstance(currentDf)` returning a non-null controller in its Apply handler.
3. Navigate the peptides grid to a deterministic row anchor for the alignment edit — set the current row to a populated peptide sequence (e.g. row 0 of the peptides grid via `DataFrame.currentRowIdx = 0`). Record the row's original `AlignedSequence` value and the per-position values from the split position columns (`col('0')`, `col('1')`, ... per the position-column split that `startAnalysis` performed). These values are the round-trip-baseline that the Reset path must restore and that the Apply path must overwrite.

## Scenarios

### Scenario 1 — Manual Alignment panel surfaces for a Monomer cell; Apply rewrites the sequence and per-position columns and refreshes downstream SAR viewers

Select a Monomer cell, open the Manual Alignment context panel, edit the aligned sequence in the text area, and click Apply — the underlying sequence column should reflect the change and downstream SAR viewers should refresh.

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

### Scenario 2 — Reset re-binds the text area to the current column value without rewriting it, and selection broadcast keeps working after an alignment edit

Extends Scenario 1 to check two things: that clicking **Reset** never writes back to any column (it only re-populates the text area from the live column value), and that the WebLogo click → selection broadcast chain still works cleanly after the row's sequence has been edited.

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

- **Prerequisite.** The Manual Alignment panel's Apply handler depends on an active SAR model (`PeptidesModel`) on the table — without one, the panel section either doesn't render or Apply crashes. Setup launches SAR first so that prerequisite is satisfied.
- **No related bugs.** No curated bug currently anchors this panel. The closest sibling bug class is GROK-14357 (a Logo Summary Table settings change silently not taking effect) — a similar "action accepted, no visible effect" failure mode, but on a different viewer. This scenario's Apply/Reset assertions structurally defend against that same class of silent-failure bug on the Manual Alignment surface.
- **See:** `public/help/datagrok/solutions/domains/bio/peptides-sar.md#Interactive Features`.

## Original trailing metadata

```json
{
  "order": 10,
  "datasets": [
    "System:DemoFiles/bio/peptides.csv"
  ]
}
```
