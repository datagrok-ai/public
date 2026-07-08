---
feature: peptides
target_layer: playwright
coverage_type: regression
priority: p1
realizes: [collaborative-selection-sync]
produced_from: atlas-driven
realized_as:
  - collaborative-selection-spec.ts
related_bugs:
  - GROK-14298
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
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
    timestamp: 2026-05-30T21:44:00Z
    spec_runs:
      - spec: collaborative-selection-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 33
        failure_keys: []
---

# Peptides — Collaborative selection: WebLogo clicks stay in sync across viewers and panels

Clicking a monomer in a WebLogo column header selects the matching rows, and that selection must propagate everywhere: the peptides grid, the Sequence Variability Map viewer, and the Distribution and Selection side-panel widgets. This guards against GROK-14298, where a broadcast selection update could crash a listener or leave some panels stale. Scenario 1 checks a single click; Scenario 2 checks Shift/Ctrl-modified multi-selection and toggle-off.

## Setup

1. Open the linked Peptides demo dataset (`System:DemoFiles/bio/peptides.csv`) — a TableView with the peptides grid and a `Macromolecule`-semtype `AlignedSequence` column should be the active view.
2. Launch SAR analysis via top-menu `Bio | Analyze | SAR...`, accept the default config in the **Analyze Peptides** dialog (Activity column auto-detected, Scaling `-lg`), click **OK**. The default SAR layout attaches: per-position columns split out on the peptides grid with WebLogo column-headers rendered (`peptides.rendering.weblogo-header`), `Sequence Variability Map` (`MonomerPosition`) viewer, `Most Potent Residues` viewer, `Logo Summary Table` viewer.
3. Confirm the underlying DataFrame's selection BitSet is empty (no rows selected) and that the `Distribution` and `Selection` property-panel widgets reflect the empty-selection state before exercising the selection backbone.

## Scenarios

### Scenario 1 — Clicking a WebLogo monomer selects rows and updates every connected viewer

Clicking a monomer glyph in a WebLogo column header should select the matching rows and broadcast that selection to every viewer and widget wired to the table.

1. On any per-position column on the peptides grid, locate a populated monomer stack in the WebLogo column-header (the stack must visibly carry at least one monomer glyph — sparse but non-empty WebLogo).
2. Click a single monomer glyph in that WebLogo column-header. The handler invokes `modifySelection` with the (monomer, position) item; `getSelectionBitset` constructs the row mask; `PeptidesModel.fireBitsetChanged('WebLogo')` broadcasts to all subscribed viewers.
3. Confirm the underlying `DataFrame.selection` BitSet now reports a non-zero `trueCount` (at least one row is selected — exactly those rows carrying the picked monomer at the picked position).
4. Confirm the peptides grid reflects the selection: the selected rows are visibly highlighted with the standard Datagrok row-selection styling.
5. Confirm the `Sequence Variability Map` (`MonomerPosition`) viewer reflects the selection — the clicked (monomer, position) cell shows the selection highlight; the viewer is not silently inert.
6. Confirm the `Most Potent Residues` viewer remains consistent with the selection — no thrown error, no `Cannot read ... on null` console error from the cross-viewer broadcast path (the GROK-14298 crash class is null-safety in listeners; this is the principal assertion for that regression target).
7. Confirm the `Distribution` property-panel widget reflects the selection — the activity-distribution histogram updates to reflect the selected-row subset (the widget reads from the same `DataFrame.selection` BitSet that `fireBitsetChanged` updated).
8. Confirm the `Selection` property-panel widget reflects the selection — the selection summary grid shows the selected rows (mirroring the WebLogo header for the selected position, per `getSelectionWidget(table, options)`).

Expected (assertion summary):
- `DataFrame.selection.trueCount > 0` after the WebLogo header click.
- All four secondary surfaces (grid row highlight, `Sequence Variability Map`, `Distribution` widget, `Selection` widget) reflect the selection update without requiring a manual refresh.
- No null-receiver runtime errors in the console or error-balloon during the broadcast (GROK-14298 crash class).
- The selection is consistent across surfaces — surfaces do not display contradicting selected-row sets.

### Scenario 2 — Shift-click and Ctrl-click extend or toggle the selection, and every re-broadcast stays consistent

Shift-clicking a second WebLogo monomer should add to the selection; Ctrl-clicking should toggle a pick off. Each change re-broadcasts to all viewers and widgets — the regression risk (GROK-14298) is that a re-broadcast crashes one listener while leaving others showing stale data.

1. Starting from the single-monomer selection established by Scenario 1, identify a second populated monomer glyph in a different WebLogo column-header (a different position column).
2. Shift-click the second monomer glyph in that WebLogo column-header. `modifySelection` applies the Shift modifier semantics (additive selection — both the original and the new (monomer, position) picks are present in the unified Selection map).
3. Confirm the underlying `DataFrame.selection` BitSet's `trueCount` updates per the additive semantics — the resulting selection is the union (or per-position intersection-then-union, per `modifySelection`'s set rules) of the two picks, not a replacement.
4. Confirm the `Sequence Variability Map` viewer highlights both picked (monomer, position) cells.
5. Confirm the `Distribution` widget re-renders against the new selection — the histogram updates to reflect the new selected-row subset; no thrown error during the second broadcast.
6. Confirm the `Selection` widget re-renders against the new selection.
7. Now Ctrl-click the second monomer glyph (toggle off): `modifySelection` removes that (monomer, position) from the Selection map; `fireBitsetChanged` re-broadcasts the reduced selection.
8. Confirm the `DataFrame.selection` BitSet returns to the Scenario-1 single-monomer state (only the original pick remains selected); all four secondary surfaces reflect the reduced selection.

Expected:
- `modifySelection` Shift / Ctrl semantics produce additive then toggle-off transitions on the unified Selection map; `getSelectionBitset` projection to `DataFrame.selection` matches.
- All secondary surfaces (`Sequence Variability Map`, `Distribution` widget, `Selection` widget) stay consistent across both broadcasts (Shift-add and Ctrl-toggle-off).
- No null-receiver runtime errors during either re-broadcast (GROK-14298 crash class on listener mutation).

## Notes

- **Related bug.** GROK-14298 ("Peptide SAR Analysis: Filtering causes panel to crash and causes lags") is the null-safety-in-listeners bug this scenario guards against: Scenario 1 step 6 and Scenario 2 step 5 both assert that a selection broadcast never throws a null-receiver error. The filter-application path itself (as opposed to the WebLogo-click selection path) is not exercised here; a dedicated filter-broadcast repro remains a separate backlog item.
- **Setup composition.** SAR is launched via the top-menu entry path rather than the Peptides context-panel Launch SAR button, mirroring sibling scenarios (`export-mutation-cliffs.md`, `export-invariant-map.md`, `sar-similarity-threshold-matrix.md`) and keeping this scenario independent of `sar.md`, which owns the context-panel entry-path coverage.
- **Scenario split.** Scenario 1 covers the basic single-click selection broadcast. Scenario 2 covers the Shift (additive) and Ctrl (toggle-off) modifier semantics and the re-broadcast that follows a selection mutation.
- **See:** `public/help/datagrok/solutions/domains/bio/peptides-sar.md#Interactive Features`; `public/help/datagrok/solutions/domains/bio/peptides-sar.md#Sequence Variability Map`.

## Original trailing metadata

```json
{
  "order": 8,
  "datasets": [
    "System:DemoFiles/bio/peptides.csv"
  ]
}
```
