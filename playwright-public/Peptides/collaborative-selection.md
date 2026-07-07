---
feature: peptides
sub_features_covered:
  - peptides.model.fire-bitset-changed
  - peptides.rendering.weblogo-header
  - peptides.util.modify-selection
  - peptides.util.get-selection-bitset
  - peptides.widgets.distribution
  - peptides.widgets.selection
target_layer: playwright
coverage_type: regression
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

## Setup

1. Open the linked Peptides demo dataset (`System:DemoFiles/bio/peptides.csv`) — a TableView with the peptides grid and a `Macromolecule`-semtype `AlignedSequence` column should be the active view.
2. Launch SAR analysis via top-menu `Bio | Analyze | SAR...`, accept the default config in the **Analyze Peptides** dialog (Activity column auto-detected, Scaling `-lg`), click **OK**. The default SAR layout attaches: per-position columns split out on the peptides grid with WebLogo column-headers rendered (`peptides.rendering.weblogo-header`), `Sequence Variability Map` (`MonomerPosition`) viewer, `Most Potent Residues` viewer, `Logo Summary Table` viewer.
3. Confirm the underlying DataFrame's selection BitSet is empty (no rows selected) and that the `Distribution` and `Selection` property-panel widgets reflect the empty-selection state before exercising the selection backbone.

## Scenarios

### Scenario 1 — WebLogo column-header click propagates selection through PeptidesModel.fireBitsetChanged to all configured viewers and property-panel widgets

Atlas anchor: `interactions[collaborative-selection-across-viewers]` (`coverage_type: regression`, `related_bugs: [GROK-14298]`) + `critical_paths[collaborative-selection-sync]` (p1, `derived_from: public/packages/Peptides/src/utils/cell-renderer.ts#L312`). Exercises the WebLogo-header click entry point of the collaborative-selection backbone: `setWebLogoRenderer` mouse handler → `modifySelection` → `getSelectionBitset` → `PeptidesModel.fireBitsetChanged(source)` → cross-viewer + cross-widget propagation against the `DataFrame.selection` BitSet.

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

### Scenario 2 — Shift / Ctrl click on a second WebLogo monomer applies `modifySelection` set-semantics; selection mutation re-broadcasts to all subscribers

Atlas anchor: same `interactions[collaborative-selection-across-viewers]` entry + `peptides.util.modify-selection` (atlas description: "applies Shift / Ctrl semantics to update the unified Selection map"). Exercises the modifier-aware path of `modifySelection` and asserts that subsequent re-broadcasts via `fireBitsetChanged` keep the cross-surface mirror consistent on a multi-pick selection — the regression risk for GROK-14298 (filter-broadcast crash + perf combo) is that a re-broadcast can crash one listener while leaving others stale.

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

- **Atlas provenance.** Scenario 1 derives from `interactions[collaborative-selection-across-viewers]` (atlas `peptides.yaml`, `coverage_type: regression`) and from `critical_paths[collaborative-selection-sync]` (p1). The critical_path entry's `derived_from:` cites `public/packages/Peptides/src/utils/cell-renderer.ts#L312` (the `setWebLogoRenderer` installation site that wires the click handler). Scenario 2 derives from the same `interactions[]` entry, extended to cover the Shift / Ctrl modifier semantics that `peptides.util.modify-selection` (atlas description: "applies Shift / Ctrl semantics to update the unified Selection map") implements at `public/packages/Peptides/src/utils/misc.ts#L332`.
- **target_layer rationale.** Multi-step UI flow — DOM-event-driven clicks on canvas-rendered WebLogo column-headers, modifier-key combinations on subsequent clicks, observation of selection state across the peptides grid + a custom Datagrok viewer + two property-panel widgets, console / error-balloon absence assertion for the GROK-14298 crash class. `apitest` is not viable because (a) the entry point is a canvas mouse handler bound to `setWebLogoRenderer`, not a JS-API method, and (b) the cross-surface broadcast assertion fans out to UI-rendered widgets (`Distribution`, `Selection`) and a UI-rendered viewer (`Sequence Variability Map`), all of which require DOM-level observation. Same layer choice as sibling `peptides.md` (which exercises the WebLogo-header click → selection path on a narrower surface).
- **coverage_type rationale.** Atlas `interactions[collaborative-selection-across-viewers].coverage_type: regression` is canonical per the rule "When picking an `interactions[]` entry as seed... `coverage_type:` = the entry's `coverage_type:` (canonical)." The companion `critical_paths[collaborative-selection-sync].priority: p1` (severity axis) maps to `regression` per the p1 → regression mapping in the layer-selection heuristic, consistent with the canonical declaration.
- **Coverage contribution.** Per the Gate F coverage audit at `scenario-chains/peptides.yaml :: gate_f_verdict.gaps[1].proposal` (cycle 03 revision 9), this scenario realizes proposal sub-clause **(c)** (`collaborative-selection-across-viewers` atlas interaction — "the highest atlas-traceability after retired (b)"). New sub_features added to the section's union: `peptides.model.fire-bitset-changed`, `peptides.util.modify-selection`, `peptides.widgets.selection` (all previously uncovered, on the `model.*`, `util.*`, and `widgets.*` uncovered families enumerated in `gaps[1].uncovered_sub_feature_families[]`). `peptides.util.get-selection-bitset` was already added to the section's union by `sar-similarity-threshold-matrix.md`; re-coverage here is intentional (it is the projection step that connects `modifySelection` to the broadcast path the scenario asserts). `peptides.rendering.weblogo-header` was already in the union via `peptides.md` and `sar-similarity-threshold-matrix.md`; re-coverage here is the entry-point click handler being exercised. `peptides.widgets.distribution` was already in the union via `peptides.md` and `sar.md`; re-coverage here is the cross-surface mirror assertion.
- **Related bugs.** `GROK-14298` ("Peptide SAR Analysis: Filtering causes panel to crash and causes lags") is the single curated bug whose `affects_sub_features[]` intersects this scenario's `sub_features_covered[]` (`peptides.model.fire-bitset-changed` + `peptides.widgets.distribution`) per `bug-library/peptides.yaml#GROK-14298` and atlas `known_issues[GROK-14298]`. The bug's class (null safety in listeners + lack of debounce/incrementality in stats recompute on broadcast) is asserted by Scenario 1 step 6 (cross-viewer broadcast must not throw a null-receiver error) and Scenario 2 step 5 (re-broadcast on selection mutation must not throw on the second / third broadcast). This scenario does NOT exercise the filter-application path — the chain analyzer's `bug_match_attempts_skipped[GROK-14298]` notes that the existing migrated scenarios in this section do not declare a filter step; the dedicated filter-broadcast cross-cutting repro remains a backlog item per the chain.
- **Setup composition.** Setup step 2 launches SAR via the top-menu entry path to keep the scenario's selection-backbone assertions cleanly separated from the alternative context-panel `Launch SAR` entry path (covered by sibling `sar.md`). The selection-backbone surface is identical between the two entry paths; choosing the top-menu path here mirrors the sister atlas-driven scenarios (`export-mutation-cliffs.md`, `export-invariant-map.md`, `sar-similarity-threshold-matrix.md`) and avoids duplicating `sar.md`'s entry-flow coverage.
- **Cross-scenario coverage shape.** Scenario 1 exercises the single-click → first-broadcast path (the basic backbone assertion). Scenario 2 exercises the modifier-key path → re-broadcast on selection mutation (the regression-class assertion for GROK-14298's listener-mutation crash mode). Together they cover the documented `peptides.util.modify-selection` semantics (Shift additive, Ctrl toggle-off) without requiring a separate dedicated scenario per modifier.
- **Deferrals.** None. Scenario steps map to existing Playwright-driveable surfaces (canvas-rendered WebLogo header clicks with modifier-key sequences, `DataFrame.selection.trueCount` readback via the JS API, viewer / widget DOM-element rendering assertions, console / error-balloon absence assertion).
- **See:** `public/help/datagrok/solutions/domains/bio/peptides-sar.md#Interactive Features` (the atlas help_docs entry for `peptides-sar.md` carries the `Interactive Features` heading with `covers_sub_features: []` — feature-level navigation pointer); `public/help/datagrok/solutions/domains/bio/peptides-sar.md#Sequence Variability Map` (maps to `peptides.viewers.monomer-position`, the viewer whose selection-mirror is asserted in step 5 of Scenario 1 and step 4 of Scenario 2).

## Original trailing metadata

```json
{
  "order": 8,
  "datasets": [
    "System:DemoFiles/bio/peptides.csv"
  ]
}
```
