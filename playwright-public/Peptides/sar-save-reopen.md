---
feature: peptides
sub_features_covered:
  - peptides.workflow.start-analysis
  - peptides.model.init
  - peptides.model
  - peptides.workflow.sar-dialog
  - peptides.viewers.monomer-position
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
related_bugs:
  - GROK-14461
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - sar-save-reopen-spec.ts
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
    timestamp: 2026-05-30T21:52:52Z
    spec_runs:
      - spec: sar-save-reopen-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 210
        failure_keys: []
---

## Setup

1. Open the linked Peptides demo dataset (`System:DemoFiles/bio/peptides.csv`) — a TableView with the peptides grid and a `Macromolecule`-semtype `AlignedSequence` column should be the active view.
2. Launch SAR analysis via top-menu `Bio | Analyze | SAR...`, accept the default config in the **Analyze Peptides** dialog (Activity column auto-detected, Scaling `-lg`), click **OK**. The default SAR layout attaches: per-position columns split out on the peptides grid with WebLogo column-headers rendered, `Sequence Variability Map` (`MonomerPosition`) viewer, `Most Potent Residues` viewer, `Logo Summary Table` viewer. Confirm `PeptidesModel.init(settings)` has executed (the `Distribution` and `Selection` property-panel widgets render against an active model on the table).
3. Configure a non-default analysis state to round-trip across project save/reopen — click a single populated monomer glyph in any WebLogo column-header to establish a non-empty `DataFrame.selection` BitSet (the `peptides.model.fire-bitset-changed` broadcast wires the selection into all SAR viewers), and verify the `Sequence Variability Map` viewer reflects the selection (the picked (monomer, position) cell shows the standard Datagrok selection highlight). This selection state, the `-lg` activity scaling chosen at the dialog, and the default-set of viewers attached by `startAnalysis` together form the round-trip surface.

## Scenarios

### Scenario 1 — Save Peptides project with SAR layout + selection + activity scaling, then reopen and confirm layout is fully restored (or cleanly re-initialized) without the GROK-14461 missing-layout regression

Atlas anchor: `interactions[sar-project-save-and-reopen]` (`coverage_type: regression`, `related_bugs: [GROK-14461]`, atlas description: "Save a Peptides project with viewers + selection + activity scaling configured, close it, reopen — layout, model state, viewer attachments, and MonomerPositionStats cache must restore (or cleanly re-initialize). Cross-feature Peptides ↔ Projects ↔ datasync ↔ d4 layout. Sister pattern with Bio GROK-19928.") + `edge_cases[GROK-14461]` (`coverage_type: regression`, atlas description: "IL-4R Peptide (IDP 6557) project doesn't restore its layout on reopen; PeptidesModel state and viewer configuration are not reconstructed from project save data. Same architectural class as Bio GROK-19928, Chem GROK-17595, PowerPack GROK-17451 / GROK-17109 — analysis state persistence across project save/reopen with datasync."). Exercises the full save/reopen lifecycle of a SAR-configured project from the launch entry path (`peptides.workflow.sar-dialog` → `peptides.workflow.start-analysis` → `peptides.model.init`) through datasync persistence to layout-application on reopen.

1. With the SAR analysis active and the Setup-step-3 selection established, save the current state as a Peptides project — from the top menu invoke `File | Save Project...` (or use the toolbar `[name="button-Save"]` Save Project entry); in the Save Project dialog assign a recognizable project name (e.g. `test-peptides-save-reopen-<timestamp>`) and click **Save**.
2. Confirm the save completes without error (no error-balloon, no console null-receiver exceptions, project appears in `Projects` browser or referenced view).
3. Close the active TableView (and the underlying SAR layout) — either close the TableView tab directly, or navigate away to the `Home` view and let Datagrok release the active project. The `DataFrame.temp['peptidesModel']` PeptidesModel singleton goes out of scope with the view.
4. Reopen the saved project — navigate to the `Projects` browser, locate the project saved in step 1, and double-click it (or right-click → `Open Project`). The project-open lifecycle fires: layout JSON is deserialized, datasync rehydrates the DataFrame, and `peptides.lifecycle.init` runs (or has run) so `MonomerWorks` / `TreeHelper` / `PeptideUtils` singletons are available before viewers re-attach.
5. Confirm the TableView opens with the peptides grid populated (the `AlignedSequence` column is present with `Macromolecule` semtype) — datasync rehydration succeeded.
6. Confirm the SAR layout is applied (the principal GROK-14461 regression assertion). The `Sequence Variability Map` (`MonomerPosition`) viewer must be visible and attached to the reopened TableView; the `Most Potent Residues` viewer must be visible and attached; the `Logo Summary Table` viewer must be visible and attached. WebLogo column-headers must be rendered on the per-position columns of the peptides grid (the SAR-layout signature). NONE of these surfaces silently disappears (the GROK-14461 failure mode is "the project opens without the expected layout").
7. Confirm the `PeptidesModel` singleton is either restored on the reopened DataFrame OR cleanly re-initialized — the `Distribution` and `Selection` property-panel widgets must render against an active model (no `PeptidesModel is null` / `Cannot read ... on null` console error from a viewer trying to bind to a missing model). This is the atlas-stated "model state restored OR cleanly re-initialized" acceptance per the `interactions[sar-project-save-and-reopen]` description.
8. Confirm the activity scaling choice from the saved project (`-lg`) is reflected on reopen — open the SAR settings dialog (via the `peptides.widgets.settings-dialog` surface), confirm the `Scaling` field reads `-lg`, and close the dialog without changing the value. (If the model was cleanly re-initialized rather than fully restored, the scaling choice must still match the saved setting, since `PeptidesSettings` is part of the project's saved JSON tag `TAGS.SETTINGS`.)

Expected (assertion summary):
- Save completes without error; project lands in the `Projects` browser.
- Reopen succeeds; `AlignedSequence` Macromolecule column is present on the rehydrated DataFrame.
- All three default SAR viewers (`Sequence Variability Map`, `Most Potent Residues`, `Logo Summary Table`) are attached to the reopened TableView and WebLogo column-headers are rendered on the per-position columns — the GROK-14461 missing-layout regression does NOT occur.
- `PeptidesModel` is present on the reopened DataFrame (`DataFrame.temp['peptidesModel']` non-null) and both property-panel widgets that bind to the model render without a null-receiver console error.
- Saved-project `Scaling` setting (`-lg`) round-trips intact.

### Scenario 2 — Reopened project preserves selection state OR exposes the documented "cleanly re-initialized" empty-selection path; no null-receiver crash on the post-reopen broadcast path

Atlas anchor: same `interactions[sar-project-save-and-reopen]` entry, extended to assert the selection-state branch of "restored OR cleanly re-initialized" acceptance. The atlas description names selection alongside layout and model as round-trip surfaces; the bug class is "PeptidesModel state and viewer configuration are not reconstructed from project save data" — selection is part of model state. The regression-class assertion is that the post-reopen `fireBitsetChanged` path (the first user-driven re-broadcast after reopen) does not crash on a partially-rehydrated or cleanly-re-initialized PeptidesModel.

1. Starting from the reopened project state established by Scenario 1, inspect the current `DataFrame.selection` BitSet `trueCount` — record whether the saved selection was restored (`trueCount > 0`, matching the Setup-step-3 single-monomer pick) or whether the model was cleanly re-initialized to an empty selection (`trueCount == 0`).
2. Confirm the observed selection branch matches one of the two atlas-permitted outcomes (restored exactly, OR cleanly empty) — if the BitSet `trueCount` is non-zero but does NOT match the saved selection, the round-trip is partially-broken and the regression is asserted as failed.
3. Confirm the `Sequence Variability Map` viewer reflects the observed selection state (selected cell highlighted if restored; no selection highlight if cleanly empty) — the viewer's selection-mirror must be consistent with the underlying BitSet on the reopened state.
4. Exercise the post-reopen broadcast path — click any populated monomer glyph in a WebLogo column-header on the reopened TableView. This fires the same `setWebLogoRenderer` → `modifySelection` → `getSelectionBitset` → `PeptidesModel.fireBitsetChanged('WebLogo')` chain established by the Setup, but now against the reopened (or re-initialized) PeptidesModel.
5. Confirm the broadcast does not crash — no error-balloon, no `Cannot read ... on null` / `Cannot read 'fire' on undefined` console error. This is the regression-class assertion for the GROK-14461-sister architectural class (broken-model state after reopen surfacing as a null-receiver crash on the first user broadcast).
6. Confirm `DataFrame.selection.trueCount > 0` after the click (selection updated) AND the `Sequence Variability Map` viewer reflects the new selection (cross-surface mirror consistency holds on the post-reopen broadcast path).
7. Confirm the `Distribution` and `Selection` property-panel widgets re-render against the new selection without throwing — the activity-distribution histogram updates to reflect the new selected-row subset; the selection summary grid lists the selected rows.

Expected:
- Selection state on reopen matches one of the two atlas-permitted outcomes (restored exactly OR cleanly empty); a mismatched non-empty selection is a round-trip failure.
- The first post-reopen user broadcast (Step 4) does not throw a null-receiver error on `PeptidesModel.fireBitsetChanged` or its listeners.
- All secondary surfaces (`Sequence Variability Map`, `Distribution` widget, `Selection` widget) stay consistent with the post-reopen selection update.

## Notes

- **Atlas provenance.** Both scenarios derive from `interactions[sar-project-save-and-reopen]` (atlas `peptides.yaml`, `coverage_type: regression`, `related_bugs: [GROK-14461]`) and from `edge_cases[GROK-14461]` (atlas `peptides.yaml`, `coverage_type: regression`, `derived_from: .claude/skills/grok-orchestrate-test-cycle/references/bug-library/peptides.yaml#GROK-14461`). The interaction entry's `sub_features:` list is the canonical seed (`peptides.workflow.start-analysis`, `peptides.model.init`, `peptides.model`); this scenario extends that seed with `peptides.workflow.sar-dialog` (the launch entry exercised by Setup step 2) and `peptides.viewers.monomer-position` (the layout-surface viewer whose reopened-attachment is the principal GROK-14461 assertion in Scenario 1 step 6) for healthy density and a tighter regression-class assertion. Sister architectural class — Bio GROK-19928, Chem GROK-17595, PowerPack GROK-17451 / GROK-17109 — analysis state persistence across project save/reopen with datasync; the scenario's broadcast-path assertion in Scenario 2 step 5 also catches the listener-mutation class shared with GROK-14298 in the post-reopen broken-state mode.
- **target_layer rationale.** Multi-step UI flow with cross-dialog state, project save + project reopen via the Datagrok Projects browser, DOM-level observation of SAR layout reattachment (three viewers + WebLogo column-headers), property-panel widget re-rendering on the reopened state, console / error-balloon absence assertion on the post-reopen broadcast. `apitest` is not viable because (a) the principal GROK-14461 regression assertion is that the saved-layout JSON is applied during project-open lifecycle and the SAR viewer DOM elements are reattached, which has no JS-API surrogate, (b) the Projects-browser open-project step is a cross-feature UI surface (`Peptides ↔ Projects ↔ datasync ↔ d4 layout`, per the atlas description), and (c) the cleanly-re-initialized branch requires DOM-level viewer-presence assertion to distinguish "model null, viewers gone" (regression) from "model fresh, viewers attached on default layout" (acceptable). Same layer choice as sibling regression scenarios (`peptides.md`, `sar.md`, `peptide-space.md`, `collaborative-selection.md`).
- **coverage_type rationale.** Atlas `interactions[sar-project-save-and-reopen].coverage_type: regression` is canonical per the rule "When picking an `interactions[]` entry as seed... `coverage_type:` = the entry's `coverage_type:` (canonical)." The companion `edge_cases[GROK-14461].coverage_type: regression` aligns. This is the standard regression-class shape for cross-feature persistence regressions per the sister-bug pattern.
- **Coverage contribution.** Per the Gate F coverage audit at `scenario-chains/peptides.yaml :: gate_f_verdict.gaps[1].proposal` (cycle 03 revision 11/12), this scenario realizes proposal sub-clause **(d)** (`sar-project-save-and-reopen` atlas interaction — "covers model.init plus addresses the GROK-14461 persistence regression"; Test Designer priority order names this as the top-line `.md` anchor due to cross-feature persistence backbone — sister-bug pattern with Bio GROK-19928 / Chem GROK-17595 / PowerPack GROK-17451 / GROK-17109). New sub_features added to the section's union: `peptides.model.init` (previously uncovered, on the `model.*` uncovered family enumerated in `gaps[1].uncovered_sub_feature_families[]`). `peptides.workflow.start-analysis` was already in the union via `peptides.md`, `sar.md`, `peptide-space.md`, `export-mutation-cliffs.md`, `sar-similarity-threshold-matrix.md`, `export-invariant-map.md`; re-coverage here exercises the launch-leg of the save/reopen round trip. `peptides.model` was already in the union via `peptides.md`, `sar.md` (model bind surfaces); re-coverage here is the model-rebind assertion on the reopened state. `peptides.workflow.sar-dialog` was already in the union via `peptides.md`, `sar.md`; re-coverage exercises the launch entry. `peptides.viewers.monomer-position` was already in the union via `sar.md`, `peptide-space.md`, `export-mutation-cliffs.md`, `export-invariant-map.md`, `sar-similarity-threshold-matrix.md`; re-coverage here is the principal "viewer reattached on reopen" GROK-14461 assertion (the layout-surface anchor). Net addition to section's coverage union: one new sub_feature (`peptides.model.init`), bringing the union toward 27/72 ≈ 37.5% — incremental progress toward the 70% threshold; cycle-03-revision-11 gap[1].proposal sub-clauses (e)–(f) remain outstanding for subsequent Test Designer dispatches.
- **Related bugs.** `GROK-14461` ("PDS | Layout is not applied in IL-4R Peptide (IDP 6557) project") is the single curated bug anchored here per `bug-library/peptides.yaml#GROK-14461` (`affects: [peptides.model.init, peptides.model, peptides.workflow.start-analysis]`, `priority: p2`, `status: regression-risk`, `test_coverage: needed`) and atlas `known_issues[GROK-14461]`. The bug's class (saved layout not applied on reopen, PeptidesModel not reconstructed from project save data) is asserted by Scenario 1 step 6 (SAR layout surfaces are attached on reopen — the direct regression check) and Scenario 1 step 7 (PeptidesModel is present and property-panel widgets bind without null-receiver error — the model-rebind check). The sister-bug architectural pattern (Bio GROK-19928 / Chem GROK-17595 / PowerPack GROK-17451 / GROK-17109) means a fix discipline applied to one of these tickets typically retires the others; this scenario empirically anchors the Peptides side of the pattern. The atlas-permitted "cleanly re-initialized" branch is preserved (Scenario 1 step 7 and Scenario 2 step 2) because the atlas description names "restored OR cleanly re-initialized" as both acceptable outcomes — the assertion is not "selection state survives exactly" but "the post-reopen state is not broken".
- **Setup composition.** Setup step 2 launches SAR via the top-menu entry path to keep the scenario's save/reopen assertions cleanly separated from the alternative context-panel `Launch SAR` entry path (covered by sibling `sar.md`). The save/reopen surface is identical between the two entry paths; choosing the top-menu path here mirrors the sister atlas-driven scenarios (`export-mutation-cliffs.md`, `export-invariant-map.md`, `sar-similarity-threshold-matrix.md`, `collaborative-selection.md`) and avoids duplicating `sar.md`'s entry-flow coverage. Setup step 3 (selection state setup) is required by the atlas-stated round-trip surface ("Save a Peptides project with viewers + selection + activity scaling configured") — without it the save/reopen exercise would not cover the selection-state branch that Scenario 2 asserts.
- **Cross-scenario coverage shape.** Scenario 1 exercises the layout + model + scaling round-trip (the principal GROK-14461 missing-layout regression assertion). Scenario 2 exercises the selection-state round-trip and the post-reopen broadcast path (the regression-class assertion that the rehydrated or re-initialized model is functional for the first user-driven broadcast). Together they cover the documented atlas surfaces ("layout, model state, viewer attachments, and MonomerPositionStats cache must restore") plus the implicit "post-reopen broadcast must not crash" follow-on that the sister-bug pattern (GROK-14298 listener-mutation crash class) and the GROK-14461 broken-model failure mode jointly anchor.
- **Deferrals.** None. Scenario steps map to existing Playwright-driveable surfaces: top-menu Bio | Analyze | SAR launch, default-config dialog accept, WebLogo column-header monomer click, top-menu File | Save Project / toolbar Save dialog, project name input, project close via TableView close / navigate-home, Projects-browser open via right-click / double-click, viewer DOM-presence assertion via `findViewer`-class lookup, `DataFrame.selection.trueCount` readback via the JS API, WebLogo column-header re-rendering observation, console / error-balloon absence assertion. The `MonomerPositionStats cache must restore` clause of the atlas description is implicitly covered by the SAR-layout reattachment assertion (the WebLogo column-headers and `Sequence Variability Map` viewer rely on `MonomerPositionStats`; their successful re-rendering implies the cache is rehydrated or recomputed) — no separate dedicated step is required because the cache is an implementation detail of the viewers, not an end-user-observable surface.
- **See:** `public/help/datagrok/solutions/domains/bio/peptides-sar.md#Saving and reopening Peptide SAR projects` (the atlas `help_docs` entry for `peptides-sar.md` includes save/reopen coverage at the workflow level; feature-level navigation pointer per the atlas `help_docs[].sections_relevant[]` flat-string form — no rich-object `covers_sub_features[]` mapping is present so this citation is feature-level, not sub_feature-specific).

## Original trailing metadata

```json
{
  "order": 9,
  "datasets": [
    "System:DemoFiles/bio/peptides.csv"
  ]
}
```
