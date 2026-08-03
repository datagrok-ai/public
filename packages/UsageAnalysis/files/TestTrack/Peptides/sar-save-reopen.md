---
feature: peptides
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [sar-project-save-and-reopen]
realizes: [bio.analyze.sar]
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

# Peptides — Save a SAR project and reopen it: layout, model state, and selection round-trip

Saving a Datagrok project with a SAR analysis configured (viewers, a row selection, and a non-default activity scaling) and reopening it should restore the full layout — or at least cleanly re-initialize it — without silently dropping viewers or crashing. This guards against GROK-14461, where a saved Peptides project's layout wasn't applied on reopen.

## Setup

1. Open the linked Peptides demo dataset (`System:DemoFiles/bio/peptides.csv`) — a TableView with the peptides grid and a `Macromolecule`-semtype `AlignedSequence` column should be the active view.
2. Launch SAR analysis via top-menu `Bio | Analyze | SAR...`, accept the default config in the **Analyze Peptides** dialog (Activity column auto-detected, Scaling `-lg`), click **OK**. The default SAR layout attaches: per-position columns split out on the peptides grid with WebLogo column-headers rendered, `Sequence Variability Map` (`MonomerPosition`) viewer, `Most Potent Residues` viewer, `Logo Summary Table` viewer. Confirm `PeptidesModel.init(settings)` has executed (the `Distribution` and `Selection` property-panel widgets render against an active model on the table).
3. Configure a non-default analysis state to round-trip across project save/reopen — click a single populated monomer glyph in any WebLogo column-header to establish a non-empty `DataFrame.selection` BitSet (the `peptides.model.fire-bitset-changed` broadcast wires the selection into all SAR viewers), and verify the `Sequence Variability Map` viewer reflects the selection (the picked (monomer, position) cell shows the standard Datagrok selection highlight). This selection state, the `-lg` activity scaling chosen at the dialog, and the default-set of viewers attached by `startAnalysis` together form the round-trip surface.

## Scenarios

### Scenario 1 — Save a Peptides project with SAR layout + selection + activity scaling, then reopen and confirm the layout is restored

Saves a project with viewers, a row selection, and activity scaling configured, closes it, and reopens it — layout, model state, and viewer attachments must restore (or cleanly re-initialize). This is the same architectural class of persistence bug seen elsewhere on the platform (Bio GROK-19928, Chem GROK-17595, PowerPack GROK-17451/GROK-17109): analysis state not surviving a project save/reopen round trip.

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

### Scenario 2 — Selection state either restores or is cleanly emptied on reopen, and the first post-reopen selection broadcast doesn't crash

Extends Scenario 1 to check the selection branch specifically: the reopened project's selection must either match the saved state exactly, or be cleanly empty — never partially or inconsistently restored. It also checks that clicking a WebLogo monomer right after reopening (the first user-driven selection broadcast) works cleanly against the rehydrated or re-initialized model.

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

- **Related bug — GROK-14461.** ("Layout is not applied in IL-4R Peptide project on reopen"). Scenario 1 step 6 (SAR viewers attached on reopen) is the direct regression check; step 7 (model present, widgets bind without a null-receiver error) is the model-rebind check. This is the same architectural bug class seen on other platform features (Bio GROK-19928, Chem GROK-17595, PowerPack GROK-17451/GROK-17109) — a fix to one of those tends to fix the others too. The acceptable outcomes are "selection restored exactly" or "cleanly re-initialized to empty" — not "state survives exactly," just "the post-reopen state isn't broken."
- **Setup composition.** SAR is launched via the top-menu entry path, keeping this scenario independent of the context-panel `Launch SAR` entry path covered by `sar.md`. The save/reopen behavior itself doesn't depend on which entry path was used.
- **See:** `public/help/datagrok/solutions/domains/bio/peptides-sar.md#Saving and reopening Peptide SAR projects`.

## Original trailing metadata

```json
{
  "order": 9,
  "datasets": [
    "System:DemoFiles/bio/peptides.csv"
  ]
}
```
