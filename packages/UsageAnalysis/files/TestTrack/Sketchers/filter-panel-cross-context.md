---
feature: chem
target_layer: playwright
coverage_type: regression
priority: p1
realizes: [chem.cp.substructure-filter, chem.x.sketcher-backend-switch-propagation, chem.x.filter-panel-sketcher-reopen]
produced_from: atlas-driven
pyramid_layer: bug-focused
migration_date: 2026-06-02
original_path: public/packages/UsageAnalysis/files/TestTrack/Sketchers/filter-panel-cross-context.md
authored_date: 2026-06-02
realized_as: [ filter-panel-cross-context-spec.ts ]
related_bugs: [ GROK-12581, GROK-12903, GROK-12905, GROK-14028 ]
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-03-sketchers-migrate-03
    timestamp: 2026-06-03T17:00:00Z
    review_round: 1
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-04-sketchers-filter-panel-revalidation
    timestamp: 2026-06-04T14:35:29Z
    spec_runs:
      - spec: filter-panel-cross-context-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 106
        failure_keys: []
  e:
    verdict: PASS
    cycle_id: cycle-2026-06-02-chem-filter-panel-cross-context
    timestamp: 2026-06-02T00:00:00Z
    failure_keys: []
  f:
    verdict: PASS
    cycle_id: 2026-06-03-sketchers-migrate-03
    timestamp: 2026-06-03T16:00:00Z
    failure_keys: []
---

# Chem | Sketcher — Filter Panel: substructure filter & backend-switch sync

Tests the substructure filter card on the Filter Panel: sketching a structure to filter the grid,
then using Reset to restore all rows and clear the sketch. Also checks that switching the sketcher
backend (OpenChemLib, Ketcher, ChemDraw) from the Filter Panel's dialog keeps the rest of the UI —
the column's hamburger menu and the Context Pane — in sync, and that the filter's sketcher dialog
still reopens correctly after switching backends back and forth.

## Setup

1. Dataset with a molecule column: `System:AppData/Chem/tests/smiles-50.csv`.
2. **Open the Filter Panel only after semantic-type detection** — otherwise the molecule column
   does not get the substructure filter (filters.md timing note). The spec forces
   `currentSketcherType = OpenChemLib` before opening filters and sets
   `showFiltersIconsConstantly = true`.

## Scenarios

### Block A — apply & clear the substructure filter (baseline + GROK-14028)

1. Open the Filter Panel; the molecule column shows a **substructure** filter card with a
   **Sketch** link.
2. Click **Sketch** → the sketcher dialog opens. Type `c1ccccc1` (benzene) into the SMILES
   input, press **Enter**, click **OK**.
   - **Expected:** the grid is filtered — visible row count drops below the total (benzene is a
     substructure of most, not all, rows).
3. On the Filter Panel click Reset.
   - **Expected (GROK-14028):** the input line in Sketcher is cleared **and** all rows are restored, on the Filter Panel .

### Block B — backend switch on the Filter Panel propagates (GROK-12581 / GROK-12903)

1. Open the filter sketcher dialog.
2. In its hamburger ≡ menu, switch the sketcher to **Ketcher**.
   - **Expected (GROK-12581):** the current sketcher changes everywhere — the shared
     `currentSketcherType` (read by the column hamburger menu and the Context Pane) becomes
     Ketcher. **Actual (original bug):** it stayed unchanged while the radio toggled.
3. Switch again to **ChemDraw** — the shared current sketcher follows again (GROK-12903: the
   selection made on the Filter Panel is the one that takes effect, not a stale previous one).

### Block C — reopen the filter sketcher after backend churn (GROK-12905)

1. From the dialog, switch back to **OpenChemLib**, then **close** the dialog.
2. Reopen the sketcher **from the Filter Panel**.
   - **Expected (GROK-12905):** the sketcher opens again. **Actual (original bug):** it didn't
     open after switching backends back and forth.

Each automated block also asserts zero sketcher `console.error`.
