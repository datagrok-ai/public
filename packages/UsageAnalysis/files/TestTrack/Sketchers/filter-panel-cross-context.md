---
feature: chem
sub_features_covered: [chem.sketcher, chem.sketcher.ocl]
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
pyramid_layer: bug-focused
migration_date: 2026-06-02
original_path: public/packages/UsageAnalysis/files/TestTrack/Sketchers/filter-panel-cross-context.md
authored_date: 2026-06-02
realized_as: [filter-panel-cross-context-spec.ts]
related_bugs: [GROK-12581, GROK-12903, GROK-12905, GROK-14028, GROK-12505]
gate_verdicts:
  b:
    verdict: PASS
    cycle_id: cycle-2026-06-02-chem-filter-panel-cross-context
    timestamp: 2026-06-02T00:00:00Z
    runs: 3
    runs_passed: 3
    note: "UI-driven via Filter Panel: Block A apply (50->47) + clear (input empties, 50/50 restored, GROK-14028); Block B backend-switch propagation in the filter sketcher dialog (GROK-12581/12903); Block C reopen after churn (GROK-12905). Green 3/3 (1.4m, 2.0m, 1.6m). Requires forcing Chem init (probe Sketcher) before getFiltersGroup so the molecule column gets the substructure filter. Block D (GROK-12505) manual; GROK-12803 excluded."
  e:
    verdict: PASS
    cycle_id: cycle-2026-06-02-chem-filter-panel-cross-context
    timestamp: 2026-06-02T00:00:00Z
    failure_keys: []
---

# Chem | Sketcher — Filter Panel cross-context

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
   - **Expected (GROK-14028):** the input line in Sektcher is cleared **and** all rows are restored, on the Filter Panel .

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

### Block D — `.structure-filter-type = Categorical` (GROK-12505) 


1. Right-click the molecule (Structure) column header → **Column Properties**.
2. Add the `.structure-filter-type` tag, set it to **Categorical**.
3. Open the Filter Panel, right-click the panel and select Remove all
4. Close the Filter and and open it again.
   - **Expected:** a **categorical** filter tab for Structure. 
