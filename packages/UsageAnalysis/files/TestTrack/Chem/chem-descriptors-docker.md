---
feature: chem
realizes_atlas:
  - chem.cp.calculate-descriptors-docker
realizes:
  - chem.calculate.descriptors
priority: p1
target_layer: playwright
pyramid_layer: ui-smoke
ui_companion:
  - chem-descriptors-docker-ui.md
coverage_type: smoke
related_bugs:
  - id: GROK-17621
    status: fixed
expected_results:
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      After clicking OK in the Descriptors dialog, the active table gains the
      two selected descriptor columns, named exactly `MolWt` and `MolLogP`,
      appended to the grid. The row count is unchanged from before the
      calculation. No console errors fire during dialog open, OK click, or
      column-append.
  - anchor: "Scenario 1 Step 6"
    expectation: >-
      Each appended descriptor column contains non-empty numeric values for
      valid SMILES rows and is not an all-null or all-zero column (a column of
      all zeros or nulls indicates the chem-chem container did not return
      results).
realized_as:
  - chem-descriptors-docker-spec.ts
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-08-19-chem-new-05
    timestamp: 2026-08-19T00:00:00Z
    failure_keys: []
    review_round: 1
    claims:
      - check_id: A-STRUCT-MECH-01
        status: PASS
      - check_id: A-STRUCT-MECH-02
        status: PASS
      - check_id: A-STRUCT-MECH-03
        status: PASS
      - check_id: A-STRUCT-MECH-04
        status: PASS
      - check_id: A-STRUCT-MECH-05
        status: PASS
      - check_id: A-STRUCT-MECH-06
        status: PASS
      - check_id: A-STRUCT-03
        status: PASS
      - check_id: A-STRUCT-04
        status: PASS
      - check_id: A-LAYER-ALIGN-01
        status: PASS
      - check_id: A-CONT-01
        status: PASS
      - check_id: A-BUG-01
        status: PASS
      - check_id: A-MERIT-01
        status: PASS
      - check_id: A-MERIT-02
        status: PASS
  e:
    verdict: PASS
    cycle_id: direct-gate-e-2026-08-22-chem-descriptors-docker-r2
    timestamp: 2026-08-22T15:41:07Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: direct-gate-b-2026-08-22-chem-descriptors-docker-r2
    timestamp: 2026-08-22T15:35:00Z
    spec_runs:
      - spec: chem-descriptors-docker-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 17
        failure_keys: []
        run_mode: headless-cold
---

# Chem — Calculate Descriptors via Docker

Covers the `chem.cp.calculate-descriptors-docker` critical path: opening the Descriptors
dialog (Chem > Calculate > Descriptors), selecting descriptors from the descriptor tree,
submitting to the chem-chem Docker container, and asserting that descriptor columns are
appended to the active table.

The bounded-error contract (`chem.int.descriptors-docker-bounded`) — when the container is
unavailable the UI must resolve to a clear error within bounded time, never an infinite
spinner (GROK-17621, fixed in 1.27.0) — is **not** automated here. It is Scenario 2 below,
and it is executed by hand from `chem-descriptors-docker-ui.md`, which realizes that atlas id
as a manual-only scenario. This file's automated coverage is therefore the happy path only,
which is why `coverage_type` is `smoke` rather than `regression`.

## Setup

1. Open `System:DemoFiles/chem/smiles.csv`. Confirm the table has 1000 rows and that the
   `canonical_smiles` column is detected as a Molecule column (a manual executor sees this
   as structure images in the cells rather than raw SMILES text; the paired spec reads the
   row count and the Molecule semType, and waits for the grid canvas to attach — it does not
   inspect the painted pixels).
   Actuation channel: a manual run opens the file through Browse; the paired spec opens it
   programmatically (`grok.dapi.files.readCsv`) because file-open actuation is not the
   subject of this scenario — it is covered by `chem-import-export-formats`. Everything
   from the Chem top menu onward is driven through the UI in both channels.
2. Confirm the Chem package is loaded and the **Chem > Calculate** top-menu tree is
   registered. For Scenario 1 (happy path), the `chem-chem` Docker container must be
   provisioned and reachable — verify that the container is running (e.g. by checking the
   Datagrok Docker panel or observing that a prior Calculate > Descriptors call returns
   successfully in a sanity run). Scenario 2 exercises the container-unavailable path and
   requires either a stopped container or a test environment where the container is absent.

## Scenarios

### Scenario 1: Descriptors happy path — columns appended

Steps:

1. With `smiles.csv` open, navigate to **Chem > Calculate > Descriptors** in the top menu.
   The Descriptors dialog opens, presenting a descriptor category tree.
2. In the Descriptors dialog:
   - Click the **None** link first. The descriptor selection is sticky — it carries the
     previous run's choices forward — so clearing it is what makes the appended-column set
     predictable; without this, extra descriptor columns land in the grid. Confirm the
     dialog's selection summary reads **0 checked**.
   - Expand the **Descriptors** group and check **MolWt**, then expand the **Crippen** group
     and check **MolLogP**. These are the real tree nodes: there is no "Lipinski /
     Physicochemical" category and no descriptor named plain "LogP".
   - Confirm the selection summary now reads **2 checked**. Read the selection back from the
     summary, not from the checkbox elements — a tree-node click can leave the checkbox
     looking checked while the dialog model has not registered it.
3. Click **OK** to submit the descriptor calculation request.
4. Wait for the calculation to complete (a progress indicator may appear while the chem-chem
   container processes the batch). The dialog closes.
5. Verify that exactly the two descriptor columns selected in Step 2 — named `MolWt` and
   `MolLogP` — are appended to the active table view's grid. The table's row count must
   equal the original 1000 rows: descriptor calculation must not drop or duplicate rows.
6. Verify that each appended descriptor column contains non-empty numeric values across
   sampled rows, including the first row and the last row of the grid. An all-null or
   all-zero descriptor column indicates the container did not return results and is a
   failure, not a pass.

Expected:

- The Descriptors dialog opens cleanly with no console errors.
- After OK, the `MolWt` and `MolLogP` columns are appended to the grid; the grid row count
  is unchanged from the pre-calculation count.
- Each descriptor column carries numeric values, first and last row included; no console
  errors fire during dialog open, dialog submission, container call, or column-append.

### Scenario 2: Container-unavailable bounded error (GROK-17621 regression guard)

**Delegated to `chem-descriptors-docker-ui.md` — manual-only, not automated here.** The
chem-chem container is provisioned and reachable on the dev validation stand, and the
Playwright harness cannot stop a platform-managed container to establish the "unavailable"
precondition; recon on 2026-08-19 saw the call succeed in ~4.5 s with no error balloon.

Run this scenario from the companion, which is authoritative for its steps and expected
results; they are deliberately not restated here. In outline: with the container down the
failure surfaces when the Descriptors dialog **opens**, not when OK is pressed. The
descriptor tree is fetched from the container while the dialog builds, so the dialog comes
up with "Could not load descriptors. The Chem service may be unavailable." where the tree
should be, and OK never enables — with no tree nothing is selected, so the dialog never
becomes valid. That bounded resolution, together with an unchanged table column count, is
what the GROK-17621 guard asserts.
