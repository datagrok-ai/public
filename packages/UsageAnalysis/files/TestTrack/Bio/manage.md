---
feature: bio
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: [bio.cp.manage-monomer-libraries]
realizes: [bio.manage.monomer-libraries]
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/bio/manage.md
migration_date: 2026-05-31
source_text_fixes:
  - dialog-opens-corrected-to-view-opens
  - its-typo-fixed
  - dataset-renamed-sample-helm-csv-to-filter-helm-csv-canonical-bio-fixture
  - step-3-expanded-from-check-dialog-output-into-structured-scenario-2
candidate_helpers: []
realized_as:
  - manage-spec.ts
unresolved_ambiguities:
  - ui-smoke-rule-1-trigger-fit-poor
scope_reductions:
  - id: SR-01
    check: E-SCENARIO-RUNTIME-ALIGNMENT
    rationale: |
      Full library CRUD (upload library JSON -> edit metadata -> delete ->
      AppData reflects) is out of scope for this scenario; it is delegated to
      the sibling scenario `bio-manage-libraries-crud.md`. This scenario
      covers the management/smoke surface only.
    verdict_status: SCOPE_REDUCTION
related_bugs: []
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T12:00:00Z
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T12:30:00Z
    review_round: 1
    failure_keys: []
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T23:45:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T22:25:00Z
    spec_runs:
      - spec: manage-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 75
        failure_keys: []
---

# Bio | Manage | Monomer Libraries — open view & toggle library checkboxes

A quick smoke check of the **Bio | Manage | Monomer Libraries**
top-menu. Opens a HELM dataset so the Macromolecule detector wires
the sequence column, then dispatches the top-menu to open the
**Manage Monomer Libraries** view; verifies the per-library checkbox
listing and that each library checkbox toggles independently.

## Setup

Dataset: `System.AppData/Bio/tests/filter_HELM.csv` (canonical HELM
fixture for this section). The dataset is opened so the Macromolecule
detector classifies the sequence column synchronously — required to
exercise the Bio context on the active table view.

The library set is **server-state dependent**: assert structure
(per-library checkbox host present, Search box present, checkbox is
toggleable), not a fixed library count. Recon (`bio.md:223`) observed
`HELMCoreLibrary.json`, `NH2.json`, `polytool-lib.json`,
`sample-lib-Aca-colored.json` on the recon server (4 checkboxes), but
this count is not load-bearing on the assertion.

## Scenarios

### Scenario 1 — Open dataset and dispatch top-menu

1. Open `System.AppData/Bio/tests/filter_HELM.csv`. The Macromolecule
   detector classifies the sequence column (atlas `bio.detector`);
   the table view opens.

2. On the menu ribbon, go to **Bio** > **Manage** > **Monomer
   Libraries**. A **View** named `Manage Monomer Libraries`
   opens (`grok.shell.v.type === 'view'`, root `.ui-panel.grok-view`).
   Verify the View opens with the `Manage Monomer Libraries` section
   header and the `Manage Duplicate Monomer Symbols` section header
   present.

### Scenario 2 — Verify per-library checkbox listing and toggle

3. Verify the View renders:
   - at least one per-library checkbox host
     (`[name="input-host-<libfile>.json"]`); selector pattern per
     helpers-registry `bio.flow.monomer-libraries` /
     `.claude/skills/grok-browser/references/bio.md#L223`,
   - the library Search box (`[name="input-host-Search"]`),
   - a checkbox element (`input[type="checkbox"].ui-input-editor`)
     inside each library host.

   Toggle one library checkbox: verify the checkbox transitions from
   checked → unchecked (or unchecked → checked) on click without
   error. Toggle a second library checkbox independently to verify
   per-row state isolation.

## Notes

- Full CRUD (upload library JSON → edit metadata → delete → AppData
  reflects) is out of scope here — covered by
  `bio-manage-libraries-crud.md`.

---
{
  "order": 2
}
