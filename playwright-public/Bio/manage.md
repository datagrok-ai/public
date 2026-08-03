---
feature: bio
sub_features_covered:
  - bio.manage.libraries-view
  - bio.detector
target_layer: playwright
coverage_type: smoke
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
scope_reductions: []
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

# Bio | Manage | Monomer Libraries — ui-smoke

Single-action ui-smoke for the `Bio | Manage | Monomer Libraries`
top-menu (atlas `bio.manage.libraries-view`,
`manageLibrariesView` — `package.ts#L1359`). Opens a HELM dataset so
the Macromolecule detector wires the sequence column, then dispatches
the top-menu to surface the `Manage Monomer Libraries` View; verifies
the per-library checkbox listing and that each library checkbox
toggles independently.

## Setup

Dataset: `System.AppData/Bio/tests/filter_HELM.csv` (canonical HELM
fixture per chain analyzer for this section). The dataset is opened so
the Macromolecule detector (atlas `bio.detector`) classifies the
sequence column synchronously — required to exercise the Bio context
on the active table view.

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
   opens (`grok.shell.v.type === 'view'`, root `.ui-panel.grok-view`;
   atlas `bio.manage.libraries-view`). Verify the View opens with the
   `Manage Monomer Libraries` section header and the `Manage Duplicate
   Monomer Symbols` section header present.

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

- **Source-text fixes silently applied** during migration:
  - The original step 2 said "A dialog opens" — corrected to "A
    **View** ... opens". Per recon (`bio.md:223`,
    MCP-validated against Bio 2.27.9), the live build opens a full
    View (`grok.shell.v.type === 'view'`, root
    `.ui-panel.grok-view`), not a `.d4-dialog`. The dialog form
    `manageMonomerLibraries` (atlas `bio.manage.libraries-dialog`,
    `package.ts#L1351`) is a separate, non-top-menu function and is
    NOT what `Bio | Manage | Monomer Libraries` dispatches.
  - Original step 3 contained the typo "it's checkboxes" — corrected
    to "its checkboxes" (possessive, not contraction). Step also
    expanded from "check dialog output and it's checkboxes
    functionality" into the structured Scenario 2 above so the
    assertion surface is decidable.
- **Sub-features covered:**
  - `bio.manage.libraries-view` (atlas L680) — primary surface
    (`Bio | Manage | Monomer Libraries` top-menu →
    `manageLibrariesView`).
  - `bio.detector` (atlas L122) — touched in setup (Macromolecule
    classification of the HELM sequence column on dataset open).
- **No `related_bugs`:** `bug-library/bio.yaml` curated bugs all
  affect analyze / search / rendering / transform surfaces; none
  intersect `bio.manage.*`. Confirmed by chain-analyzer skip
  decisions (`bio.yaml`, chain `bug_match_attempts_skipped` block).
- **Helper coverage:** `bio.flow.monomer-libraries` already
  registered in `helpers-registry.yaml:243` →
  `.claude/skills/grok-browser/references/bio.md:223`. No new
  candidate helpers proposed.
- **See atlas** `bio.cp.manage-monomer-libraries` (p1) and
  `bio.cp.monomer-library-crud` (p1) for the broader critical-path
  framing; this scenario is the entry-point ui-smoke. Full CRUD
  (upload library JSON → edit metadata → delete → AppData reflects)
  is out of scope for this 3-step ui-smoke.
- **Unresolved ambiguity** (carried in frontmatter
  `unresolved_ambiguities`): chain analyzer noted that this section
  ui-smoke selection is a poor fit for Rule 1's canonical trigger
  ("create + share + delete a single feature-entity"); manage.md is
  the shortest non-bug/non-matrix/non-integration scenario in the
  Bio section but its surface is dialog/view affordance rather than
  entity lifecycle. Operator may want to author a dedicated
  `bio-ui-smoke.md` seeded by the smoke-eligible p0 critical paths
  (`bio.cp.detect-open-fasta`, `bio.cp.detect-open-helm`,
  `bio.cp.bio-service-surface-init`, etc.). Flagged here so Critic F
  coverage-mode picks up the gap.

---
{
  "order": 2
}
