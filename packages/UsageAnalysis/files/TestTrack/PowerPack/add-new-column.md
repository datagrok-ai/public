---
feature: powerpack
sub_features_covered:
  - powerpack.dialogs.add-new-column-func
  - powerpack.dialogs.add-new-column
target_layer: playwright
coverage_type: smoke
pyramid_layer: ui-smoke
ui_coverage_responsibility:
  - add-new-column-dialog
  - add-new-column-toolbar-icon
  - add-new-column-recent-activities
  - add-new-column-autocomplete
  - add-new-column-drag-n-drop-columns
ui_coverage_delegated_to: null
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/PowerPack/add-new-column.md
migration_date: 2026-05-23
source_text_fixes:
  - testtrack-star-icon-replaced-with-explicit-demog-open
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
related_bugs:
  - GROK-17109
  - GROK-17004
realized_as:
  - add-new-column-spec.ts
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: 2026-05-23-powerpack-migrate-02
    timestamp: 2026-05-23T00:00:00Z
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-05-23-powerpack-migrate-02
    timestamp: 2026-05-23T00:00:00Z
    review_round: 1
    failure_keys: []
    claims:
      - check: A-STRUCT-MECH-01
        status: PASS
        evidence: |
          Frontmatter parses as YAML and carries all four required
          fields with valid values. feature is powerpack. The
          sub_features_covered list has two entries (the
          add-new-column-func and add-new-column ids). target_layer
          is playwright. coverage_type is smoke. Both sub-feature
          ids resolve to atlas entries at lines 615 and 622 of
          feature-atlas/powerpack.yaml.
      - check: A-STRUCT-MECH-02
        status: PASS
        evidence: |
          Body carries the H2 Scenarios container and one H3
          scenario heading (Add a column via the dialog and
          autofill from Recent Activities). Single-scenario file
          with no missing headings.
      - check: A-STRUCT-MECH-03
        status: PASS
        evidence: |
          Eight numbered steps (1. through 8.) appear under the
          scenario heading, covering dialog open, UI sanity, resize,
          name entry, autocomplete plus drag-n-drop formula
          composition, OK click, reopen, and Recent Activities
          autofill verification.
      - check: A-STRUCT-MECH-04
        status: PASS
        evidence: |
          The single scenario is non-empty; all 8 numbered steps
          carry concrete actions and verifications, not placeholder
          bullets or empty markers.
      - check: A-STRUCT-MECH-05
        status: PASS
        evidence: |
          target_layer value playwright is within the allowed enum
          set (playwright, apitest, manual-only) per verdict-enums
          derived_enums.target_layer.canonical_values.
      - check: A-STRUCT-MECH-06
        status: PASS
        evidence: |
          coverage_type value smoke is within the allowed enum set
          (smoke, regression, edge, perf) per verdict-enums
          derived_enums.coverage_type.canonical_values; not a
          severity-axis (p0..p3) value.
      - check: A-STRUCT-03
        status: PASS
        evidence: |
          coverage_type label set at file-level frontmatter
          (coverage_type smoke) applies to the file single scenario.
          No conflicting per-heading marker. Value comes from the
          unified test-kind enum, not the severity axis.
      - check: A-STRUCT-04
        status: PASS
        evidence: |
          Setup factored into a Setup section (opening
          System:DemoFiles/demog.csv) that the scenario implicitly
          consumes. No duplicated setup steps inside the numbered
          scenario step list. Single-scenario file so the
          across-scenarios duplication concern is also trivially
          satisfied.
      - check: A-LAYER-ALIGN-01
        status: PASS
        evidence: |
          pyramid_layer ui-smoke is paired with coverage_type smoke
          as required by the alignment rule (ui-smoke MUST be smoke).
          No mismatch.
      - check: A-CONT-01
        status: PASS
        evidence: |
          Body references concrete artifacts only. Demog dataset
          path System:DemoFiles/demog.csv; real columns HEIGHT and
          WEIGHT; real Round platform function; real Add new column
          toolbar icon; real Recent Activities dialog control. The
          created column name New is a concrete literal. No
          angle-bracket, square-bracket, or generic stand-in
          placeholders appear. The parametric-marker exception is
          inapplicable here.
      - check: A-BUG-01
        status: PASS
        evidence: |
          Atlas known_issues (feature-atlas/powerpack.yaml lines
          1393-1414) carry GROK-17109 and GROK-17004 with
          test_coverage.exists false (semantically equivalent to
          the legacy needed literal) and both affect the scenario
          sub-features add-new-column and add-new-column-func.
          Both bugs are addressed under clause (a) by appearing in
          this scenario related_bugs frontmatter list. Clause (b)
          is additionally satisfied since the Notes section names
          both bug ids and explains how the smoke covers the
          dialog-headline path while delegating the full
          save+datasync+reopen invariant (GROK-17109) and
          complex-paste invariant (GROK-17004) to sibling
          AddNewColumn scenarios and chain-level
          bug_focused_candidates (confirmed at
          scenario-chains/powerpack.yaml lines 308-320). Other
          atlas known_issues entries (GROK-17451, GROK-17269,
          GROK-18656, etc.) affect formula-lines, power-search and
          widgets surfaces, not this scenario sub-features, and
          are owned by sibling scenarios per the chain
          bug_focused_candidates block.
      - check: A-MERIT-01
        status: PASS
        evidence: |
          No step or scenario opted out for effort or complexity
          reasons. The scope-delegation discussion in Notes cites
          a real atlas-level dependency (sibling AddNewColumn
          scenarios and chain-level bug_focused_candidates that
          own the deferred bug-reproduction invariants), not
          effort.
      - check: A-MERIT-02
        status: PASS
        evidence: |
          No TODO add later or to be done in next phase deferral
          markers in body or frontmatter. The Notes section
          explicitly references the prerequisite sibling scenarios
          (AddNewColumn/add-new-column.md,
          AddNewColumn/formula-refreshing.md,
          AddNewColumn/highlight.md) and chain-level
          bug_focused_candidates as the coverage paths owning the
          deferred invariants.
  e:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T09:42:13Z
    spec_runs:
      - spec: add-new-column-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 31
        failure_keys: []
---

# Add New Column — Demog smoke (toolbar + UI sanity + Recent Activities autofill)

Single-source UI smoke for the Add New Column dialog on the Demog
dataset. Exercises the toolbar entry path (`addNewColumnDialog`
editor-for `AddNewColumn`), the dialog's built-in UI sanity (no
overlapping content, resizable, tooltips), the autocomplete hint
mechanic, the column drag-n-drop into the formula field, and the
Recent Activities autofill after the dialog is reopened.

Chain witness role: this scenario is the PowerPack chain's smoke
witness (`smoke_scenario: add-new-column.md`) for the dialog UI
surface; the sibling AddNewColumn/ scenarios delegate their basic
dialog-open-and-add flows here per chain
`ui_coverage_plan.delegated_scenarios`.

## Setup

1. Open `System:DemoFiles/demog.csv` so the Demog table view is the
   active view (replaces the TestTrack "star" / "Open test data"
   icon — same dataset, explicit source per the JSON metadata block
   on the original scenario).

## Scenarios

### Add a column via the dialog and autofill from Recent Activities

1. **Open the Add New Column dialog.** Click the "Add new column"
   toolbar icon on the Demog table view ribbon. A dialog opens.
2. **Verify dialog UI sanity.** Hover across the entire dialog and
   verify:
   - No overlapping text anywhere in the dialog.
   - No unnecessary scrollbars, icons, or list contents extending
     beyond the dialog boundaries.
   - Tooltips on dialog controls render with clear descriptions.
3. **Verify the dialog resizes.** Resize the dialog larger and then
   smaller; the layout adjusts appropriately to both directions
   without clipping or overlap.
4. **Add a new column named "New" with formula
   `Round(${HEIGHT} + ${WEIGHT})`.**
   - Type the column name `New` in the name input.
   - Compose the formula in the editor using BOTH mechanics:
     - **Autocomplete hint.** Start typing `Rou` in the formula
       editor and accept the `Round` function from the
       autocomplete tooltip suggestions.
     - **Column drag-n-drop.** Drag the `HEIGHT` column header
       from the columns panel into the formula editor at the
       argument position; do the same for `WEIGHT` so the
       formula becomes `Round(${HEIGHT} + ${WEIGHT})`.
5. **Click OK.** A new column named `New` is added to the Demog
   dataset.
6. **Reopen the Add New Column dialog.** Click the "Add new
   column" toolbar icon again. The dialog reopens.
7. **Open Recent Activities and select the most recent entry.**
   Locate the Recent Activities icon in the dialog and click it;
   select the most recent activity (the entry corresponding to the
   `New` column added at Step 5).
8. **Verify the form autofills.** The dialog's name input is
   prefilled with `New` and the formula editor is prefilled with
   the formula composed at Step 4 (`Round(${HEIGHT} + ${WEIGHT})`).

## Notes

- **Chain context.** This is the PowerPack chain's
  `smoke_scenario` per
  `scenario-chains/powerpack.yaml :: ui_coverage_plan`. Owns all
  five flows in `ui_coverage_responsibility:` directly; no
  delegation.
- **Related bugs.** `GROK-17109` (calculated columns persist
  across save+datasync+reopen) and `GROK-17004` (paste handler
  crash on complex formulas) both touch
  `powerpack.dialogs.add-new-column-func` +
  `powerpack.dialogs.add-new-column`. This smoke exercises the
  dialog's headline path but does NOT walk the full
  save+datasync+reopen invariant (covered by
  `AddNewColumn/add-new-column.md` + `AddNewColumn/formula-refreshing.md`)
  nor the complex-paste invariant (covered by
  `AddNewColumn/highlight.md`). Cross-cutting candidates are
  emitted at the chain level (`bug_focused_candidates`).
- **Sibling spec.** A Playwright spec already exists at
  `public/packages/PowerPack/src/tests/add-new-column.ts` (see
  `existing-test-index.yaml`); house-style anchor for Automator
  when authoring the migrated scenario's `-spec.ts`.
- **Source-text fix.** The original referenced opening the Demog
  dataset via the TestTrack "star" / "Open test data" icon
  (TestTrack-runner-specific UI). The migrated scenario references
  the dataset explicitly as `System:DemoFiles/demog.csv` per the
  original's trailing JSON metadata block. Same dataset, explicit
  source.
- **Original trailing JSON metadata.** The original scenario ended
  with `{"order": 2, "datasets": ["System:DemoFiles/demog.csv"]}`.
  The `order` field is captured in chain
  `order_from_files`; the `datasets` field is captured by the
  explicit `System:DemoFiles/demog.csv` reference in Setup.
