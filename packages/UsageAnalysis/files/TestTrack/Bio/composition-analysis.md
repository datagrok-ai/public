---
feature: bio
sub_features_covered:
  - bio.analyze.composition
  - bio.viewers.web-logo
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/bio/composition-analysis.md
migration_date: 2026-05-31
realized_as:
  - composition-analysis-spec.ts
source_text_fixes:
  - heading-normalized-to-bio-analyze-composition-integration
  - dataset-list-formatting-trailing-commas-cleaned
  - menu-path-suffix-analysis-dropped-per-atlas-canonical-form
  - implicit-three-dataset-matrix-made-explicit
  - step-2-viewer-presence-promoted-to-explicit-assertion
  - step-3-click-letter-rephrased-as-viewer-to-grid-bridge-assertion
  - dataset-names-substituted-sample-to-filter-per-chain-canonical-bio-fixtures
  - dataset-folder-substituted-samples-to-tests-per-chain-canonical-bio-fixtures
candidate_helpers:
  - bio.flow.composition
unresolved_ambiguities:
  - context-pane-property-checklist-not-defined-in-atlas
  - weblogo-canvas-interactivity-settle-window-not-atlas-declared
  - click-letter-canvas-hit-test-coordinates-not-atlas-declared
scope_reductions:
  - id: SR-01
    check: A-CONT-01
    rationale: |
      Source Step 5 ("change arbitrary properties") is non-deterministic —
      atlas does not enumerate a canonical Context-Pane property checklist
      for the Composition WebLogo viewer (parallels `analyze.md`
      SR-02). The Gear → Context-Pane wiring is preserved as a step
      (the Gear icon opens the Composition viewer's property surface and
      at least one property accepts a user-driven edit) but the
      correctness assertion on the edited-property state is deferred
      until atlas or operator supplies a concrete checklist of
      Context-Pane properties to verify; otherwise the verification
      would silently accept no-op edits.
    verdict_status: SCOPE_REDUCTION
related_bugs: []
gate_verdicts:
  a:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T16:30:00Z
    review_round: 1
    failure_keys: []
    scope_reduction_proposal: |
      Independent fresh-review verdict (cycle 2026-06-01-bio-migrate-02).

      All mechanical checks PASS. Frontmatter parses cleanly as YAML and
      carries the four required fields (feature=bio,
      sub_features_covered=[bio.analyze.composition, bio.viewers.web-logo],
      target_layer=playwright, coverage_type=regression). Three scenario
      headings present (Scenario 1, Scenario 2, Scenario 3), each with at
      least one numbered step (Scenario 1: steps 1-2; Scenario 2: step 3;
      Scenario 3: steps 4-5). No empty scenarios. target_layer=playwright
      and coverage_type=regression both belong to the unified enums.

      Structural-semantic checks PASS. coverage_type=regression at
      frontmatter level satisfies A-STRUCT-03 across all three scenarios.
      The 3-dataset matrix and the multi-column-dialog-absent precondition
      are factored into the dedicated `## Setup` block rather than
      duplicated under each scenario, satisfying A-STRUCT-04.

      A-LAYER-ALIGN-01 PASS-by-vacuity. The frontmatter has no
      pyramid_layer key (the chain YAML composition-analysis.md block
      assigns pyramid_layer=integration, but the hard alignment rule
      under A-LAYER-ALIGN-01 only triggers when the SCENARIO frontmatter
      itself declares pyramid_layer=ui-smoke; absent key returns PASS).

      A-CONT-01 trigger still present in body; carry-forward returns
      SCOPE_REDUCTION (not PASS). Scenario 3 Step 5 verifies "edit at
      least one editable property of the Composition WebLogo viewer" /
      "at least one property accepts a user-driven edit" — a
      generic-stand-in pattern that does not name a concrete WebLogo
      property. frontmatter scope_reductions[] carries SR-01 with
      check=A-CONT-01 and verdict_status=SCOPE_REDUCTION (bound by a
      prior cycle). Per the carry-forward rule, PASS requires both
      verdict_status=SCOPE_REDUCTION AND that the body no longer carry
      the trigger. The body intentionally preserves Steps 4-5 as the
      assertable-surface portion of the SR (the Gear-icon opens the
      Composition viewer's property surface and at least one property
      accepts a user-driven edit) while deferring the property-name
      correctness assertion. Trigger present → carry forward returns
      SCOPE_REDUCTION this cycle; verdict_status field on the SR is
      not mutated (forensic record).

      The cited dependency is real and atlas-level. atlas
      bio.viewers.web-logo (package.ts#L449) describes the WebLogo
      property surface as sequenceColumnName, positionWidth,
      startPosition, endPosition, mode, etc. but pins no canonical
      Context-Pane checklist of which subset MUST be present after
      Gear-click; atlas bio.analyze.composition (package.ts#L1041) does
      not enumerate one either. The dependency is also surfaced in
      unresolved_ambiguities[] as
      context-pane-property-checklist-not-defined-in-atlas, and is
      acknowledged in the body Notes block as a parallel to the
      analyze.md deferral on the same atlas surface.

      A-BUG-01 PASS-by-vacuity. atlas known_issues holds 8 bug entries
      (GROK-16111, GROK-12164, GROK-18474, GROK-15176, GROK-18616,
      GROK-16596, GROK-19928, GROK-19150) and all 8 use the schema form
      `test_coverage: { exists: false, paths: [] }`. Zero entries use
      the literal `test_coverage: needed` that A-BUG-01 keys on. The
      one WebLogo-affecting bug (GROK-18474, MSA column-header click on
      FASTA) lists bio.viewers.web-logo in affects_sub_features but
      does NOT use the needed literal; it is also excluded by chain
      bug_match_attempts_skipped[GROK-18474] (this scenario does not
      click MSA column headers) and the body Notes block explicitly
      acknowledges the exclusion.

      A-MERIT-01 PASS. SR-01 cites a real atlas-level dependency
      (Context-Pane property checklist not enumerated in
      bio.viewers.web-logo or bio.analyze.composition atlas entries),
      not effort or complexity. A-MERIT-02 PASS. SR-01 is not an
      unprompted deferral; it names the blocking dependency and the
      bind condition (atlas or operator supplies a concrete
      Context-Pane property checklist).

      Cross-scenario predicates (chain-wide 3+ sub-feature interaction,
      edge/perf presence, >=70% atlas coverage, density across chain)
      are Gate F's responsibility and not evaluated here.
    claims:
      - check: A-STRUCT-MECH-01
        status: PASS
      - check: A-STRUCT-MECH-02
        status: PASS
      - check: A-STRUCT-MECH-03
        status: PASS
      - check: A-STRUCT-MECH-04
        status: PASS
      - check: A-STRUCT-MECH-05
        status: PASS
      - check: A-STRUCT-MECH-06
        status: PASS
      - check: A-STRUCT-03
        status: PASS
      - check: A-STRUCT-04
        status: PASS
      - check: A-LAYER-ALIGN-01
        status: PASS
      - check: A-CONT-01
        status: FAIL
      - check: A-BUG-01
        status: PASS
      - check: A-MERIT-01
        status: PASS
      - check: A-MERIT-02
        status: PASS
  d:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-01T16:10:00Z
    failure_keys: []
    scope_reduction_proposal: >
      Independent fresh-review verdict (cycle 2026-06-01-bio-migrate-02).


      Mechanical checks all PASS:
        - D-STRUCT-MECH-03 — all 8 required migration-output
          frontmatter fields present (feature, sub_features_covered,
          target_layer=playwright, coverage_type=regression,
          produced_from=migrated, original_path, migration_date,
          related_bugs); deprecated `migrated_from` is absent.
        - D-STRUCT-MECH-05 — `original_path` resolves on disk
          (case-insensitive on Windows; on-disk casing is capital-B
          `Bio/`). Note: this scenario is in the Migrator
          already-migrated self-referential state, with the migrator
          dispatch correctly emitting verdict: PASS, no_op: true,
          re_migration: true. There is no distinct pre-migration
          TestTrack source to diff against — the closest surviving
          pre-migration record is the run log
          `composition-analysis-run.md` (5-step manual flow across 3
          datasets), which is the canonical source I use for content
          checks.
        - D-FRONTMATTER-PHASE1-01 — all 4 Phase 1 fields
          (source_text_fixes, candidate_helpers, unresolved_ambiguities,
          scope_reductions) present as parseable YAML lists.
        - D-FRONTMATTER-PHASE1-02 — all 8 `source_text_fixes[]` slugs
          ≤ 80 chars (max observed: slug 8
          `dataset-folder-substituted-samples-to-tests-per-chain-canonical-bio-fixtures`
          at 75 chars), kebab-case, no duplicates; 1
          `candidate_helpers[]` entry (`bio.flow.composition`) is a
          qualified dotted-path with no args suffix; 3
          `unresolved_ambiguities[]` are kebab-case slugs with no
          duplicates; 1 `scope_reductions[]` entry SR-01 has the
          required {id, check, rationale, verdict_status} shape with
          check: A-CONT-01 and verdict_status: SCOPE_REDUCTION.

      Content checks identify exactly one undocumented matrix

      reduction that requires a `scope_reductions[]` entry but does

      not rise to FAIL:


      D-STRUCT-02 (Scenario 3 dataset-matrix reduction). The canonical

      pre-migration record (`composition-analysis-run.md`) reports all

      5 steps ran across a 3-dataset matrix (filter_FASTA / filter_HELM

      / filter_MSA — implicit 5 × 3 = 15 cells). The migrated scenario

      retains the 3-dataset matrix for Scenarios 1-2 (Steps 1-3,

      `for each dataset in {filter_FASTA.csv, filter_HELM.csv,

      filter_MSA.csv}` loops) but restricts Scenario 3 (Steps 4-5,

      Gear → Context Pane) to filter_FASTA.csv only (body anchor:

      "After Composition is docked in Scenario 1 on

      `filter_FASTA.csv`"). Matrix reduced from 5 × 3 = 15 cells to

      3+3+3+1+1 = 11 cells. The reduction is consistent with chain

      YAML `composition-analysis.md.notes` (filter_FASTA.csv designated

      as representative source per the Step 1 dataset list) and is

      body-visible via the Scenario 3 anchor text. However, no entry

      in `scope_reductions[]` cites `check: D-STRUCT-02` with

      rationale — per the D-STRUCT-02 vocabulary ("Matrix scenario

      combinations lost … without documented reduction reasoning in

      migrated frontmatter's scope_reductions[]"), the reduction is

      undocumented at the frontmatter audit surface.


      This is the same condition surfaced by the prior cycle

      (2026-06-01-bio-migrate-01). The Migrator's default-mode no_op

      this cycle declined to apply the SR-02 proposal (per the

      migrator dispatch's contract_violations_avoided note: applying a

      prior critic's SR-population proposal is the

      `apply_scope_reduction` mode contract, not default-mode

      authoring). The missing SR-02 therefore persists into this cycle

      and the SCOPE_REDUCTION verdict carries forward.


      Proposed remediation (binds via apply_scope_reduction Migrator

      mode or operator-side edit): add SR-02 to `scope_reductions[]`:
        - id: SR-02
          check: D-STRUCT-02
          rationale: |
            Scenario 3 (Steps 4-5, Gear → Context Pane) restricted to
            filter_FASTA.csv rather than the 3-dataset matrix
            exercised by Scenarios 1-2 and the original run log.
            Reasoning: chain YAML composition-analysis.md.notes
            designates filter_FASTA.csv as the representative source
            for this scenario per the Step 1 dataset list. The Gear →
            Context-Pane wiring is a cross-subsystem invariant
            (Bio + Datagrok property-panel framework) that does not
            vary across the 3 Macromolecule notations — exercising it
            on the canonical FASTA branch is sufficient to validate
            the wiring. The click-letter → grid-selection bridge
            (Scenarios 1-2) does exercise all 3 notations to validate
            per-notation canvas hit-handling. Matrix reduced from
            5 × 3 = 15 cells to 3+3+3+1+1 = 11 cells.
          verdict_status: SCOPE_REDUCTION

      All other content checks PASS:
        - D-STEP-01 PASS — all 5 original steps mapped (Steps 1-2 →
          Scenario 1; Step 3 → Scenario 2; Steps 4-5 → Scenario 3).
        - D-STEP-02 PASS — Expected results preserved: WebLogo docks
          (verified against `grok.shell.tv.viewers`), no multi-column
          dialog, ≥1 row selected on letter click, Context Pane opens
          with property surface, property accepts user-driven edit
          (with SR-01 deferring the property-name-specific correctness
          assertion).
        - D-EDGE-01 PASS — Multi-column dialog absent-case explicitly
          noted under Setup (single Macromolecule column per dataset).
          Settle-window (~3s) and x-offset probes (18/30/50/80/120/160)
          from the run log are flagged in `unresolved_ambiguities[]`,
          not silently dropped.
        - D-STRUCT-01 PASS — Cross-scenario relationship to
          `analyze.md` (umbrella runner) documented in body Notes
          block and the chain YAML.
        - D-UI-DELEGATION-01 PASS — SR-01 is an A-CONT-01 deferral
          (property-name correctness), not a UI-to-JS-API
          substitution; the check does not apply.
        - D-SAN-02 PASS — Acknowledgment-contract slugs cover heading
          normalization, dataset-list formatting cleanup, menu-path
          suffix drop ("analysis" trailing word), implicit matrix
          surfacing, Step 2/Step 3 rephrasings, and the sample→filter
          / samples→tests dataset substitutions. Run-log empirical
          observations (settle window, x-offset probes, specific
          property toggles `Show Position Labels` /
          `Skip Empty Positions`) are flagged in
          `unresolved_ambiguities[]`.
        - D-MERIT-01 PASS — SR-01 cites a real atlas dependency
          (Context-Pane property checklist not enumerated in
          `bio.viewers.web-logo` or `bio.analyze.composition` atlas
          entries), not effort/complexity.
        - D-MERIT-02 PASS — SR-01 cites the blocking dependency
          ("until atlas or operator supplies a concrete checklist of
          Context-Pane properties to verify").

      Verdict binds when SR-02 is recorded in `scope_reductions[]`.

      Body prose does not need to change — the reduction is already

      body-visible via the Scenario 3 anchor text and consistent with

      chain-YAML representative-source guidance.
  f:
    verdict: PASS
    cycle_id: 2026-06-01-bio-migrate-02
    timestamp: 2026-06-02T04:15:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-02-bio-automate-01
    timestamp: 2026-06-02T13:25:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-02-bio-automate-01
    timestamp: 2026-06-02T16:02:00Z
    spec_runs:
      - spec: composition-analysis-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 140
        failure_keys: []
---

# Bio | Analyze | Composition — composition analysis integration

Integration scenario for the `Bio | Analyze | Composition` top menu
(atlas `bio.analyze.composition` — `package.ts#L1041`) across three
canonical Macromolecule notations (FASTA, HELM, MSA). Verifies the
multi-subsystem path: top-menu dispatch → Composition WebLogo viewer
docking (`bio.viewers.web-logo`, atlas `package.ts#L449`) →
viewer-to-grid selection bridge (click-letter selects rows) → Gear
icon → Context Pane property editor surface.

The umbrella runner that exercises Composition alongside Sequence
Space and Activity Cliffs lives in `analyze.md`; this scenario is the
deeper Composition-specific scenario. WebLogo paint correctness
(letter-height ordering, color, conservation track) is atlas
`manual_only` (`bio.viewers.web-logo`) — this scenario asserts the
user-visible behaviours that drive that viewer (docking, click→select,
Gear→Context-Pane wiring) without claiming pixel-level WebLogo
correctness.

## Setup

Composition is exercised on each of three test datasets to verify
notation-agnostic behaviour across the canonical Macromolecule
notations:

- `System.AppData/Bio/tests/filter_FASTA.csv`
- `System.AppData/Bio/tests/filter_HELM.csv`
- `System.AppData/Bio/tests/filter_MSA.csv`

The implicit matrix is **3 datasets × 1 Composition flow = 3 cells**.
Each dataset has a single Macromolecule column so the multi-column
choice dialog (atlas `bio.analyze.composition` — "Multi-column choice
dialog when more than one Macromolecule column is present") does NOT
appear in this scenario; the Composition viewer docks directly on
dispatch.

## Scenarios

### Scenario 1 — Open dataset and dispatch Composition

For each dataset in `{filter_FASTA.csv, filter_HELM.csv, filter_MSA.csv}`:

1. Open the dataset from `System.AppData/Bio/tests/`. The Macromolecule
   detector (atlas `bio.detector`) classifies the sequence column
   synchronously; the table view opens.

2. On the menu ribbon, open **Bio** > **Analyze** > **Composition**
   (atlas `bio.analyze.composition`, `compositionAnalysis` —
   `package.ts#L1041`). Verify:
   - A WebLogo viewer docks (atlas `bio.viewers.web-logo`,
     `WebLogoViewer` — `package.ts#L449`); viewer presence verified
     against `grok.shell.tv.viewers` rather than a pixel-level paint
     assertion (paint correctness is atlas `manual_only`).
   - No multi-column choice dialog appears (the dataset has a single
     Macromolecule column).

### Scenario 2 — Click letter in WebLogo selects rows in the grid

For each dataset (same matrix as Scenario 1):

3. Click a letter cell in the docked Composition WebLogo viewer. Verify:
   - At least one row becomes selected in the source grid (viewer →
     grid selection bridge — the `bio.analyze.composition` →
     `bio.viewers.web-logo` click-handler routes the canvas hit to
     a grid-row selection set on the source dataframe).
   - The selection-row count is positive and consistent with the
     letter's frequency at that position (assertion checks "≥1 row
     selected"; the exact count is dataset-position-dependent and
     not atlas-canonical, so the row-count exact match is not
     asserted).

### Scenario 3 — Gear icon opens Context Pane; properties are editable

After Composition is docked in Scenario 1 on `filter_FASTA.csv`:

4. On the docked Composition WebLogo viewer, click the **Gear** icon.
   Verify:
   - The Context Pane opens with the Composition WebLogo viewer's
     property surface bound (the WebLogo `WebLogoViewer` —
     `package.ts#L449` — properties: `sequenceColumnName`,
     `positionWidth`, `startPosition`, `endPosition`, `mode`, etc.).
   - The property-editor surface is present (presence of the
     property-grid is asserted; the concrete property checklist is
     deferred — see SR-01).

5. In the Context Pane, edit at least one editable property of the
   Composition WebLogo viewer. Verify:
   - At least one property accepts a user-driven edit (the
     property-grid input commits the edited value back to the
     viewer; the specific property and its expected post-edit value
     are deferred — see SR-01).
   - No error balloon is surfaced from the property edit (the
     baseline assertion: the edit path does not crash; the property
     correctness assertion is deferred per SR-01).

## Notes

- **Source-text fixes silently applied** during migration (recorded
  in frontmatter `source_text_fixes`):
  - The original heading `#### Composition Analysis manual test` was
    normalized to `# Bio | Analyze | Composition — composition
    analysis integration` to match the sister scenarios' heading
    convention (`analyze.md`, `sequence-space.md`) and to surface
    the atlas top-menu form `Bio | Analyze | Composition` in the
    scenario title.
  - The original Step 1 listed datasets as `System.AppData/Bio/tests/
    filter_FASTA.csv,` / `filter_HELM.csv,` / `filter_MSA.csv.` with
    trailing commas and a trailing period on the last item; the
    migrated Setup uses a clean bullet list.
  - The original Step 2 said "open **Bio** > **Analyze** >
    **Composition** analysis" with a trailing "analysis" word;
    atlas registers the menu item as `Bio | Analyze | Composition`
    (no "analysis" suffix) per `bio.analyze.composition`
    (`package.ts#L1041`). The migrated scenario uses the
    atlas-canonical form.
  - The original Step 1 listed datasets without explicitly framing
    them as a per-dataset loop; surfaced here under Setup as a
    3 × 1 matrix and reflected in each Scenario as a "for each
    dataset" loop.
  - The original Step 2 said "A viewer opens" as a single trailing
    phrase; promoted here to an explicit verification step
    (WebLogo-viewer presence against `grok.shell.tv.viewers` +
    multi-column choice dialog absence).
  - The original Step 3 said "Check that you select elements from
    table if you click any letter in the viewer"; rephrased here as
    the viewer-to-grid selection bridge assertion with an explicit
    "≥1 row selected" baseline.
- **Sub-features covered:**
  - `bio.analyze.composition` (atlas L449, `package.ts#L1041`) —
    top-menu `Bio | Analyze | Composition` dispatch + viewer
    docking, with multi-column dialog when more than one
    Macromolecule column is present (not exercised here).
  - `bio.viewers.web-logo` (atlas L229, `package.ts#L449`) —
    canvas-based sequence-logo viewer; pixel-level paint
    correctness is atlas `manual_only` and is NOT asserted here.
    The user-visible behaviours that drive the viewer (docking
    presence, click-letter-selects-rows bridge, Gear-icon →
    Context-Pane wiring, property-editor edit acceptance) are
    asserted at the assertable surface.
  - Maps directly onto atlas critical path
    `bio.cp.composition-analysis` (p1,
    `derived_from: public/packages/Bio/src/package.ts#L1041`).
- **Manual-only surface touched but not asserted in pixel form:**
  `bio.viewers.web-logo`. Atlas `manual_only` declares the WebLogo
  paint (letter-height ordering, color, conservation track) is
  human-inspection only; this scenario asserts viewer-docking
  presence + click-select bridge + Gear→Context-Pane wiring +
  property-edit acceptance only. Per the chain note
  (`composition-analysis.md` notes block): "this scenario
  exercises the user-visible behaviors that drive that viewer
  without claiming pixel-level WebLogo correctness — assertable
  surface (click select, Gear panel) is automatable even though
  the underlying WebLogo paint is manual-only."
- **No related bugs.** Bug-library `bio.yaml` was consulted and no
  curated bug intersects the Composition / WebLogo-viewer-via-
  Composition surface. GROK-18474 (MSA column-header click crash
  on FASTA, `bio.rendering.column-header` + `bio.viewers.web-logo`)
  is explicitly excluded by chain
  `bug_match_attempts_skipped[GROK-18474]` — this scenario does not
  exercise the MSA column-header click path. `related_bugs: []`
  in frontmatter reflects this.
- **Cross-scenario context:** `analyze.md` exercises Composition
  alongside Sequence Space and Activity Cliffs as an umbrella
  runner (its Scenario 2 also verifies Composition viewer docking;
  its Scenario 3 also exercises the Composition Gear→Context-Pane
  wiring on `filter_FASTA.csv`). `composition-analysis.md` is the
  deeper Composition-specific scenario that additionally asserts
  the click-letter-selects-rows bridge (Scenario 2 here) — the
  unique-to-this-scenario assertion.
- **Unresolved ambiguities** (carried in frontmatter
  `unresolved_ambiguities`):
  - The Context-Pane property checklist for the Composition WebLogo
    viewer's Gear surface is not declared in atlas
    (`bio.viewers.web-logo` describes the WebLogo property surface
    as `sequenceColumnName`, `positionWidth`, `startPosition`,
    `endPosition`, `mode`, etc., but does not pin a canonical
    "which subset must be present in the Context Pane after Gear
    click" checklist). The prior run log
    (`composition-analysis-run.md`) used `Show Position Labels`
    and `Skip Empty Positions` toggles as concrete picks, but
    neither atlas nor an operator decision binds those as the
    canonical edit set; see SR-01.
  - The canvas-interactivity settle window for the Composition
    WebLogo viewer (the time between viewer-docked and
    click-letter-routes-to-grid-selection) is not declared in
    atlas. The prior run log notes a ~3s settle is required on
    the 2nd/3rd dataset before the canvas hit handlers wire
    cleanly — empirical observation, not an atlas contract.
    Operator clarification is needed before Automator can encode
    a typed assertion on the settle window.
  - The click-letter canvas hit-test coordinates are not declared
    in atlas. The prior run log probes a set of x-offsets (`18,
    30, 50, 80, 120, 160`) so narrower letter columns on HELM/MSA
    still land on a letter cell; this is an empirical robustness
    pattern, not an atlas-declared selector. Operator
    clarification is welcomed if a canonical canvas hit-test
    procedure should be encoded.
- **Candidate helpers** (carried in frontmatter `candidate_helpers`):
  - `bio.flow.composition` — registered in
    `helpers-registry.yaml` as a grok-browser procedure for the
    `Bio | Analyze | Composition` flow. Surfaced as a candidate
    Automator-side helper for the top-menu dispatch portion of
    Scenario 1 Step 2.
- **Deferred via SCOPE_REDUCTION** (see frontmatter `scope_reductions`):
  - SR-01 — the correctness assertion on the edited-property state
    in Step 5 (source Step 5, "change arbitrary properties") is
    deferred until atlas or operator supplies a concrete
    Context-Pane property checklist for the Composition WebLogo
    viewer. The Gear→Context-Pane wiring + edit-acceptance flow
    itself is preserved as a step (the Gear icon opens the
    property-grid surface and at least one property accepts a
    user-driven edit).
- **`derived_from:` provenance:** atlas entries cited above derive
  from code anchors only — `bio.analyze.composition` from
  `package.ts#L1041`, `bio.viewers.web-logo` from
  `package.ts#L449`, `bio.cp.composition-analysis` from
  `public/packages/Bio/src/package.ts#L1041`. No help-doc
  derivations were used (per the binding sourcing rule).

---
{
  "order": 10
}
