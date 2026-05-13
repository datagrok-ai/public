# Migration Report — info-panels.md

First-cycle migration of `info-panels.md` (TestTrack/Chem section) per
chain YAML `scenario-chains/chem.yaml` rev 1
(`output_plan.info-panels.md`: target_layer `playwright`, strategy
`simple`, classification `medium`, `pyramid_layer: integration`,
`depends_on: []`).

## Step mapping

The original is a two-block scenario separated by `***`. Block 1 (5
numbered steps + 4 nested sub-steps under step 4 + 2 implicit
verification bullets under step 5) walks the column-header Context
Panel surface on `smiles.csv`, drives Rendering + Highlight panels on
`chembl-scaffolds.csv`, then walks the cell-level Context Panel on the
first molecule of `smiles.csv`. Block 2 (6 numbered steps) is a multi-
format coverage walk across SMILES / molV2000 / molV3000 / SMARTS,
exercising cell rendering and the Context Panel tab set per format.
JSON footer carries `order: 1` + `datasets[]` (5 entries) — metadata,
not a step; preserved in chain YAML `order_from_files` and migrated
`## Setup` step 1.

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| Block 1 Step 1 — "Open linked datasets" | Setup step 1 (5 datasets enumerated) + Phase A step 1 (open `smiles.csv`) + Phase A step 5 (open `chembl-scaffolds.csv`) + Phase B per-row step 1 (open per-format dataset) | preserved (split for clarity — original "open linked datasets" expands across Setup + Phase A + Phase B opens) |
| Block 1 Step 2 — "Use smiles.cvs" | Phase A step 1 (open `smiles.csv`) + Phase A step 2 (verify SMILES auto-detection) | preserved (split for clarity — typo "cvs" → "csv" silently corrected) |
| Block 1 Step 3 — "Click the canonical_smiles column's header." | Phase A step 3 | preserved |
| Block 1 Step 4 — "Go to the Context Pane and check all info panels (tabs)." | Phase A step 4 (walk all Chemistry / Biology / Structure tabs) | preserved |
| Block 1 Step 4 sub-step 1 — "for testing panel Chemistry -> Rendering use chembl_scaffolds.cvs" | Phase A step 5 (open `chembl-scaffolds.csv`) + step 6 (open Rendering panel on it) | preserved |
| Block 1 Step 4 sub-step 2 — "Choose 'Scaffold' as a Scaffold column and check 'Highlight scaffold'." | Phase A step 6 (Scaffold column choice + Highlight scaffold checkbox) | preserved |
| Block 1 Step 4 sub-step 3 — "Molecules in 'Smiles' column should be aligned by scaffold and scaffold should be highlighted" | Phase A step 7 (verification) | preserved as verification |
| Block 1 Step 4 sub-step 4 — "Chemistry -> Hightlight. Add the structure to Hightlight and choose the color different from Scaffold hightlightnings. Both of the hightlightnings should be visible on Smiles column." | Phase A step 8 (add structure to Highlight at distinct color) + Phase A step 9 (verify both layered highlights visible) | preserved as verification (split for clarity — action + verification separated; original typos "Hightlight"/"hightlightnings" silently corrected to "Highlight"/"highlights") |
| Block 1 Step 5 — "Click the first molecule in the canonical_smiles column." | Phase A step 10 (return to `smiles.csv` + click first molecule cell) | preserved |
| Block 1 Step 5 bullet 1 — "Context Panel switches to molecule data." | Phase A step 11 (verify Context Panel mode-switch column→molecule) | preserved as verification |
| Block 1 Step 5 bullet 2 — "Expand all tabs on molecule Context Panel. No errors should appear." | Phase A step 12 (expand every molecule-context tab, no errors) | preserved as verification |
| Block 2 Step 1 — "Open linked datasets (need to check smiles / molV2000 / molV3000 / smarts formats)" | Phase B (data-driven table over 4 formats) + per-row step 1 (open dataset for current format) | preserved (split for clarity — multi-format scope made explicit as a 4-row parameterized data-driven walk) |
| Block 2 Step 2 — "Make sure that structures of molecules are rendered." | Phase B per-row step 2 (verify cell rendering for format) | preserved as verification |
| Block 2 Step 3 — "Click a cell with a structure." | Phase B per-row step 3 (click molecule cell) | preserved |
| Block 2 Step 4 — "Check that all necessary panels are displayed on the Context Panel." | Phase B per-row step 4 (verify Chem panels appear) | preserved as verification |
| Block 2 Step 5 — "Expand each tab on the Context Panel." | Phase B per-row step 5 (expand each tab) | preserved |
| Block 2 Step 6 — "Make sure the content for each info panel is displayed correctly." | Phase B per-row step 6 (verify content renders) | preserved as verification |
| JSON footer `{ "order": 1, "datasets": [...5 entries...] }` | (dropped from body) | metadata-not-step (chain analysis convention; captured in `scenario-chains/chem.yaml` `order_from_files` + migrated Setup step 1 dataset enumeration) |

No original numbered step or implicit verification bullet is silently
dropped.

## Decisions

- **Why this `target_layer`:** chose `playwright` per
  `scenario-chains/chem.yaml` `output_plan.info-panels.md.target_layer
  = playwright`. The scenario requires DOM-level Context Panel mode-
  switching (column-header click vs. cell click), tab expansion, and
  visual layering of two highlight types on canvas-rendered cells —
  all UI behaviors that the playwright layer drives natively and the
  apitest layer cannot exercise. The chain YAML reason field:
  "Standalone Context Panel walk across 4 format variants and
  Rendering+Highlight tabs. Sibling info-panels-spec.ts already exists
  at playwright layer per existing-test-index (category: Chem). The
  Rendering panel scaffold-alignment and Highlight color-layering flows
  require real DOM context-panel rendering that playwright provides
  natively."
- **Why this `coverage_type`:** chose `regression`. The scenario is not
  a single-bug repro (no `related_bugs` populated; Phase A + Phase B
  walk the canonical Chem Context Panel surface across 4 formats), not
  a one-tab smoke (12 Chem panel ids exercised across Phase A + B), and
  not a perf test (no timing or volume thresholds). It IS a multi-
  format regression walk over the Chem panel surface — verifying that
  every Chem info panel renders correctly for every supported molecule
  notation format. Regression matches: protects against feature
  surface erosion as Chem panels are added / refactored / extended.
- **Why this `strategy`:** `simple` per chain YAML
  `output_plan.info-panels.md.strategy = simple`. Single scenario, no
  cross-file fixture; chain analyzer classified it `medium`, not
  `complex-standalone`. Pattern 1 (decomposition) is NOT triggered.
  Phase B uses a data-driven parameterized loop (the format table in
  the migrated `### Phase B` block — single sequence of 6 steps
  iterated over 4 format rows) — that is a within-scenario shape
  choice for downstream Automator to render as a parameterized
  Playwright `test.describe.parallel(...).forEach(...)` or equivalent;
  it does NOT change the scenario's `strategy: simple` classification.
- **Sibling tests consulted (READ-ONLY per Invariant 2):**
  - `public/packages/UsageAnalysis/files/TestTrack/Chem/info-panels-spec.ts`
    — existing playwright-layer test ("Chem: Info Panels on smiles.csv",
    category Chem, helpers_called `spec-login`, patterns
    `uses-fixture-runners` / `uses-grok.dapi` / `uses-grok.shell` /
    `uses-page.evaluate` / `uses-page.locator` / `uses-playwright-test`).
    Per `existing-test-index.yaml` line 32499-32514. Confirms the
    section convention: playwright + `spec-login` helper + page.evaluate
    for JS API setup + page.locator for DOM assertions. Read-only;
    not modified per Migrator boundary.
  - `public/packages/UsageAnalysis/files/TestTrack/Peptides/info-panels-spec.ts`
    — adjacent section "Peptides — Info Panels" — same shape pattern;
    cross-section parity confirms the info-panels naming convention.
    Read-only.
  - `public/packages/UsageAnalysis/files/TestTrack/Projects/complex.md`
    + `complex-migration-report.md` — shape anchor for migrated `.md`
    + report structure (frontmatter, Setup/Scenarios/Notes order, Step
    mapping table, Decisions / Opt-outs / Deferred / Edge cases /
    Unresolved ambiguities). Read-only.
- **Helpers reused / candidate helpers:**
  - **Reused (registered in `helpers-registry.yaml`):**
    - `loginToDatagrok` (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:49`)
      — section-standard fixture login; Automator will use in spec
      `beforeEach`.
    - `softStep` (`public/packages/UsageAnalysis/files/TestTrack/spec-login.ts:18`)
      — per-step soft-fail wrapper; Automator wraps each Phase A / Phase
      B verification in a `softStep`.
    - `expandAccordionPane`
      (`public/packages/UITests/playwright/scripts/helpers.ts:141`) —
      Context Panel tabs are accordion panes; this helper expands a pane
      by title (e.g. `expandAccordionPane(page, 'Descriptors')`).
      Applicable in Phase A step 4 + Phase A step 12 + Phase B step 5.
    - `closeAllViews`
      (`public/packages/UITests/playwright/scripts/helpers.ts:129`) —
      cleanup between Phase A and Phase B + between Phase B format
      iterations.
  - **Candidate helpers (NOT yet in registry — flagged for addition
    via helpers-registry curator; per migration-prompt Helpers
    discipline §):**
    - `helpers.playwright.chem.openMoleculeColumnContextPanel(page,
      dataset, columnName)` — opens dataset + clicks the molecule
      column header to drive Context Panel into column-context mode.
      Surfaced by Phase A steps 1-3 and Phase B per-row steps 1-3
      (cell variant). Not invented in the migrated body; spec author
      can inline equivalent code until helper lands.
    - `helpers.playwright.chem.openMoleculeCellContextPanel(page,
      dataset, columnName, rowIndex)` — opens dataset + clicks the
      first molecule cell to drive Context Panel into molecule-context
      mode. Surfaced by Phase A steps 10-11 and Phase B per-row
      steps 3-4. Companion helper to the above.
    - `helpers.playwright.chem.assertChemPanelTabsPresent(page,
      expectedTabSet)` — asserts that the Chem-authored tab set
      (parameter `expectedTabSet`) is present on the Context Panel,
      tolerant of additional non-Chem tabs injected by other packages.
      Surfaced by Phase A step 4 and Phase B per-row step 4. Encodes
      the "tab inventory caveat" from migrated Notes.
- **Bug library consulted:** yes — `bug-library/chem.yaml` rev 1
  (12 curated_bugs, generated_date 2026-05-05). Three bugs intersect
  this scenario's flows but do NOT trigger inclusion in `related_bugs`
  per the chain YAML's frontmatter delta — all three are
  parallel-coverage candidates owned by dedicated bug-focused specs in
  the chain's `bug_focused_candidates[]`:
  - **GROK-16870** (RDKit cell renderer crashes on null DataFrame in
    non-Chem viewer tooltip context). `affects` intersects
    `chem.rendering.rdkit-renderer`, `chem.rendering.molecule-cell`,
    `chem.rendering` — all three appear in this scenario's
    `sub_features_covered`. However the bug's repro requires hovering
    over a data point in a non-Chem viewer (Box Plot, Scatter Plot,
    Histogram, Grid) — info-panels.md does NOT add such viewers; it
    walks the Context Panel ON the grid view, which is a different
    rendering context. Chain emits `chem-grok-16870-spec.ts` as the
    bug-focused candidate; awareness only here. Cited in migrated
    Notes section.
  - **GROK-17964** (Convert Notations duplicate after handler error,
    persists across project reload). `affects` `chem.notation.action`,
    `chem.notation`, `chem.notation.convert-mol`. info-panels.md does
    NOT trigger a Convert Notations error; chain emits
    `chem-grok-17964-spec.ts` as the bug-focused candidate. Awareness
    only here. Cited in migrated Notes section.
  - **GROK-17621** (Calculate Descriptors Docker timeout) was reviewed
    and rejected for inclusion — this scenario triggers the
    **Chemistry > Descriptors** Context Panel tab (single-molecule
    OCL-based panel, `chem.panels.descriptors`), NOT the top-menu
    `Chem | Calculate | Descriptors` workflow which is the Docker-
    container-backed surface that this bug affects
    (`chem.calculate.descriptors`). Different sub_features; no
    overlap.
- **Decision log queried:** yes — `decision-log.yaml` (10049 lines,
  revision tracked) grepped for `feature: chem` — **zero matches**.
  This is the FIRST CYCLE for the chem section. No prior
  `migration_decisions`, `layer_decisions`, `manual_only`, or
  `failed_attempts` entries exist for chem. No "approaches off the
  table" constraints apply.
- **Cross-cutting bug citations (chain YAML `bug_focused_candidates[]`).**
  Per `migration-prompt.md` "Cross-cutting bug citations from chain
  YAML" §: this scenario is a span in two candidates — GROK-16870
  (proposed_spec `chem-grok-16870-spec.ts`) and GROK-17964
  (proposed_spec `chem-grok-17964-spec.ts`). Both are cited in
  migrated Notes section ("Cross-cutting bug awareness — GROK-16870"
  + "Cross-cutting bug awareness — GROK-17964"). Per migration-prompt
  this is RECOMMENDED, not mandatory; F-BUG-COVERAGE-01 at section-
  complete is authoritative.
- **UI delegation status.** Per chain YAML
  `ui_coverage_plan.delegated_scenarios.info-panels.md`:
  `delegated_to: null` (UI coverage owned, not delegated to
  Advanced/scaffold-tree-functions.md ui-smoke). All
  `ui_coverage_responsibility` flows (12 panel UI flows +
  `context-panel-column-header-click` + `context-panel-cell-click`)
  are owned by this scenario. No UI coverage gap introduced; no
  SCOPE_REDUCTION proposal moves UI driving to a sibling — therefore
  the UI delegation discipline of migration-prompt § "UI delegation
  in SCOPE_REDUCTION proposals" does not apply here (no SR proposal
  to delegate).
- **No invented sub_features / helpers.** Every id in
  `sub_features_covered` is present in
  `feature-atlas/chem.yaml` (verified per atlas lines 60-80 for
  rendering + 891-892 for notation.detect-smiles + 1088-1171 for
  panels). Every reused helper is in `helpers-registry.yaml`;
  candidate helpers are flagged as candidates, NOT referenced by name
  in the migrated body.

## Opt-outs (SCOPE_REDUCTION proposals)

### SR-01: A-STRUCT-02 carryforward (chain-level edge/perf coverage)

**Cited technical dependency:** A-STRUCT-02 ("at least one scenario is
edge case or negative path — i.e. its frontmatter `coverage_type` is
`edge` or `perf`") is a **section-level structural invariant**. The
chem `info-panels.md` scenario covers the Context Panel rendering
surface across multiple notation formats (smiles / molV2000 / molV3000
/ SMARTS) as a regression walk — it exercises the multi-format
invariant but is not itself an edge-case probe or a negative-path test
of a specific failure mode. Forcing `coverage_type: edge` onto
`info-panels.md` would mis-classify it; forcing a synthetic Phase C
with `coverage_type: edge` would invent content not present in the
original (D-SAN-02 violation risk).

**Resolution path** (mirrors radar.md 2026-05-07 carryforward
precedent, decision-log entry around line 8528-8530 — "radar.md:
SCOPE_REDUCTION (3 waivers — A-STRUCT-02 + A-STRUCT-05 carryforward +
A-MERIT-02 step 13 new)"): A-STRUCT-02 will be satisfied **chain-wide**
by the 10 `bug_focused_candidates[]` already emitted under
`scenario-chains/chem.yaml` rev 2 (each chem-bug-focused spec
naturally carries `coverage_type: edge` because it tests a specific
bug repro path — sister-pattern to charts `radar-save-reopen-bug.md`,
`sunburst-date-column-bug.md`, `sunburst-scatterplot-color-pollution-
bug.md`, `tree-rowsource-onclick-state-bug.md`, all of which carry
`coverage_type: edge` per `grep -rE "coverage_type:"
public/packages/UsageAnalysis/files/TestTrack/Charts/*.md` 2026-05-11
section inventory).

**Why this is a real dependency, not effort:** the per-scenario
edge|perf marker cannot be invented mid-migration without violating
D-SAN-02 (no silent additions). The chain-wide satisfaction path is
already wired (10 candidate specs to be authored downstream by Test
Designer / Automator). info-panels.md must not block the migrate
cycle on a structural property that is the **chain's** responsibility,
not the per-scenario responsibility.

**Scope of waiver:** A-STRUCT-02 verdict on `info-panels.md` is
**deferred to chain-level evaluation**. Critic A should return
SCOPE_REDUCTION (carryforward), not FAIL, citing this proposal.

_All other content checks (A-STRUCT-01, A-STRUCT-03 through
A-STRUCT-06, A-COVERAGE-*, A-MERIT-*) are expected to PASS on
`info-panels.md` without any other SR proposals._

## Deferred items (NOT opt-outs)

(none)

_No step is deferred awaiting a prerequisite that does not exist yet.
The migrated body's steps are realizable against the current
playwright + helpers-registry + atlas state. (The candidate helpers
surfaced in Decisions are convenience abstractions that the
Automator may inline pending registration — they do NOT block spec
realization.)_

## Edge cases

The original lists no explicit "edge case" keyword. Implicit edge
cases derivable from the scenario content:

- **SMARTS format edge case.** Block 2's `smarts` format slice
  (`SMARTS_example_temp.csv`) exercises the SMARTS-only branch of the
  RDKit cell renderer + the SMARTS-only branch of the Context Panel
  tab content. SMARTS cells cannot round-trip through SMILES-only
  parsing — `chem.notation.is-smarts` regex-detection path applies.
  PRESERVED as scenario step (Phase B per-row step 2 verification for
  the `smarts` row, plus implicit-format-edge-case note in the
  migrated step 2 body).
- **molV2000 vs. molV3000 format dichotomy.** Block 2 walks both
  `molV2000` (mol1K.sdf) and `molV3000` (ApprovedDrugs2015.sdf) — the
  RDKit cell renderer + MolfileHandler must handle both V2000 and V3K
  blocks. The V3K branch is regression-prone (V3K is the newer block-
  format extension, exercised less in legacy datasets). PRESERVED as
  scenario step (Phase B rows molV2000 and molV3000, per-row step 2
  rendering verification).
- **Context Panel mode-switch (column-context vs. molecule-context).**
  Phase A step 3 clicks the column header (column-context); Phase A
  step 10 clicks a cell (molecule-context). The Context Panel must
  switch modes cleanly — different tab sets visible per mode
  (column-context shows Rendering / Highlight panel which apply to
  the column; molecule-context shows Descriptors / Drug Likeness /
  etc. which apply to the single molecule). PRESERVED as scenario
  steps (Phase A step 3 + Phase A step 11 verification — explicit
  mode-switch assertion).
- **Two-layer highlight stacking on canvas-rendered cells.** Phase A
  step 8-9 layers a custom-structure Highlight panel highlight ON TOP
  OF the Scaffold-alignment Rendering panel highlight. Both must
  render at distinct colors simultaneously without clobbering — this
  is the layering invariant. PRESERVED as scenario step (Phase A
  step 9 verification).
- **Tab inventory caveat (cross-package tab injection).** Atlas-
  derived: other packages (Bio, Peptides, ChemDraw) may inject
  additional Context Panel tabs not authored by Chem. Migrated body's
  Phase A step 4 / Phase B step 4 verify "all applicable Chem info
  panels" — bounded to the Chem-authored set per atlas
  `chem.panels.*` (12 panel ids), tolerant of additional non-Chem
  tabs. PRESERVED as a Notes-level caveat in migrated body (not a
  separate scenario step — the caveat governs how the verification
  is interpreted, not what is performed). Flag for atlas curator: a
  future atlas extension could enumerate "Chem-authored tabs only"
  as a sub_feature `chem.panels.set` for cleaner machine reading.
- **SMILES auto-detection threshold (≥80% via `detectSmiles`).** Atlas
  `chem.notation.detect-smiles` documents `detectSmiles(col, min)` —
  "populates units/semType/cellRenderer when ≥80% of sampled
  categories parse as SMILES." `smiles.csv` (Phase A) cleanly meets
  the threshold; molV2000/molV3000 SDF files arrive with explicit
  units (block format, not SMILES) so the detection path differs.
  PRESERVED as an implicit verification (Phase A step 2 — verify
  auto-detection populated units / semType / cellRenderer). Format-
  specific detection paths are not flagged as a separate scenario
  step (over-specification for the regression intent); the
  verification is bounded to "molecule cells render" not "specific
  notation-detection telemetry."

No edge case is moved to atlas, manual_only, deferred, or a separate
scenario. All implicit edge cases are PRESERVED as scenario steps OR
Notes-level caveats in the migrated body.

## Unresolved ambiguities

- **Exact tab inventory across Phase A step 4 and Phase A step 12.**
  Step 4 (column-header click) and step 12 (cell click) drive
  different Context Panel modes — column-context exposes Rendering +
  Highlight panels (which apply to the column as a whole), and
  molecule-context exposes single-molecule panels (Descriptors / Drug
  Likeness / etc.). The original lists "all info panels (tabs)" once
  in step 4 without distinguishing modes. Migrated body resolves by
  walking the same Chem panel id set (11-12 panels per atlas
  `chem.panels.*`) at step 4 AND at step 12 — but the Automator at
  spec time should pin the exact tab set per mode by reading
  `chem.panels.*` parents at runtime, NOT by hardcoding a tab name
  list. Flag for QA pair review at spec realization.
- **"Choose 'Scaffold' as a Scaffold column" — Scaffold column name
  case-sensitivity.** Original references `Scaffold` (capital S) as
  the column name in `chembl-scaffolds.csv`. The actual column name
  in the file may be `Scaffold` exactly OR a variant
  (`scaffold` / `scaffolds`). Migrated body preserves the original
  casing `Scaffold`. Automator verifies at spec time against the
  actual dataset column header (resilient match preferred over
  case-sensitive hardcode).
- **"a structure" for the Highlight panel (Phase A step 8).** Original
  says "Add the structure to Hightlight" — unclear which structure.
  Most likely: some substructure SMARTS / SMILES that is distinct
  from the Scaffold column's scaffold (so the two highlights overlay
  visibly at distinct colors). Migrated body says "add a structure to
  the Highlight list and choose a color distinct from the Scaffold
  highlight color" — preserves the ambiguity (any non-Scaffold
  structure is acceptable). Automator picks a deterministic structure
  at spec time (e.g. a small SMARTS `c1ccccc1` for the aromatic ring
  pattern, distinct from the Scaffold column entries).
- **"Color different from Scaffold hightlightnings" — concrete color
  pair.** Original asks for "the color different from Scaffold" —
  unclear which specific colors. The Highlight panel's color picker
  is the binding mechanism (color picker widget). Migrated body
  preserves the ambiguity ("a color distinct from the Scaffold
  highlight color"). Automator picks deterministic colors at spec
  time (e.g. Scaffold red `#FF0000` + Highlight blue `#0000FF`).
- **Number of cells to verify in Phase B step 2 ("Make sure structures
  are rendered").** Original is silent on "verify the first cell?
  all visible cells? a sample?". Most performant: verify the first N
  visible cells render. Migrated body says "structures of molecules
  appear in cells via the RDKit cell renderer" — bounded by
  visibility. Automator picks N at spec time (e.g. first 5 visible
  cells).
- **Block 1 vs. Block 2 dataset reuse for the `smiles` format.** Block 2
  step 1 explicitly enumerates `smiles (smiles.csv)` as one of the 4
  formats — but Block 1 already walked the Context Panel on
  `smiles.csv` (steps 1-5). Whether Block 2's SMILES format slice is a
  re-walk of Block 1's coverage OR a distinct walk (e.g. a different
  cell, different tabs) is silent in the original. Migrated body treats
  the SMILES Phase B row as a per-format-parity walk over the same
  cell-context Chem panels — symmetric with the molV2000 / molV3000 /
  smarts rows. Automator may de-duplicate Phase A coverage at spec
  realization (e.g. mark Phase B SMILES row as `.skip` if Phase A
  coverage is canonical for SMILES), OR keep it as a parity assertion;
  spec-time decision.
