---
feature: chem
sub_features_covered: [chem.analyze.scaffold-tree, chem.analyze.scaffold-tree.viewer, chem.analyze.scaffold-tree.add, chem.analyze.scaffold-tree.generate, chem.analyze.scaffold-tree.filter]
target_layer: playwright
coverage_type: smoke
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/Advanced/scaffold-tree-functions.md
migration_date: 2026-05-11
source_text_fixes: []
candidate_helpers:
  - helpers.playwright.chem.addScaffoldTreeViewer
  - helpers.playwright.chem.generateScaffoldTreeViaMagicWand
  - helpers.playwright.chem.clickScaffoldTreeNode
unresolved_ambiguities:
  - check-them-step-5-property-panel-verification-depth
  - magic-wand-icon-dom-locator
  - scaffold-tree-viewer-empty-state-placeholder-text
scope_reductions:
  - id: SR-01
    check: A-STRUCT-02
    rationale: "A-STRUCT-02 carryforward (chain-level edge/perf coverage)"
    verdict_status: SCOPE_REDUCTION
related_bugs: []
---

# Scaffold Tree functions

Canonical Scaffold Tree viewer create + configure smoke walk. Four
beats end-to-end on `smiles-50.csv`:

- **Add** the Scaffold Tree viewer via **Top menu > Chem > Analyze >
  Scaffold Tree** — empty-tree state opens.
- **Generate** the scaffold tree via the magic-wand icon (tooltip
  *Generate from molecular column*).
- **Click** the first scaffold node — the table filters by that
  scaffold.
- **Inspect** the viewer toolbox (top of viewer) and the viewer
  property panel (Context Panel).

Per chain YAML (`scenario-chains/chem.yaml` rev 2): independent
scenario (`depends_on: []`), `classification: simple`,
`pyramid_layer: ui-smoke`, `target_layer: playwright`, strategy
`simple`. Selected as the section's UI smoke per Rule 1 — canonical
viewer-create-and-configure pattern for Scaffold Tree, the chain's
most-tested viewer (3 of 10 atlas critical paths and 4 of 12 known
bugs touch it). Realizes atlas critical path
`chem.cp.scaffold-tree-add-filter` (p0) end-to-end. UI coverage
owned (`ui_coverage_delegated_to: null`).

## Setup

1. **Provision linked dataset.** The scenario consumes a single
   bundled Datagrok file (no external provisioning required):
   - `System:AppData/Chem/tests/smiles-50.csv` — 50-row SMILES dataset
     (canonical SMILES molecule column, the source of the scaffold
     tree).
2. **Confirm Chem package is loaded** so that **Chem > Analyze >
   Scaffold Tree** is registered as a top-menu entry (per atlas
   `chem.analyze.scaffold-tree.add`, package source
   `Chem/src/package.ts#L2237`) and the `ScaffoldTreeViewer`
   (`Chem/src/widgets/scaffold-tree.ts`, atlas
   `chem.analyze.scaffold-tree.viewer`) is registered.
3. **No fixture consumed.** Per chain YAML `depends_on: []`. The
   dataset is opened fresh inside the scenario.

## Scenarios

### Scaffold Tree add + generate + filter + inspect

End-to-end walk: open dataset → add viewer (empty state) → magic-wand
generate → first-node click filter → toolbox + property-panel
inspection.

1. Open `System:AppData/Chem/tests/smiles-50.csv`. Verify the molecule
   column renders via the RDKit cell renderer and the SMILES column
   auto-detects as a Molecule column (semType / units / cell renderer
   populated). The table view becomes the active view.
2. From the top menu, run **Chem > Analyze > Scaffold Tree**.
3. Verify: the Scaffold Tree viewer is added to the active table view
   in its empty state — the viewer displays a placeholder reading
   **Scaffold Tree is empty** (no nodes rendered yet, since the tree
   has not been generated).
4. In the empty Scaffold Tree viewer, press the **magic-wand** icon
   (tooltip: **Generate from molecular column**). Scaffold tree
   generation starts (the Scaffold Tree blob is produced via the
   `getScaffoldTree` script — atlas
   `chem.analyze.scaffold-tree.generate`).
5. Verify: scaffold tree generation completes — the empty-state
   placeholder is replaced by the rendered scaffold hierarchy
   (scaffold nodes appear in the viewer; the first scaffold node is
   visible and interactable). No console errors fire during
   generation.
6. Click the first scaffold node in the rendered tree.
7. Verify: the `smiles-50.csv` table is filtered by the selected
   scaffold (rows whose molecule contains the scaffold remain;
   non-matching rows are filtered out via the `Scaffold Tree Filter`
   — atlas `chem.analyze.scaffold-tree.filter`). The grid row count
   reflects the filter. No console errors.
8. With the Scaffold Tree viewer focused, inspect the viewer
   **toolbox** at the top of the viewer (the in-viewer toolbar
   surface — icons / actions exposed by the Scaffold Tree viewer for
   in-viewer operations). Verify: the toolbox renders without errors
   and its action set is visible.
9. With the Scaffold Tree viewer still focused, open the viewer's
   property panel on the **Context Panel** (e.g. via the viewer's
   gear icon, or by clicking the viewer header so that the Context
   Panel switches to viewer-properties mode). Verify: the property
   panel renders the Scaffold Tree viewer's property set
   (configurable knobs of the viewer) without errors. Expand the
   property groups and confirm each property control is visible and
   interactable.

## Notes

- **`coverage_type: smoke`** — `pyramid_layer: ui-smoke` per chain
  rev 2; the scenario is the canonical Scaffold Tree create + configure
  smoke walk. `regression` / `edge` / `perf` classifications do not
  apply: there is no parameter-axis multi-format walk (single dataset,
  single format), no specific bug-repro path being locked in (the
  scenario walks the happy path only), and no timing or volume
  thresholds. The section-level A-STRUCT-02 edge|perf coverage path
  is satisfied chain-wide via the 10 `bug_focused_candidates[]` in
  `scenario-chains/chem.yaml` rev 2 (see migration report SR-01).
- **No JS API substitution.** Every step in this scenario is a
  declared UI flow per chain `ui_coverage_responsibility`
  (`chem-add-scaffold-tree`, `chem-scaffold-tree-magic-wand-generate`,
  `chem-scaffold-tree-node-click-filter`, `chem-scaffold-tree-toolbox`,
  `chem-scaffold-tree-property-panel`). Because the scenario is the
  section's UI smoke (`pyramid_layer: ui-smoke`, Rule 1 election with
  documented step-count deviation cited in chain
  `unresolved_ambiguities`), UI driving via Playwright is mandatory.
  Substituting `grok.dapi.*` for the top-menu add, the magic-wand
  click, the scaffold-node click, the toolbox inspection, or the
  property-panel inspection would defeat the smoke's intent and
  violate the chain's UI-coverage-responsibility contract.
- **Cross-cutting bug awareness — GROK-12758.** Per chain
  `bug_focused_candidates[]`, the Scaffold Tree × Sketcher ×
  Substructure-Search integration bug (opening a scaffold node in
  the Sketcher corrupts subsequent `searchSubstructure` state) is
  owned by the dedicated `chem-grok-12758-spec.ts` candidate spec
  (spans `Advanced/scaffold-tree.md:Step 1`, `filter-panel.md:Step
  1`, `sketcher.md:Step 2`). `scaffold-tree-functions.md` does NOT
  open scaffold nodes in the Sketcher — it walks add → magic-wand
  generate → node-click filter → toolbox + property panel only.
  Awareness only; `related_bugs: []` because the bug-focused candidate
  owns the invariant.
- **Cross-cutting bug awareness — GROK-17595.** Per chain
  `bug_focused_candidates[]`, the Scaffold Tree blob deserialization
  bug on save+reopen with Substructure Search + datasync is owned
  by the dedicated `chem-grok-17595-spec.ts` candidate spec.
  `scaffold-tree-functions.md` does NOT save / close / reopen a
  project; the auto-generated tree blob is in-session only per chain
  `produces` field. Awareness only; not added to `related_bugs`.
- **Cross-cutting bug awareness — github-3004.** Per chain
  `bug_focused_candidates[]`, the multi-table Scaffold Tree
  wrong-binding bug (opening Scaffold Tree from the second table
  binds to the first-opened table) is owned by the dedicated
  `chem-github-3004-spec.ts` candidate spec (spans
  `Advanced/scaffold-tree.md:Step 1` and
  `Advanced/scaffold-tree-functions.md:Step 1`).
  `scaffold-tree-functions.md` opens only one table
  (`smiles-50.csv`); the multi-table state invariant is not
  exercised here. Awareness only; not added to `related_bugs` —
  the bug-focused candidate owns the invariant.
- **Helpers usage.** No registered helper currently abstracts the
  Scaffold Tree create + magic-wand + node-click pattern. The spec
  Author (Automator) will use existing
  `helpers-registry.yaml` entries — `loginToDatagrok`, `softStep`,
  `closeAllViews` for the standard section harness, and may inline
  the Scaffold-Tree-specific DOM driving until a dedicated helper
  lands. Candidate helpers surfaced in migration report Decisions
  section.
- **Order in chain.** `order: 1` per source JSON; tie-broken
  lexicographically AFTER `info-panels.md`. Section-level smoke
  scenario per chain `ui_coverage_plan.smoke_scenario`.
- **Magic-wand icon discovery.** The magic-wand icon is the
  in-viewer affordance for **Generate from molecular column**
  (tooltip text verbatim from the original scenario). Discovery is
  via tooltip text inside the empty Scaffold Tree viewer; Automator
  locates by tooltip role/text rather than by a fragile class
  selector.
- **Toolbox surface.** The "toolbox of the scaffold tree viewer" in
  the original refers to the in-viewer toolbar (top of the viewer),
  not the application-level **Toolbox** sidebar
  (`.d4-toolbox`). The scaffold tree viewer exposes a localized
  action surface that should render with the viewer; the verification
  is "is the action set visible and free of rendering errors", not
  an enumeration of every icon (the icon set is implementation-defined
  in `Chem/src/widgets/scaffold-tree.ts` and may evolve).
