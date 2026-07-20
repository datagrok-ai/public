---
feature: chem
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: [chem.cp.scaffold-tree-add-filter]
realizes: [chem.analyze.scaffold-tree, chem.scaffold-tree]
realized_as:
  - scaffold-tree-functions-spec.ts
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

This is the section's canonical smoke test, selected because Scaffold Tree is the most heavily used and
most bug-prone viewer in the Chem package.

## Setup

1. **Provision linked dataset.** The scenario consumes a single
   bundled Datagrok file (no external provisioning required):
   - `System:AppData/Chem/tests/smiles-50.csv` — 50-row SMILES dataset
     (canonical SMILES molecule column, the source of the scaffold
     tree).
2. **Confirm Chem package is loaded** so that **Chem > Analyze >
   Scaffold Tree** is registered as a top-menu entry (package source
   `Chem/src/package.ts#L2237`) and the `ScaffoldTreeViewer`
   (`Chem/src/widgets/scaffold-tree.ts`) is registered.
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
   `getScaffoldTree` script).
5. Verify: scaffold tree generation completes — the empty-state
   placeholder is replaced by the rendered scaffold hierarchy
   (scaffold nodes appear in the viewer; the first scaffold node is
   visible and interactable). No console errors fire during
   generation.
6. Click the first scaffold node in the rendered tree.
7. Verify: the `smiles-50.csv` table is filtered by the selected
   scaffold (rows whose molecule contains the scaffold remain;
   non-matching rows are filtered out via the `Scaffold Tree Filter`).
   The grid row count reflects the filter. No console errors.
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

- **Cross-cutting bugs not covered here.** Known Scaffold Tree bugs are covered by dedicated bug-focused
  specs, not this smoke walk: opening a scaffold node in the Sketcher corrupting subsequent
  substructure-search state (GROK-12758), Scaffold Tree blob deserialization on save+reopen with
  Substructure Search (GROK-17595), and a multi-table wrong-binding bug where Scaffold Tree binds to the
  first-opened table instead of the active one (github-3004). This scenario only opens one table and
  never opens a scaffold node in the Sketcher, so none of the three surfaces are exercised here.
- **"Toolbox" clarification.** In this scenario, "the viewer's toolbox" means the in-viewer toolbar at
  the top of the Scaffold Tree viewer, not the application-level Toolbox sidebar.
