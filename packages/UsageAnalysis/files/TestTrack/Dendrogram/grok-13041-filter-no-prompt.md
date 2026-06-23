---
feature: dendrogram
sub_features_covered:
  - dendrogram.clustering.inject-tree-for-grid
  - dendrogram.event.selection-changed
target_layer: playwright
coverage_type: edge
produced_from: atlas-driven
related_bugs:
  - GROK-13041
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - grok-13041-filter-no-prompt-spec.ts
gate_verdicts:
  e:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T20:15:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T19:45:30Z
    spec_runs:
      - spec: grok-13041-filter-no-prompt-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 73
        failure_keys: []
---

# Dendrogram — GROK-13041 filter does not trigger remove/revert prompt

Bug-focused regression for [GROK-13041](https://reddata.atlassian.net/browse/GROK-13041)
("Bio | Tools: Dendrogram needs reset after filtering"). With a dendrogram
attached as a grid neighbor, applying a column filter to the host grid MUST
NOT trigger the remove-dendrogram / revert-sort prompt that is intended only
for sort-driven row reordering. The bug is fixed; this scenario is the
regression guard that pins the post-fix behavior. Anchors atlas edge case
`dendrogram.ec.filter-does-not-trigger-remove-revert-prompt`
(`source_bug: GROK-13041`, `coverage_type: regression`).

## Setup

- A clean Datagrok view (no preloaded tables).
- The Chem package is installed and registered (`Chem | Analyze` top-menu
  surface and `Chem:getMorganFingerprints` available).
- The Dendrogram package is installed and registered.
- Dataset: `System:AppData/Chem/mol1K.csv` (1000-row dataset with a
  `molecule` column carrying SMILES).
- A dendrogram is built and injected as a grid neighbor via
  `Chem | Analyze | Hierarchical Clustering...` (the canonical entry point
  for `dendrogram.clustering.inject-tree-for-grid`). The magic-wand button
  `.dendrogram-assign-clusters-bttn` is the mounted-and-ready signal.

## Scenarios

### Scenario 1 — Filtering a column does NOT trigger the remove/revert prompt

Steps:

1. Open `System:AppData/Chem/mol1K.csv` and wait for the `molecule` column
   to render structures.
2. Run `Chem | Analyze | Hierarchical Clustering...`. In the
   `[name="dialog-Hierarchical-Clustering"]` dialog set `Features` =
   `molecule`, `Distance` = `euclidean`, `Linkage` = `ward`, click
   `[name="button-OK"]`. Wait for the dendrogram neighbor to appear (poll
   for `.dendrogram-assign-clusters-bttn`).
3. Confirm the dendrogram neighbor is attached:
   `!!document.querySelector('.dendrogram-assign-clusters-bttn')` is
   `true`, and the renderer registers itself on the grid via
   `grok.shell.tv.grid.temp['__dendrogram_neighbor_temp__']`.
4. Open the Filter Panel for the active table view
   (`grok.shell.tv.getFiltersGroup()`).
5. Apply a column filter on a numeric column (e.g. set the `pIC50_HIV_Integrase`
   filter to a range that excludes a non-empty subset of rows). The filter
   restricts the visible row set without changing the underlying row order.

   * Expected result: no `.d4-dialog` carrying a "remove dendrogram" or
     "revert sort" prompt appears. No console error is emitted. The
     dendrogram neighbor remains attached to the grid
     (`.dendrogram-assign-clusters-bttn` still in DOM;
     `grok.shell.tv.grid.temp['__dendrogram_neighbor_temp__']` still set).
     No `Revert columns sort order to see Dendrogram Tree` overlay is
     shown over the tree neighbor.
6. Clear the filter.

   * Expected result: dendrogram remains attached and unchanged; no
     additional dialog or prompt was raised during the filter / clear
     round-trip.

### Scenario 2 — Sorting a column DOES trigger the prompt (positive contrast)

Establishes the symmetric positive contrast — the trigger condition is
genuinely guarded by the prompt for the sort path, ensuring Scenario 1's
absence assertion is not vacuously trivial.

Steps:

1. Continuing from the Setup state (mol1K loaded, dendrogram attached,
   filter cleared).
2. Click the header of a numeric column (e.g. `pIC50_HIV_Integrase`) to
   sort it ascending. Row order in the grid changes; the dendrogram's
   leaf order no longer matches the grid's row order.

   * Expected result: the remove-dendrogram / revert-sort prompt is
     raised, and the `Revert columns sort order to see Dendrogram Tree`
     overlay text is shown over the tree neighbor
     (`inject-tree-for-grid2.ts#L326`). One of two outcomes is acceptable:
     either a `.d4-dialog` containing the remove/revert choice text
     appears, OR the overlay is shown and the tree is hidden until the
     sort is reverted. The exact selectors for the prompt dialog are not
     yet captured (ref doc section `## bug-grok-13041-filter-no-prompt`
     flags them as a follow-up MCP recon item). The spec author SHOULD
     pin both the prompt dialog selector and the overlay text node in a
     follow-up recon round; until then, the spec MAY assert the overlay
     by its known text content via
     `getHTMLElementbyInnerText('div', 'Revert columns sort order')`
     (or the equivalent UITests helper) and skip the prompt-dialog
     selector.
3. Click the sorted column header again to clear the sort, or revert via
   the overlay action.

   * Expected result: dendrogram is reattached and the overlay is
     dismissed.

## Notes

- target_layer rationale: filter-vs-sort prompt-triggering is a UI
  surface (Filter Panel input + grid header click + prompt dialog +
  overlay text), not assertable through the JS API alone. `playwright`
  is the canonical layer per the STEP D heuristic (multi-step UI flow
  with cross-dialog state).
- coverage_type rationale: `edge` per the dispatch bug brief and the
  Critic-F SR proposed_action — the bug brief's `coverage_type:
  regression` is overridden to `edge` to also discharge
  F-STRUCT-NEGATIVE-01 (the section has no `coverage_type: edge` or
  `coverage_type: perf` scenario; this one closes that gap as the
  Critic-F SR explicitly recommends). The atlas edge case itself is
  recorded as `coverage_type: regression` in
  `feature-atlas/dendrogram.yaml#edge_cases`, but the scenario's
  surface here is the boundary-condition trigger discrimination (filter
  vs sort), which is the canonical `edge` shape.
- related_bugs: [GROK-13041] is load-bearing — it is the signal Critic
  F's F-BUG-COVERAGE-01 branch (ii) re-checks on the next round to
  clear the `bug-uncovered` gap (anchored on the realized .md, not on
  the chain `bug_focused_candidates[]` proposal).
- Net-new vs atlas-resolved coverage union: 0 breadth (both `affects`
  ids are already in the union via `assign-clusters.md` and
  `hierarchical-clustering-bio.md`). Value here is non-breadth per the
  Test Designer net-new-differencing rule's escape clause for
  justified-non-breadth scenarios — F-BUG-COVERAGE-01 closure +
  F-STRUCT-NEGATIVE-01 closure are the value axes.
- See: `feature-atlas/dendrogram.yaml#edge_cases[dendrogram.ec.filter-does-not-trigger-remove-revert-prompt]`
  (`source_bug: GROK-13041`, `coverage_type: regression`,
  `derived_from: bug-library:dendrogram.yaml#GROK-13041`).
- See: `bug-library/dendrogram.yaml#GROK-13041` (`status: fixed`,
  `priority: p2`, `test_coverage: needed`).
- See: `.claude/skills/grok-browser/references/dendrogram.md`
  `## bug-grok-13041-filter-no-prompt` (DOM-validated 2026-06-03;
  sort-prompt selectors flagged as follow-up MCP recon).
- atlas entry derived from public/packages/Dendrogram/src/viewers/inject-tree-for-grid2.ts#L326
  (sort-overlay registration site; `Revert columns sort order to see
  Dendrogram Tree` text).
- Deferrals: the sort-prompt-dialog selectors are explicitly deferred
  per the ref doc's Forward TODO — pinning them requires a follow-up
  MCP recon round against `dev.datagrok.ai` (real dependency:
  dialog-selectors not in the 2026-06-03 recon's selector validation
  matrix). Per Lattice Rule 13 / A-MERIT-02 this is a cited real
  technical dependency, not a TODO-add-later.
