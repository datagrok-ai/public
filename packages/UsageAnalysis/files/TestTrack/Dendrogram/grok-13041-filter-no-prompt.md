---
feature: dendrogram
target_layer: playwright
coverage_type: edge
priority: p1
realizes: [GROK-13041]
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
regression guard that pins the post-fix behavior.

## Setup

- A clean Datagrok view (no preloaded tables).
- The Chem package is installed and registered (`Chem | Analyze` top-menu
  surface and `Chem:getMorganFingerprints` available).
- The Dendrogram package is installed and registered.
- Dataset: `System:AppData/Chem/mol1K.csv` (1000-row dataset with a
  `molecule` column carrying SMILES).
- A dendrogram is built and injected as a grid neighbor via
  `Chem | Analyze | Hierarchical Clustering...`. The magic-wand button
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

- See: `bug-library/dendrogram.yaml#GROK-13041` — marked `status:
  fixed`, `priority: p2`.
- See: `.claude/skills/grok-browser/references/dendrogram.md`
  `## bug-grok-13041-filter-no-prompt`.
- **Deferred.** The exact selectors for the remove/revert prompt
  dialog haven't been pinned against the live UI yet; until they are,
  the spec asserts the overlay text instead of the dialog selector
  (see Scenario 2 step 2).
