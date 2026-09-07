---
feature: caching
realizes_atlas: [caching.cp.computed-cell-render-cache, caching.cp.query-result-cache]
priority: p0
target_layer: manual-only
coverage_type: smoke
manual_only_reason: |
  Both checks are perceived-latency judgments (first computation vs cached
  re-render) with no DOM state to assert on.
related_bugs: []
---

# Local caching — manual checklist

## Cached rendering of computed cell widgets

Use the full SPGI (not spgi-100): the cache contrast needs enough rows for the
first computation pass to be noticeably slow.

1. Close all, open SPGI
2. Add a computed widget column: open **Add New Column** and enter the formula
   `Chem:ChemistryGasteigerPartialCharges(${Structure})`, name it `gasteiger`
3. Scroll the grid down through unvisited rows — new cells compute visibly
   (placeholders while the depiction is being calculated)
4. Scroll back up to already-visited rows — those cells render immediately
   from cache, with no recomputation delay
5. Scroll down again over the previously computed range — same instant render
6. Cleanup: remove the `gasteiger` column and close the view

## Cached query results

Use a Chembl query not yet run in this session — a pre-warmed server cache
would mask the effect.

1. In Browse, open Databases > Chembl and run **ChemblNumberOfStructures** —
   note the time to results
2. Run **BrowseAllChEMBLStructures**, then switch back to
   **ChemblNumberOfStructures** — the previously requested results appear
   almost instantly (no full re-execution)
3. Switch between the two queries a few times — every revisit renders from
   cache with no visible wait

---
{
  "order": 1,
  "datasets": ["System:DemoFiles/chem/SPGI.csv"]
}
