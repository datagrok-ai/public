---
feature: charts
target_layer: playwright
coverage_type: edge
priority: p1
realizes: [tree-rowsource-onclick-state-machine]
pyramid_layer: bug-focused
ui_coverage_responsibility: []
ui_coverage_delegated_to: null
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Charts/tree-rowsource-onclick-state-bug.md
date_created: 2026-05-07
authored_by: orchestrator-test-designer-charts-migrate-2026-05-07
related_bugs:
  - github-3245
---

# Tree viewer — Row Source × On Click state machine (github-3245)

Bug-focused regression scenario for github-3245: the Tree viewer's
`rowSource` × `onClick` combination used to behave non-orthogonally —
the same logical state could produce different visual results
depending on the order transitions happened in. The fix landed in
Charts 1.24.0; this scenario locks in state-machine consistency across
all 9 combinations of `rowSource` (All/Filtered/Selected) × `onClick`
(Select/Filter/None), regardless of the path taken to reach them. The
canonical `tree.md` only exercises collaborative filtering at default
settings, not these combinations.

## Setup

A clean Datagrok session. Single scenario; demog.csv as the source.
Programmatic exercise — no canvas Shift+Click (canvas synthesis
AMBIGUOUS per tree.md Notes precedent).

## Scenarios

### Scenario 1: Tree rowSource × onClick combinations are state-consistent (github-3245 invariant)

Steps:

1. Open `System:DemoFiles/demog.csv` and add a Tree viewer with
   hierarchy `['CONTROL', 'SEX', 'RACE']`. Wait 3000ms.

2. **Iterate over 9 combinations** of `rowSource` ∈ {`All`, `Filtered`,
   `Selected`} × `onClick` ∈ {`Select`, `Filter`, `None`}:
   - For each combination, call
     `tree.setOptions({rowSource: <rs>, onClick: <oc>})`.
   - Wait 800ms for state to commit.
   - Read back `tree.props.get('rowSource')` and `tree.props.get('onClick')`
     (race-tolerant try/catch).
   - **Invariant 1 (state-consistency):** if read-back succeeds, the
     values match the setOptions input.
   - **Invariant 2 (visual stability):** viewer's root DOM is
     non-empty and has non-zero size (no white screen / DOM tear-down).
   - **Invariant 3 (no console error):** no new console error during
     the transition (filtering benign 404 / favicon noise).
   - Log `{rs, oc, readbackRs, readbackOc, hasContent}` for diagnostic
     visibility.

3. **Transition history independence:** after iterating 9 combinations
   in one order (e.g., All/Select → All/Filter → ... → Selected/None),
   re-iterate in REVERSE order and verify all 3 invariants still hold.
   This is the "transition history" check — the same logical state must
   produce the same visual result regardless of path taken.

4. **Final ESC-equivalent state-clear check:** call
   `df.selection.setAll(false)` + `df.selection.fireChanged()` to
   simulate ESC clearing table selection. Verify the Tree viewer's
   root DOM remains visually stable (no stuck highlights — programmatic
   equivalent of the bug ticket's "ESC does not clear viewer highlight"
   regression check).

## Notes

- The github-3245 invariant: all 9 `rowSource` × `onClick`
  combinations satisfy state-consistency, visual stability, and
  no-console-error — and the same holds when the combinations are
  iterated in reverse order (transition-history independence).
- The original bug ticket describes pressing ESC to clear the table
  selection while the Tree viewer's highlight stays stuck. This
  scenario uses the programmatic equivalent (`df.selection.setAll(false)`
  + `fireChanged()`) as a substitute — an actual DOM ESC keypress
  isn't exercised yet.

## Dataset metadata

```json
{
  "order": 35,
  "datasets": ["System:DemoFiles/demog.csv"]
}
```
