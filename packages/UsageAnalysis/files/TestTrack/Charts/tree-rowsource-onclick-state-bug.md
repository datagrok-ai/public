---
feature: charts
sub_features_covered:
  - charts.tree
  - charts.tree.on-click
  - charts.tree.color-palette
target_layer: playwright
coverage_type: edge
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

Bug-focused regression scenario for github-3245: Tree's `rowSource` ×
`onClick` state machine exhibits non-orthogonal behavior — the same
logical state can produce different visual results depending on
transition history. Fix landed in Charts 1.24.0 — this scenario locks
the state-machine consistency invariant.

`pyramid_layer: bug-focused`. The canonical `tree.md` does NOT exercise
the rowSource × onClick combinations — it only does collaborative
filtering at default settings. This scenario covers each of the 9
combinations (3 rowSource × 3 onClick) for state consistency.

`related_bugs: [github-3245]` — bug-library reproduction class:
combinations of `rowSource` (All/Filtered/Selected) × `onClick`
(Select/Filter/None) must produce consistent viewer state regardless of
transition path.

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

- **github-3245 invariant carrier:** all 9 rowSource × onClick
  combinations satisfy the 3 invariants (state-consistency, visual
  stability, no console error) AND the transition-history independence
  check (reverse-order iteration).
- **ESC behavior:** the bug ticket cites ESC clearing table selection
  but viewer highlight remaining stuck. Programmatic equivalent
  (`df.selection.setAll(false)` + `fireChanged()`) is the
  selector-pending substitute. Real DOM ESC awaits
  `references/charts.md` UI registry per cycle deferral.
- **Authority:** atlas-driven; closes the bug coverage gap for
  github-3245 surfaced in chain rev 2.

## Dataset metadata

```json
{
  "order": 35,
  "datasets": ["System:DemoFiles/demog.csv"]
}
```
