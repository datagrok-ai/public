---
feature: peptides
target_layer: playwright
coverage_type: edge
priority: p1
realizes_atlas: [GROK-19145]
realizes: []
produced_from: atlas-driven
related_bugs:
  - GROK-19145
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - sar-similarity-threshold-matrix-spec.ts
gate_verdicts:
  b:
    verdict: PASS
    cycle_id: 2026-05-30-peptides-automate-02
    timestamp: 2026-05-31T00:30:00Z
    spec_runs:
      - spec: sar-similarity-threshold-matrix-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 112
        failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-05-30-peptides-automate-02
    timestamp: 2026-05-31T00:00:00Z
    failure_keys: []
  f:
    verdict: PASS
    cycle_id: 2026-05-29-peptides-migrate-02
    timestamp: 2026-05-29T00:00:00Z
    failure_keys: []
---

# Peptides — SAR completes cleanly across Similarity threshold values, including the edge case that crashed before

Launching SAR with the Similarity threshold set to low, medium, high, and an extreme value should never crash with a null-reference error, even when the extreme value produces very sparse or no qualifying pairs. This guards against GROK-19145, where a high Similarity threshold could throw on a null BitSet; on sparse output the expected behavior is either a rendered (possibly sparse) result or an explicit "no qualifying pairs" message — never a silent blank WebLogo with no error.

## Setup

1. Open the linked Peptides demo dataset (`System:DemoFiles/bio/peptides.csv`) — a TableView with the peptides grid and a `Macromolecule`-semtype `AlignedSequence` column should be the active view.
2. Confirm the active table has at least one numerical activity column (the demo dataset's `Activity` column satisfies the `peptidesDialog` validator gate — "active table has a Macromolecule column and at least one numerical activity column"); no warning toast should appear when SAR is launched.

## Scenarios

### Scenario 1 — Similarity threshold matrix: SAR completes without a crash across low / medium / high / extreme values

Runs the SAR launch four times with increasingly aggressive Similarity threshold values, including `90` — the value that originally triggered the GROK-19145 crash. At every value, SAR must either complete and show the configured viewers, or surface an explicit informative message — never an uncaught null-reference crash.

1. Launch the SAR configuration dialog via top-menu `Bio | Analyze | SAR...`. The `Analyze Peptides` dialog opens (`peptidesDialog` → `analyzePeptidesUI`).
2. Set the **Similarity** threshold to `10` (low value — sparse filter, most pairwise comparisons qualify) leaving all other inputs at their defaults (Activity column auto-detected, Scaling `-lg`).
3. Click **OK**. SAR analysis runs to completion; the `Sequence Variability Map` viewer attaches to the TableView and renders with populated cells; the WebLogo column-headers on per-position columns render monomer stacks (non-empty WebLogo).
4. Reopen the dialog via `Bio | Analyze | SAR...`, set **Similarity** to `50` (medium value), click **OK**. SAR re-launches without crash; viewers render; WebLogo column-headers render monomer stacks.
5. Reopen the dialog via `Bio | Analyze | SAR...`, set **Similarity** to `75` (high value — aggressive filter, fewer pairwise comparisons qualify), click **OK**. SAR re-launches without crash; viewers render; WebLogo column-headers render (may show fewer monomers per stack if stats are sparse, but must not be silently empty without an accompanying empty-state indicator).
6. Reopen the dialog via `Bio | Analyze | SAR...`, set **Similarity** to `90` (extreme value — the original GROK-19145 repro value), click **OK**. SAR analysis must NOT throw an uncaught `method not found: 'setTrue' on null` runtime error in the JS console / error balloon; the workflow either completes (viewers present, WebLogo header may show sparse but rendered output), or surfaces an explicit informative balloon / dialog message about the filter being too aggressive. The Datagrok error-balloon mechanism must not fire on a null-receiver TypeError.

Expected (assertion summary across all four Similarity values):
- No uncaught error balloon, no `Cannot read properties of null` / `method not found ... on null` console error, no thrown TypeError from `setTrue` / `getRawData` / similar BitSet operations.
- For each Similarity value, either the `Sequence Variability Map` viewer is attached AND WebLogo column-headers are rendered (populated or sparsely populated but visibly drawn), OR an informative platform message is surfaced explaining that the filter produced no qualifying pairs.
- The TableView remains responsive after each parametrized launch (no frozen UI, no lingering modal blocking subsequent dialog opens).

### Scenario 2 — Similarity threshold = 90 preserves the WebLogo rendering contract even with sparse output

Isolates the "WebLogo silently renders blank on sparse stats" failure mode described in the bug report, using a single minimal-viewer launch so the sparse-output risk isn't diluted by the broader matrix in Scenario 1.

1. Launch SAR via top-menu `Bio | Analyze | SAR...`, set **Similarity** to `90`, leave Activity / Scaling at defaults, and additionally toggle off any optional viewers (Most Potent Residues, Logo Summary Table, etc.) in the dialog so the launch produces only the minimal Sequence Variability Map configuration — minimizing the surface area on which to assert.
2. Click **OK** and wait for SAR launch to settle.
3. Verify: the `Sequence Variability Map` viewer is present in the TableView (attached, not skipped).
4. Verify: the per-position column headers on the peptides grid render a WebLogo glyph (`peptides.rendering.weblogo-header` is installed) — either populated with monomer stacks OR drawn as an explicit empty-state placeholder. A header that is silently blank (no glyph, no placeholder, no error) is a regression.
5. Verify: clicking a monomer on any populated WebLogo header invokes the selection backbone (`modifySelection` → `fireBitsetChanged` → DataFrame BitSet update) without crashing — the WebLogo header click handler remains live even when the underlying stats are sparse.

Expected:
- No null-receiver runtime errors regardless of how sparse the underlying `MonomerPositionStats` is at Similarity = 90.
- `Sequence Variability Map` viewer attaches; WebLogo column-headers draw either populated stacks or an explicit empty-state placeholder.
- WebLogo click handler responds without crashing (selection backbone survives sparse stats).

## Notes

- **Related bug — GROK-19145.** Both scenarios exist to reproduce and guard against this bug (a high Similarity threshold could throw on a null BitSet). This `.md` is the primary anchor scenario for that bug; a dedicated regression spec extends the same threshold-parameter surface.
- **See:** `public/help/datagrok/solutions/domains/bio/peptides-sar.md#Analysis options`; `public/help/datagrok/solutions/domains/bio/peptides-sar.md#Table View`.

## Original trailing metadata

```json
{
  "order": 6,
  "datasets": [
    "System:DemoFiles/bio/peptides.csv"
  ]
}
```
