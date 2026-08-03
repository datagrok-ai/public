---
feature: peptides
sub_features_covered:
  - peptides.workflow.start-analysis
  - peptides.util.get-selection-bitset
  - peptides.compute.calculate-monomer-position-statistics
  - peptides.rendering.weblogo-header
  - peptides.workflow.sar-dialog
target_layer: playwright
coverage_type: edge
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

## Setup

1. Open the linked Peptides demo dataset (`System:DemoFiles/bio/peptides.csv`) — a TableView with the peptides grid and a `Macromolecule`-semtype `AlignedSequence` column should be the active view.
2. Confirm the active table has at least one numerical activity column (the demo dataset's `Activity` column satisfies the `peptidesDialog` validator gate — "active table has a Macromolecule column and at least one numerical activity column"); no warning toast should appear when SAR is launched.

## Scenarios

### Scenario 1 — Similarity threshold matrix: SAR completes without null-receiver crash across low / medium / high / extreme values

Atlas anchor: `edge_cases[0]` (`coverage_type: edge`, `source_bug: GROK-19145`). Exercises the post-OK compute path with the Similarity threshold ranging across representative values, including the extreme value (`90`) that triggered the original `setTrue` on null-BitSet crash. The scenario asserts graceful behavior: SAR either completes and shows the configured viewers, or surfaces an explicit informative error — but never throws a null-receiver runtime crash.

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

### Scenario 2 — Similarity threshold = 90 on a sparser-input dataset preserves WebLogo rendering contract

Atlas anchor: same `edge_cases[0]` entry; isolates the WebLogo-empty-on-sparse-stats failure mode described in the bug body. Uses a single-launch path against the dataset configuration most likely to trip the empty-output branch, so the empty-WebLogo regression risk surfaces independently of the multi-launch matrix.

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

- **Atlas provenance.** Both scenarios derive from atlas `edge_cases[0]` (`feature-atlas/peptides.yaml` edge_cases block, `source_bug: GROK-19145`, `coverage_type: edge`, `derived_from: bug-library/peptides.yaml#GROK-19145`). The `sub_features_covered:` frontmatter list mirrors the atlas entry's `sub_features:` verbatim, preserving the canonical mapping required by Critic E spec-mode cross-checks.
- **target_layer rationale.** Multi-step UI dialog flow (open `Analyze Peptides` dialog, set Similarity input, click OK, observe TableView mutation, assert against the rendered Sequence Variability Map + WebLogo headers, optionally inspect Datagrok error-balloon surface). The Similarity threshold is a dialog input that must be driven through the UI to faithfully reproduce the post-OK compute-path crash described in the bug — `apitest` cannot exercise the same code path because the dialog's `analyzePeptidesUI` → `startAnalysis` callback wraps the same code with UI-bound validation and toast routing. The failure mode is also asserted against the UI-bound `peptides.rendering.weblogo-header` renderer and the Datagrok error-balloon surface, both of which are UI-only.
- **coverage_type rationale.** Atlas `edge_cases[0].coverage_type: edge` is canonical per the rule "When the scenario being authored maps onto an atlas `edge_cases[]` entry, the frontmatter `coverage_type:` MUST match that entry's `coverage_type:` verbatim. Do NOT re-derive." Scenario carries `coverage_type: edge` to match the atlas entry. This scenario simultaneously retires the section-wide `F-STRUCT-NEGATIVE-01` gap (no edge / perf coverage_type scenario across the section) flagged in `scenario-chains/peptides.yaml :: gate_f_verdict.gaps[type: negative-coverage-missing]`.
- **Coverage contribution.** Sub_features added to the section's union: `peptides.util.get-selection-bitset` (previously uncovered, on the `util.*` uncovered family from the Gate F gaps), and `peptides.compute.calculate-monomer-position-statistics` (previously uncovered, on the `compute.*` uncovered family). `peptides.workflow.start-analysis`, `peptides.rendering.weblogo-header`, and `peptides.workflow.sar-dialog` were already in the section's union via `sar.md` / `peptides.md` — re-coverage here is intentional (the edge-class anchoring is on the same workflow shape but exercises the post-OK compute path that the smoke / regression scenarios do not).
- **Related bugs.** `GROK-19145` is the single curated bug whose `affects[]` matches this scenario's `sub_features_covered[]` verbatim per `bug-library/peptides.yaml#GROK-19145`. The chain's existing `bug_focused_candidates[0]` (`peptides-grok-19145-spec.ts`) proposes a dedicated cross-scenario regression spec covering the same threshold-parameter surface; this `.md` is the migrated-scenario anchor those candidates extend.
- **Deferrals.** None. Scenario steps map to existing Playwright-driveable surfaces (top-menu invocation, dialog interaction, viewer attachment assertion, WebLogo header rendering assertion, error-balloon absence assertion).
- **See:** `public/help/datagrok/solutions/domains/bio/peptides-sar.md#Analysis options` (atlas help_docs rich-object form maps the `Analysis options` heading to `peptides.workflow.analyze-ui` + `peptides.workflow.sar-dialog`); `public/help/datagrok/solutions/domains/bio/peptides-sar.md#Table View` (maps to `peptides.workflow.start-analysis`).

## Original trailing metadata

```json
{
  "order": 6,
  "datasets": [
    "System:DemoFiles/bio/peptides.csv"
  ]
}
```
