# Migration Report — color-consistency

## Step mapping

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| 1. "Open SPGI" | Setup step 1 | preserved |
| 2. "Add six viewers: Histogram, Line chart, Bar chart, Pie chart, Trellis plot, Box plot (omit Scatter plot — it's covered by scatterplot.md)" | Setup step 2 | preserved |
| 3. "Set the categorical legend on each viewer to Stereo Category, using the viewer's own legend-source property" | Setup step 3 | preserved |
| 4. "In the grid, enable Categorical color coding for Stereo Category, then change at least two category colors — every other viewer's legend must reflect the new colors" | Scenarios > 1 steps 1–3 | preserved (split for clarity) |
| 5. "On any single viewer (e.g. Bar chart), open the legend color picker for one category and pick a new color — the change must propagate back to the grid and to every other viewer (this verifies the column is the single source of truth, not the viewer)" | Scenarios > 2 steps 1–4 | preserved (split for clarity) |
| 6. "Save layout → re-apply layout → verify every viewer still shows the customized palette" | Scenarios > 3 steps 1–3 | preserved |
| 7. "Save project → reopen project → verify the customized palette persists across the round-trip" | Scenarios > 4 steps 1–3 | preserved |

## Decisions

- **Why this `target_layer`:** chose `playwright` per chain `output_plan` — color picker dialog widget + grid color-coding setup + save-layout / reapply / save-project / reopen persistence cycles require Playwright UI driving.
- **Why this `strategy`:** `simple` per chain — single linear flow, no chained sub-state across scenarios.
- **Why this `pyramid_layer`:** `integration` per chain Rule 4 — multi-viewer co-existence (6 viewer types verifying bidirectional sync grid↔viewers↔legend). Not smoke (delegates standard UI to visibility-and-positioning.md), not source-matrix, not bug-focused.
- **Why this `priority`:** `p0` — the "column is single source of truth" invariant is foundational; per-viewer customization workflows break if it regresses. Customer signal NIBR (github-3132) reinforces priority.
- **Why this `coverage_type`:** `edge` (post-SR application 2026-05-07; was initially `regression`). Critic A returned SCOPE_REDUCTION at gate A round 1 noting that Scenarios 3 + 4 substantively test the GROK-17278 persistence edge case from atlas `edge_cases[]`. The file's persistence-across-save invariant (across BOTH layout and project save paths) is the more specific characterization than generic regression. Sections 1 + 2 retain regression-class semantics within the file but the dominant atlas-anchored characterization is edge. SR applied without re-invoking Critic A per orchestrator routing rule (SR does not count as review round).
- **Sibling tests consulted:**
  - `visibility-and-positioning-spec.ts` (section UI smoke) — color picker dialog pattern + helper imports.
  - `color-consistency-spec.ts` (existing) — same scenario; informs spec-style anchoring for downstream Automator rewrite.
- **Helpers reused:** `loginToDatagrok`, `softStep`, `specTestOptions`, `stepErrors` from `spec-login` (registered). No new helpers proposed.
- **Bug library consulted:** yes — 3 bugs intersect sub_features_covered:
  - GROK-17438 affects `legend.item.color-picker`, `legend.visibility`, `legend.splitter-resize` — color-picker intersection. Cross-cutting bug-focused spec proposed.
  - github-3132 affects `legend.item.color-picker`, `legend.allow-item-coloring` — direct match. Sequential color-change additivity is the bug; this scenario tests one-color-change baseline.
  - GROK-17278 affects `legend.item.color-picker`, `legend.column` — color-picker intersection. Layout/project persistence is the bug surface; Scenarios 3 + 4 test positive baseline.
- **Decision log queried:** yes — empty for `legend` feature.
- **UI delegation:** standard legend UI flows (visibility / position / save-dialog widgets) delegated to `visibility-and-positioning.md`. Color-picker dialog and grid color-coding setup are owned by this scenario (not delegated).
- **Cross-cutting bug citations from chain YAML:** chain `bug_focused_candidates[]` lists GROK-17438 (Step 5), github-3132 (Step 4), GROK-17278 (Step 6) with spans into this scenario. Awareness only; spec-level downstream.
- **Scatter plot deliberately omitted from setup:** original step 2 explicitly excludes Scatter plot ("omit Scatter plot — it's covered by `scatterplot.md`"). Preserved verbatim — section's color-consistency story splits across these two scenarios.

## Opt-outs (SCOPE_REDUCTION proposals)

(none)

## Deferred items (NOT opt-outs)

(none)

## Edge cases

- **Single-source-of-truth invariant:** preserved as Scenario 2 step 4 ("verify the change propagates to every other viewer's legend... column is single source of truth, not a per-viewer state"). This is the bug-class invariant for github-3132 and GROK-17438.
- **Persistence across BOTH layout AND project save:** two separate scenarios (3 and 4). These are independent serialization paths — both must round-trip color customizations. Bug GROK-17278 documents the regression on linechart specifically; this scenario verifies the positive baseline across all 6 host viewers.
- **Categorical color coding mode requirement:** Scenario 1 step 1 explicitly enables `Categorical color coding` in the grid before changing colors. Without categorical mode, color customizations apply differently (per the atlas `legend.use-custom-color-coding` sub-feature). Edge handled by enforcing the mode in the setup.

## Unresolved ambiguities

- **Scenario 2 step 1 — "any viewer (e.g. Bar chart)":** original step says "On any single viewer (e.g. Bar chart)" — Bar chart is the example, but the test could legitimately use any of the 6 viewers. Migrated picks Bar chart deterministically; rationale: Bar chart's stack semantics (Scenario 9 in filtering.md) and explicit citation in source make it the canonical choice. Flag for downstream Automator if a different viewer turns out to have flake characteristics in the color picker hover.
