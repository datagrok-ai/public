# Migration Report — line-chart

## Step mapping

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| 1. "Open SPGI" | Setup step 1 | preserved |
| 2. "Add a Line chart" | Setup step 2 | preserved |
| 3. "Set Split = Series — the legend must list every distinct value of Series (typically 4–10 categories), each rendered with a distinct color from the categorical palette (no two adjacent categories share the same color)" | Scenarios > 1 steps 1–3 | preserved (split for clarity) |
| 4. "On the Context Panel > Data, enable Multi Axis (lc.props.multiAxis = true) — each Y line gets its own subplot, and each subplot has its own legend driven by Split" | Scenarios > 2 steps 1–3 | preserved (split for clarity) |
| 5. "Save layout → re-apply layout — Split, Multi Axis, and per-line legends must survive the round-trip" | Scenarios > 3 steps 1–3 | preserved (split for clarity) |
| 6. "Save project → reopen project — same state must persist" | Scenarios > 3 steps 4–5 | preserved (split for clarity) |
| 7. "Configure two Y columns: lc.props.yColumnNames = ['Average Mass', 'TPSA'] — both lines render, each with its own legend block" | Scenarios > 4 steps 1–2 | preserved (split for clarity) |
| 8. "Replace one Y column via the in-plot column selector (or lc.props.yColumnNames = ['Average Mass', 'NIBR logP'] if the in-plot selector is hover-only) — the corresponding legend block must update to reflect the new column's values" | Scenarios > 4 steps 3–4 | preserved (split for clarity) |
| 9. "Save layout → re-apply layout — the new Y column choice must persist" | Scenarios > 4 steps 5–6 | preserved (split for clarity) |
| 10. "Save project → reopen project — same state must persist" | Scenarios > 4 steps 7–8 | preserved (split for clarity) |

## Decisions

- **Why this `target_layer`:** chose `playwright` per chain `output_plan` — multiple save-layout / reapply / save-project / reopen persistence cycles + Context Panel UI driving for Multi Axis toggle + in-plot column selector. Combination requires Playwright.
- **Why this `strategy`:** `simple` per chain — single linear flow building incrementally on the same Line chart; no independent test sections.
- **Why this `pyramid_layer`:** `integration` per chain Rule 4 — multi-subsystem within Line chart (Split + Multi Axis + multi-Y configuration + persistence).
- **Why this `priority`:** `p1` — regression coverage for line-chart-specific legend behaviors. Not p0 because line chart is one of seven viewer types and the section ui-smoke is the canonical p0.
- **Why this `coverage_type`:** `regression` at file level — protects against line-chart legend regressions across Multi Axis dispatch and Y-column replacement paths. **Scenarios 3 + 4 carry per-scenario `[coverage_type: edge]` heading markers post-SR application 2026-05-07** (Critic A returned SCOPE_REDUCTION at gate A round 1: persistence-cycle scenarios 3 + 4 directly exercise GROK-17278 atlas edge_case pattern; per-scenario heading markers satisfy A-STRUCT-02 without splitting the file). SR applied without re-invoking Critic A per orchestrator routing rule (SR does not count as review round).
- **Sibling tests consulted:**
  - `visibility-and-positioning-spec.ts` (section UI smoke) — helper imports + spec shape.
  - `line-chart-spec.ts` (existing) — same scenario; informs spec-style anchoring for downstream Automator rewrite.
  - `public/packages/UsageAnalysis/files/TestTrack/Viewers/line-chart-spec.ts` (sibling Line chart test outside Legend section) — broader line chart coverage; not Legend-specific.
- **Helpers reused:** `loginToDatagrok`, `softStep`, `specTestOptions`, `stepErrors` from `spec-login` (registered). No new helpers proposed.
- **Bug library consulted:** yes — 3 bugs intersect sub_features_covered:
  - GROK-17222 (refresh.on-data-change / column) — direct title match ("Line chart: legend is not consistent with filtering"). Bug repro requires a filter step not in this scenario; bug-focused spec `legend-grok-17222-spec.ts` covers the repro.
  - GROK-17278 (color-picker / column) — direct title match ("Linechart: changed colors are not saved to the layout and project"). Bug repro changes a category color via legend not in this scenario; bug-focused spec `legend-grok-17278-spec.ts` covers the repro.
  - GROK-19083 (show-main-item-icons / column) — column intersection only; markers not relevant to Line chart. Weak match retained per migration prompt rule (intersection on any sub-feature).
- **Decision log queried:** yes — empty for `legend` feature.
- **UI delegation:** standard legend UI flows (visibility / position / color picker / save-dialog widgets) delegated to `visibility-and-positioning.md`. The Context Panel > Data toggle for Multi Axis and the in-plot column selector are owned by this scenario.
- **Cross-cutting bug citations from chain YAML:** chain `bug_focused_candidates[]` lists GROK-17222 (line-chart.md:Step 3) and GROK-17278 (line-chart.md:Step 5); spans map to migrated Scenario 1 step 1 and Scenario 3 step 1 respectively. Awareness only.

## Opt-outs (SCOPE_REDUCTION proposals)

(none)

## Deferred items (NOT opt-outs)

(none)

## Edge cases

- **Layout filter-panel rebuild settle (≥3 s):** preserved as Scenario 3 step 2 ("allow ≥3 s settle"). Layout-restore assertions are flake-prone without this wait.
- **In-plot column selector vs `yColumnNames` API path:** preserved as Scenario 4 step 3 with both paths — original explicitly notes the in-plot selector is hover-only (the JS API is the deterministic alternative). Implicit edge: the test should assert the intended UI path (in-plot selector) when possible; falls back to JS API if hover gating defeats Playwright.
- **Per-subplot legend cardinality with Multi Axis:** preserved as Scenario 2 step 3 — Multi Axis introduces N subplot legends (one per Y column), each driven independently by `Split`. Implicit edge: legend count must equal Y column count after Multi Axis enable.
- **Categorical palette no-adjacent-collision:** preserved as Scenario 1 step 3 — original explicitly asserts "no two adjacent categories share the same color". Implicit edge protecting against palette regression.

## Unresolved ambiguities

- **Scenario 4 step 3 — in-plot column selector vs `yColumnNames` JS API path:** original says "Replace one Y column via the in-plot column selector (or `lc.props.yColumnNames=[...]` if the in-plot selector is hover-only)". Migrated step retains both paths; downstream Automator decides based on whether the in-plot selector is automatable in Playwright (hover gating may defeat selectors). If selector is hover-only and not automatable, JS API path is canonical and a B14 reference-addition proposal may be needed for the in-plot selector DOM. Flag for downstream.
- **Scenario 4 step 4 — "the corresponding legend block must update":** original step 8 says "the corresponding legend block must update to reflect the new column's values". Underspecified what "values" means — is it (a) the legend categories (driven by `Split = Series`, which doesn't change with Y replacement), (b) the Y axis values shown in legend tooltip, or (c) the legend block's title showing the new Y column name? Migrated assertion checks (c) (block title reflects new Y column); (a) is a non-change with `Split = Series` constant. Flag for downstream Automator: pin the verification mode based on actual line-chart legend implementation.
