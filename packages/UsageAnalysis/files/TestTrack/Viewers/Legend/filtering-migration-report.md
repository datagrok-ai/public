# Migration Report — filtering

## Step mapping

(Original numbering: Section 1 steps 1–14, Section 2 steps 1–3.)

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| 1.1 "Open SPGI" | Setup step 1 | preserved |
| 1.2 "Add seven viewers..." | Setup step 2 | preserved |
| 1.3 "Set the categorical legend on each viewer to Stereo Category..." | Setup step 3 | preserved |
| 1.4 "Open the Filter Panel" | Setup step 4 | preserved |
| 1.5 "Apply a numerical filter Average Mass > 400 — expected ≈ 1588 / 3624 rows; every viewer's legend must show only categories present in the filtered subset" | Scenarios > 1 steps 1–3 | preserved (split for clarity) |
| 1.6 "Apply a categorical filter on Stereo Category — keep only R_ONE and S_UNKN; legend on every viewer should now list exactly two categories" | Scenarios > 1 steps 4–5 | preserved |
| 1.7 "Apply a structure filter on Core — verify viewer + legend update" | Scenarios > 1 steps 6–7 | preserved |
| 1.8 "Save layout → re-apply layout — wait at least 3 s after loadLayout for the filter panel to rebuild before checking that filter state and legends survived" | Scenarios > 2 steps 1–3 | preserved (split for clarity) |
| 1.9 "Reset all filters (df.filter.setAll(true))" | Scenarios > 3 steps 1–2 | preserved |
| 1.10 "On the Scatter plot, set the in-viewer Filter property to ${Stereo Category} in [...] — only the two categories should remain in its legend" | Scenarios > 4 steps 1–2 | preserved |
| 1.11 "While the in-viewer filter is active, add an extra Filter Panel filter Average Mass > 300 — both filters compose; legend stays at two categories" | Scenarios > 5 steps 1–3 | preserved (split for clarity) |
| 1.12 "Filter the table by interacting with each viewer (Scatter alt-drag, Bar OnClick=Filter, Pie OnClick=Filter, Trellis OnClick=Filter)" | Scenarios > 6 steps 1–8 | preserved (4 sub-bullets split for clarity) |
| 1.13 "Save layout → re-apply layout (allow ≥3 s settle) — the filter state from step 12 must be preserved" | Scenarios > 7 steps 1–3 | preserved |
| 1.14 "On the Scatter plot, set Row Source to each of All, Filtered, FilteredSelected, Selected in turn — for each value, note the legend categories shown" | Scenarios > 8 steps 1–8 | preserved (4 modes split for clarity) |
| 2.1 "Continuing from above, on the Bar chart set: Value=CAST Idea ID, Category=Stereo Category, Stack=Primary Scaffold Name" | Scenarios > 9 step 1 | preserved |
| 2.2 "Uncheck Value > Include nulls (bc.props.includeNulls = false)" | Scenarios > 9 step 2 | preserved |
| 2.3 "In the Filter Panel, deselect a few Primary Scaffold Name categories — the bar chart legend must list only the stack categories that are still drawn (no ghost entries from the deselected ones)" | Scenarios > 9 steps 3–4 | preserved (split for clarity) |

## Decisions

- **Why this `target_layer`:** chose `playwright` per chain `output_plan` — Filter Panel UI driving + in-viewer filter property setter + click-to-filter on multiple viewer types + 2 layout save/reapply cycles + Row Source enum cycling. Persistence + UI driving combination requires Playwright.
- **Why this `strategy`:** `chained_tests` per chain — Scenario 9 explicitly continues from Scenario 8's Bar chart state. Single file with shared `beforeAll` for the 7-viewer setup; Scenario 9 references the live Bar chart from the test environment.
- **Why this `pyramid_layer`:** `integration` per chain Rule 4 — multi-subsystem co-existence (3 filter source mechanisms × 7 viewer types + Row Source mode enum). Not bug-focused; not matrix in the parameter-enumeration sense; not smoke (delegated to visibility-and-positioning.md).
- **Why this `priority`:** `p1` — regression coverage for the legend's reactivity to data-change events. Filter panel + in-viewer + click-to-filter span the three primary user paths to legend updates.
- **Why this `coverage_type`:** `regression` at file level — protects against regression of legend refresh on filter-state changes (the GROK-17222 invariant). **Scenario 9 carries a per-scenario `[coverage_type: edge]` marker post-SR application 2026-05-07** (Critic A returned SCOPE_REDUCTION at gate A round 1: Scenario 9's Bar chart Stack ghost-entry behavior under `includeNulls=false` is a true edge case tied to the GROK-17222 family; per-scenario heading marker satisfies A-STRUCT-03 without splitting the file). SR applied without re-invoking Critic A per orchestrator routing rule (SR does not count as review round).
- **Sibling tests consulted:**
  - `visibility-and-positioning-spec.ts` (section UI smoke) — helper imports + spec shape.
  - `color-consistency-spec.ts` — multi-viewer cross-coordination pattern.
- **Helpers reused:** `loginToDatagrok`, `softStep`, `specTestOptions`, `stepErrors` from `spec-login` (registered). No new helpers proposed.
- **Bug library consulted:** yes — 1 bug intersects sub_features_covered (GROK-17222 affects `legend.refresh.on-data-change`, `legend.column`). Listed in `related_bugs`. Cross-cutting bug-focused spec `legend-grok-17222-spec.ts` proposed by chain — this scenario verifies the positive baseline.
- **Decision log queried:** yes — empty for `legend` feature.
- **UI delegation:** standard legend UI flows (visibility / position / color picker / save-dialog widgets) delegated to `visibility-and-positioning.md` per chain `ui_coverage_plan.delegated_scenarios[]`.
- **Cross-cutting bug citations from chain YAML:** chain `bug_focused_candidates[]` lists GROK-17222 with span `filtering.md:Step 5` (numerical filter) — the migrated body's Scenario 1 step 1 maps to this span. Awareness only; spec-level downstream.

## Opt-outs (SCOPE_REDUCTION proposals)

(none)

## Deferred items (NOT opt-outs)

(none)

## Edge cases

- **Filter panel rebuild latency after `loadLayout`:** preserved as Scenario 2 step 2 ("wait at least 3 s after loadLayout for the filter panel to rebuild") and Scenario 7 step 2 ("allow ≥3 s settle"). Without this wait, layout-restore assertions are flake-prone.
- **In-viewer filter compose with Filter Panel filter:** preserved as Scenario 5 — both filters apply simultaneously without Scatter plot legend regressing to non-filtered categories.
- **Bar chart Stack ghost entries with `includeNulls=false`:** preserved as Scenario 9 — Stack-driven legend must list only the categories actually drawn after Filter Panel deselection. Bug-class edge (legend ghost entries in stacks).
- **`FilteredSelected` Row Source semantics:** flagged in Unresolved ambiguities (behavior when no rows are selected is unclear from atlas).

## Unresolved ambiguities

- **Scenario 8 step 5–6 — `FilteredSelected` behavior with empty selection:** original step 14 says "for each value, note the legend categories shown (it should reflect only the rows in that source)". When no rows are selected, `FilteredSelected` may legitimately be empty. Atlas does not pin this; spec assertion likely needs to handle both branches (empty when no selection / scoped when selection exists). Flag for downstream Automator: assertion must guard against the empty-selection state.
- **Scenario 4 step 2 — cross-viewer behavior of in-viewer filter:** original says "only the two categories should remain in its legend" (singular — only the Scatter plot's legend). Migrated step explicitly notes "other viewers may differ" — confirms that in-viewer filter does NOT propagate to other viewers' legends. Verify with platform reference if assumption is wrong.
- **Scenario 9 step 3 — "deselect a few" ambiguity:** original says "deselect a few Primary Scaffold Name categories" — non-deterministic count. Migrated proposes Automator pick a deterministic subset (e.g. 3 categories) and document in spec header. Final subset selection is downstream.
