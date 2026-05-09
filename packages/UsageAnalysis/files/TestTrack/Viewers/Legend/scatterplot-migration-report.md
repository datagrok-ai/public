# Migration Report — scatterplot

## Step mapping

(Original numbering uses Section.Step convention: 1.X = Section 1, 2.X = Section 2, etc.)

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| 1.1 "Open SPGI" | Scenarios > 1 step 1 | preserved |
| 1.2 "Add a scatterplot" | Scenarios > 1 step 2 | preserved |
| 1.3 "Set Color and Marker to Series — legend must be combined" | Scenarios > 1 steps 3–4 | preserved (split for clarity) |
| 1.4 "Check color picker visibility, change some colors" | Scenarios > 1 step 5 | preserved |
| 1.5 "Save and apply the layout" | Scenarios > 1 steps 6–7 | preserved (split for clarity) |
| 1.6 "Save and open the project — changes should persist" | Scenarios > 1 steps 8–9 | preserved (split for clarity) |
| 1.7 "On the plot, add new Color columns: linear formula → linear legend; categorical formula → categorical legend" | Scenarios > 1 step 10 | preserved (two formulas as bullet sub-steps) |
| 1.8 "Set Color to ID, Marker to Core — verify legend is displayed" | Scenarios > 1 step 11 | preserved |
| 1.9 "Close All" | Scenarios > 1 step 12 | preserved as cleanup |
| 2.1 "Open SPGI" | Scenarios > 2 step 1 | preserved |
| 2.2 "Add new columns col1/col2 with formulas" | Scenarios > 2 step 2 | preserved (two formulas as bullet sub-steps) |
| 2.3 "Add a scatterplot" | Scenarios > 2 step 3 | preserved |
| 2.4 "Set X axis to col1, Color to Stereo Category" | Scenarios > 2 steps 4–5 | preserved (split for clarity — verification after setup) |
| 2.5 "Change X axis to col2" | Scenarios > 2 step 6 | preserved |
| 2.6 "Verify legend categories update according to data" | Scenarios > 2 step 7 | preserved as verification |
| 2.7 "Test zooming/filtering — verify legend stays consistent" | Scenarios > 2 steps 8–9 | preserved (split for clarity) |
| 2.8 "Close All" | Scenarios > 2 step 10 | preserved as cleanup |
| 3.1 "Open SPGI" | Scenarios > 3 step 1 | preserved |
| 3.2 "Add a scatterplot" | Scenarios > 3 step 2 | preserved |
| 3.3 "Set Marker to Stereo category" | Scenarios > 3 step 3 | preserved |
| 3.4 "Apply in-viewer filter ${Stereo Category} in [R_ONE, S_UNKN]" | Scenarios > 3 steps 4–5 | preserved (split for clarity — verification after setup) |
| 3.5 "Add another scatterplot with the same filter and marker" | Scenarios > 3 step 6 | preserved |
| 3.6 "Save and apply the layout — verify filtered legends" | Scenarios > 3 steps 7–8 | preserved (split for clarity) |
| 3.7 "Close All" | Scenarios > 3 step 9 | preserved as cleanup |
| 4.1 "Open SPGI" | Scenarios > 4 step 1 | preserved |
| 4.2 "Add a scatterplot" | Scenarios > 4 step 2 | preserved |
| 4.3 "Set X/Y axes to Chemical space X/Y" | Scenarios > 4 step 3 | preserved |
| 4.4 "Set Color to Primary scaffold name, Marker to Stereo category" | Scenarios > 4 step 4 | preserved |
| 4.5 "In Filter Panel, filter Primary scaffold name — verify data and legend update" | Scenarios > 4 steps 5–6 | preserved (split for clarity) |
| 4.6 "Click 'R_ONE' in the scatterplot legend — verify correct additional filtering, no reappearing filtered-out points" | Scenarios > 4 steps 7–8 | preserved (split for clarity) |
| 4.7 "Close All" | Scenarios > 4 step 9 | preserved as cleanup |
| 5.1 "Open SPGI" | Scenarios > 5 step 1 | preserved |
| 5.2 "Add a scatterplot, box plot, and PC plot" | Scenarios > 5 step 2 | preserved |
| 5.3 "Set Color to Chemical Space X" | Scenarios > 5 step 3 | preserved |
| 5.4 "In the grid, enable linear color coding for Chemical Space X — check legends on both viewers" | Scenarios > 5 steps 4–5 | preserved (split for clarity) |
| 5.5 "Change color schema, invert, apply to text — check legends" | Scenarios > 5 steps 6–7 | preserved (split for clarity) |
| 5.6 "Save and apply layout — verify color changes" | Scenarios > 5 steps 8–9 | preserved (split for clarity; see Unresolved ambiguities for assertion depth) |
| 5.7 "Change Chemical Space X coding to categorical, modify colors — verify legends update" | Scenarios > 5 steps 10–11 | preserved (split for clarity) |
| 5.8 "Save and open the project — verify changes persist" | Scenarios > 5 steps 12–13 | preserved (split for clarity) |

## Decisions

- **Why this `target_layer`:** chose `playwright` per chain `output_plan` — every section involves UI driving (color picker dialog, in-viewer filter property, Filter Panel, grid color-coding setup) and 4 of 5 sections include save-layout / save-project persistence cycles.
- **Why this `strategy`:** `chained_tests` per chain — 5 independent test blocks in one spec file (each Scenario starts fresh with Open SPGI). Shared spec-level setup (login, page setup) via beforeAll; no inter-Scenario state.
- **Why this `pyramid_layer`:** `integration` per chain Rule 4 — multi-subsystem (scatterplot × Color/Marker combined × axis change × in-viewer filter × Filter Panel × click-to-filter × grid color coding linear-vs-categorical × persistence).
- **Why this `priority`:** `p1` — regression coverage for scatterplot-specific legend behaviors. P0 is reserved for the section ui-smoke and the foundational color-consistency invariant.
- **Why this `coverage_type`:** `regression` at file level — protects against scatterplot legend regressions across multiple behavior dimensions. **Scenarios 1, 3, 5 carry per-scenario `[coverage_type: edge]` heading markers post-SR application 2026-05-07** (Critic A returned SCOPE_REDUCTION at gate A round 1: persistence-cycle scenarios 1 / 3 / 5 substantively test atlas edge_cases tied to GROK-17438 / github-3132 / GROK-17278 / GROK-17222 invariants; per-scenario heading markers satisfy A-STRUCT-02 without relabeling the file). SR applied without re-invoking Critic A per orchestrator routing rule (SR does not count as review round).
- **Sibling tests consulted:**
  - `visibility-and-positioning-spec.ts` (section UI smoke) — helper imports + spec shape.
  - `scatterplot-spec.ts` (existing) — same scenario; informs spec-style anchoring for downstream Automator rewrite.
  - `color-consistency-spec.ts` — multi-viewer cross-coordination pattern (relevant for Section 5 multi-viewer color coding).
- **Helpers reused:** `loginToDatagrok`, `softStep`, `specTestOptions`, `stepErrors` from `spec-login` (registered). No new helpers proposed.
- **Bug library consulted:** yes — 5 bugs intersect sub_features_covered:
  - GROK-17438 (color-picker / visibility / splitter-resize) — color-picker intersection (Section 1).
  - GROK-17222 (refresh.on-data-change / column) — both intersect (Section 4 covers Filter Panel filtering directly).
  - github-3132 (color-picker / allow-item-coloring) — direct match (Section 1 sequential color changes).
  - GROK-17278 (color-picker / column) — both intersect (Section 1 layout/project persistence).
  - GROK-19083 (show-main-item-icons / column) — both intersect (Section 1 step 11 Marker=Core; Section 3 Marker=Stereo).
- **Decision log queried:** yes — empty for `legend` feature.
- **UI delegation:** standard legend UI flows (visibility / position / save-dialog widgets) delegated to `visibility-and-positioning.md`. Color picker dialog usage and filter panel + click-to-filter flows are owned by this scenario.
- **Cross-cutting bug citations from chain YAML:** chain `bug_focused_candidates[]` lists 4 entries with spans into this scenario (GROK-17438:Step 4, GROK-17222:Step 5, github-3132:Step 4, GROK-17278:Step 5). Migrated body's Scenario 1 step 5 + Scenario 4 step 5 + Scenario 1 step 6 are the corresponding migrated locations. Awareness only; spec-level downstream.
- **Section independence:** original explicitly uses "Close All" between sections, signaling independent test scopes. Preserved in migrated as 5 separate Scenarios with own Open SPGI setup.

## Opt-outs (SCOPE_REDUCTION proposals)

(none)

## Deferred items (NOT opt-outs)

(none)

## Edge cases

- **Combined Color + Marker legend:** preserved as Scenario 1 step 4. The "combined" semantics (one entry per category showing both color and marker glyph) is the atlas `legend.extra-column` sub-feature behavior — implicit edge captured.
- **Derived nullable color columns (linear vs categorical legend dispatch):** preserved as Scenario 1 step 10. The `if(...== 'S_UNKN', null, ...)` formula introduces nulls; the legend type (numerical scale vs categorical items) depends on the formula's output column type — implicit edge in atlas `legend.color-scale.numerical` / `legend.use-custom-color-coding` branch.
- **In-viewer filter on multiple scatterplots in same layout:** preserved as Scenario 3 — two scatterplots with the same filter test layout-level filter persistence and per-viewer legend correctness.
- **Click-to-filter on legend composes with Filter Panel filter:** preserved as Scenario 4 step 8 — explicit assertion that legend click-filter does not "reset" Filter Panel filter; the two filters compose.
- **Grid color coding linear → categorical mode switch:** preserved as Scenario 5 step 10. The mode switch is a coloring-paradigm change; implicit edge that legend type must update accordingly.

## Unresolved ambiguities

- **Scenario 5 step 9 — assertion depth on layout-persisted color changes:** original says "Save and apply layout — verify color changes". Underspecified: is the assertion (a) the linear color schema persists, (b) specific custom colors persist, or (c) both? Migrated proposes assertion includes both schema choice + invert state + text-application setting. Flag for downstream Automator: pin a specific verifiable subset.
- **Scenario 4 step 8 — "no reappearing filtered-out points" semantics:** original step 6 in Section 4 says "verify correct additional filtering, no reappearing filtered-out points". Could mean (a) Filter Panel filter is preserved after legend click, or (b) no race condition where legend click momentarily un-applies the panel filter. Migrated assertion is the former (composition); flag for downstream if (b) is the intended invariant.
- **Scenario 5 step 5 — "both viewers" (PC plot vs box plot):** original step 4 says "check legends on both viewers". With three viewers added (scatterplot, box plot, PC plot), "both" is ambiguous. Migrated assertion checks scatterplot + box plot (the two with categorical legend behavior); PC plot's parallel-coordinates rendering may not have a categorical legend in the same form. Flag for downstream if PC plot is the intended target instead of box plot.
