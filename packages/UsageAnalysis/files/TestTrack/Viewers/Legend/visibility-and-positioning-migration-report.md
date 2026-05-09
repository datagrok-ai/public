# Migration Report — visibility-and-positioning

## Step mapping

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| 1. "Open SPGI" | Setup step 1 | preserved |
| 2. "Add seven viewers: Scatter plot, Histogram, Line chart, Bar chart, Pie chart, Trellis plot, Box plot" | Setup step 2 | preserved |
| 3. "Arrange viewers... 2×4 grid with the dock manager" | Setup step 3 | preserved |
| 4. "Set the categorical legend on each viewer to Stereo Category" | Setup step 4 | preserved |
| 5. "Verify the legend is visible on every viewer (DOM `[name='legend']` count == 7) and the swatch colors match" | Scenarios > 1 steps 1–2 | preserved as verification |
| 6. "On any viewer, change the legend source to a different column... and back to Stereo Category — must redraw with the new categories each time" | Scenarios > 2 steps 1–4 | preserved (split for clarity) |
| 7. "Adjust the legend size — verify it shrinks/grows" | Scenarios > 3 steps 1–4 | preserved (split for clarity) |
| 8. "Ctrl+click a category in the Scatter plot legend... only that category remains; click the X on a swatch to exclude it instead" | Scenarios > 4 steps 1–4 | preserved (split for clarity) |
| 9. "Color picker — change color via dialog — Cancel discards, OK commits, and the new color propagates to every viewer using Stereo Category for legend" | Scenarios > 5 steps 1–7 | preserved (split for clarity) |
| 10. "Switch the legend column to Primary Series Name — the (no value) / empty-value swatch is present and its color can be changed via the same picker" | Scenarios > 6 steps 1–5 | preserved (split for clarity) |
| 11. "Save layout → re-apply layout — legend column, custom colors, and visibility state must all persist" | Scenarios > 7 steps 1–3 | preserved |
| 12. "Right-click the legend → set Visibility = Always, Position = Auto on every viewer — legends must remain visible" | Scenarios > 8 steps 1–2 | preserved |
| 13. "Resize a viewer — with Position = Auto, the legend should reposition to whichever side has the most free space" | Scenarios > 8 steps 3–4 | preserved |
| 14. "Save layout → re-apply layout — Visibility=Always and Position=Auto must persist" | Scenarios > 8 steps 5–6 | preserved |
| 15. "Right-click the legend → uncheck auto-positioning, set Visibility = Auto on every viewer" | Scenarios > 9 step 1 | preserved |
| 16. "Reduce a viewer's size below ~250 px — with Visibility=Auto, the legend must hide; restore to ≥400 px and the legend must reappear" | Scenarios > 9 steps 2–5 | preserved (split for clarity) |
| 17. "Set legendPosition to each of LeftTop, LeftBottom, RightTop, RightBottom, then enable mini-legend mode on a few viewers — verify corner rendering and mini-variant compactness" | Scenarios > 10 steps 1–6 | preserved (split for clarity) |
| 18. "Save layout → re-apply layout — corner position and mini-legend mode must persist" | Scenarios > 10 steps 7–8 | preserved |
| 19. "Save project → reopen project — verify positioning and mini-legend mode survive the persistence round-trip" | Scenarios > 11 steps 1–3 | preserved |
| 20. "Close All" | Scenarios > 11 step 4 | preserved as cleanup |

## Decisions

- **Why this `target_layer`:** chose `playwright` because the scenario requires Right-click context-menu interaction (Visibility / Position properties), color-picker dialog widget, splitter drag, save-layout / reapply / save-project / reopen persistence cycles, and DOM `[name="legend"]` count assertions. None substitutable by JS API alone. Matches chain `output_plan.target_layer: playwright`. Sibling spec `visibility-and-positioning-spec.ts` (existing) already operates at this layer.
- **Why this `strategy`:** `simple` per chain `output_plan.strategy: simple`. Single linear flow; 11 sub-scenarios are organizational beats within the same `test()` block, not independent test entries.
- **Why this `pyramid_layer`:** `ui-smoke` per chain assignment (Rule 1 residual — broadest UI surface in section, no other scenario qualifies). Section's canonical UI smoke; all 5 other scenarios delegate standard legend UI flows to it.
- **Why this `priority`:** `p0` — golden-path smoke covering the section's most-used legend UI flows (visibility / position / color picker / persistence). All 4 atlas `critical_paths` p0 entries (filter-by-legend-click, categorical-color-coded-show, per-category-color-change, plus position-property which is p1) intersect this scenario's flows.
- **Sibling tests consulted:**
  - `public/packages/UsageAnalysis/files/TestTrack/Viewers/Legend/visibility-and-positioning-spec.ts` (existing — will be rewritten by Automator stage; informs spec-style anchoring).
  - `public/packages/UsageAnalysis/files/TestTrack/Viewers/Legend/color-consistency-spec.ts` (cross-viewer color propagation pattern).
  - `public/packages/UsageAnalysis/files/TestTrack/Projects/uploading-spec.ts` (save-dialog + reload assertion shape).
- **Helpers reused:**
  - `loginToDatagrok` (`spec-login.ts:26`) — registered in helpers-registry.yaml.
  - `softStep` (`spec-login.ts:18`) — registered.
  - `specTestOptions` and `stepErrors` from `spec-login` — registered.
  - No new helpers proposed for this scenario.
- **Bug library consulted:** yes — 6 bugs in `bug-library/legend.yaml`. All 6 have `affects` intersecting this scenario's `sub_features_covered`; all 6 listed in frontmatter `related_bugs`. Detail:
  - GROK-17438 (legend.item.color-picker, legend.visibility, legend.splitter-resize) — primary surface for color-change-hides-legend bug, fully testable in this scenario.
  - GROK-17222 (legend.refresh.on-data-change, legend.column) — column intersection only; refresh-on-filter is filtering.md's domain.
  - github-3132 (legend.item.color-picker, legend.allow-item-coloring) — Scenario 5 covers sequential color changes; bug-focused spec proposed for full invariant.
  - GROK-17278 (legend.item.color-picker, legend.column) — partial coverage in Scenario 7 (layout persistence) and Scenario 11 (project persistence).
  - GROK-19083 (legend.show-main-item-icons, legend.column) — column intersection only; markers-deselect is structure-rendering.md's domain.
  - GROK-19041 (legend.auto-position, legend.position) — relevant but regression-line overlap requires regression line setup not in this scenario; bug skipped at chain level (`reproduction_unparseable`).
- **Decision log queried:** yes — empty for `legend` feature (no `migration_decisions`, `layer_decisions`, `failed_attempts`, or `manual_only` entries; first migrate cycle for this section).
- **Cross-cutting bug citations from chain YAML:** 5 cross-cutting candidates reference this scenario in `spans`: GROK-17438 (Step 9 — color picker), GROK-17222 (Step 6 — column switch), github-3132 (Step 9 — sequential color change), GROK-17278 (Step 9 + Step 11 — color persistence), GROK-19083 (Step 4 — but column-only weak match). See chain `bug_focused_candidates[]`. This smoke surfaces awareness only; Automator will produce dedicated `legend-<bug-id>-spec.ts` files downstream.
- **Chain-vs-frontmatter divergence:** chain rev 1 `sub_features_covered` for this scenario did not include `legend.auto-position`; this migration adds it because step 13 (Position = Auto + viewer resize → legend repositions to side with most free space) explicitly tests the auto-position decision logic. Chain regeneration recommended at section-complete to align — minor divergence, no Critic blocker.
- **`pyramid_layer: ui-smoke` constraint compliance:** every flow in `ui_coverage_responsibility[]` (12 flows: pcmdLegendVisibility, pcmdLegendPosition, legend-color-picker-dialog, legend-column-property-selector, legend-resize-handle, legend-mini-icon-toggle, legend-corner-collapse-chevron, legend-item-click-filter, legend-item-cross-click, legend-no-value-swatch, layout-save-reapply, project-save-reopen) has at least one matching scenario step. No JS API substitution applied where UI driving is the contract.

## Opt-outs (SCOPE_REDUCTION proposals)

(none)

_Rationale: every original step has a UI-drivable equivalent; no flow is technically unreachable from Playwright. The smoke is the section's authoritative UI exercise — no opt-outs permitted under `pyramid_layer: ui-smoke` (per migration prompt section "Authoring patterns" Pattern 1 / "REQUIRED actions" b)._

## Deferred items (NOT opt-outs)

(none)

_All 20 original steps mapped 1-to-1 (or 1-to-many for split-for-clarity rows) into the migrated body. No prerequisite blocks any step._

## Edge cases

- **`(no value)` swatch on a column with empty values:** preserved as Scenario 6 (`Primary Series Name` is the canonical column with nulls per the original step 10). Tests `legend.show-nulls` sub-feature behavior.
- **Viewer size < ~250 px auto-hide threshold:** preserved as Scenario 9 step 2 (`viewer.root.style.width = '200px'`). Tests `legend.auto-show` boundary.
- **Color picker `Cancel` vs `OK` semantics:** preserved as Scenario 5 steps 3–4 vs 5–7. Tests that Cancel discards and the original color is retained — bug invariant from github-3132 reinforces additive color-change persistence (covered by chain bug-focused candidate, not duplicated here).
- **Position = Auto under viewer aspect change:** preserved as Scenario 8 step 4. Tests `legend.auto-position` decision logic (Right vs Top per aspect ratio + measured legend box width).
- **Mini-legend mode on a subset of viewers:** preserved as Scenario 10 step 5 with deterministic-subset note in Unresolved ambiguities.

## Unresolved ambiguities

- **Scenario 4 step 1 — `sp.legend.toggle('R_ONE')` JS API exposure:** original step says "Ctrl+click a category in the Scatter plot legend (or call `sp.legend.toggle('R_ONE')` if exposed)". JS API alternative is unconfirmed; UI Ctrl+click is the canonical path for this smoke (per `pyramid_layer: ui-smoke` requirement). Flag for QA pair review — if `sp.legend.toggle()` is exposed and reliable, future revision can use it as a deterministic alternative; for now, UI driving is the contract.
- **Scenario 10 step 5 — "a few viewers" subset for mini-legend:** original step is non-deterministic ("on a few viewers"). Migrated step proposes a deterministic subset to be selected by the Automator (recommendation: 3 viewers — Scatter plot, Bar chart, Pie chart — covering the three most-used host types). Final subset selection is downstream Automator decision; document the chosen subset in spec header comments per Pattern-2 self-doc convention.
- **Scenario 11 step 3 verification depth:** "verify positioning and mini-legend mode survive the persistence round-trip" — original is high-level. Migrated assertion deferred to Automator: at minimum verify (a) corner position values preserved per viewer, (b) mini-legend mode flag preserved per viewer, (c) Visibility / Position values preserved. If any sub-assertion is platform-flake-prone, Automator may downscope to a subset with hypothesis logged.
