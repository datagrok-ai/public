# Migration Report — structure-rendering

## Step mapping

| Original step | Migrated location | Decision |
|---------------|-------------------|----------|
| 1. "Open SPGI — wait for Molecule semantic-type detection on Core / Structure columns" | Setup step 1 | preserved |
| 2. "Add seven viewers: Scatter plot, Histogram, Line chart, Bar chart, Pie chart, Trellis plot, Box plot" | Setup step 2 | preserved |
| 3. "On each viewer, set the legend column to Core (a Molecule-semType column)... Verify that every legend renders rendered molecule thumbnails (not raw SMILES text or category indices)" | Scenarios > 1 steps 1–2 | preserved (split for clarity) |
| 4. "On the Scatter plot, set Marker = Core and Color = Core — markers must render as molecule glyphs" | Scenarios > 2 steps 1–2 | preserved |
| 5. "Change Color to Series (keep Marker = Core) — markers stay as molecule glyphs; legend now shows Series colors alongside the structure markers" | Scenarios > 2 steps 3–5 | preserved (split for clarity) |
| 6. "Save layout → re-apply layout — every legend must still render structures after the round-trip" | Scenarios > 3 steps 1–3 | preserved |
| 7. "Save project → reopen project — verify legends and structure rendering survive the persistence round-trip" | Scenarios > 4 steps 1–3 | preserved |

## Decisions

- **Why this `target_layer`:** chose `playwright` per chain `output_plan` — save-layout / reapply / save-project / reopen persistence cycles + DOM canvas thumbnail assertions require Playwright.
- **Why this `strategy`:** `simple` per chain — single linear flow, no chained sub-state.
- **Why this `pyramid_layer`:** `integration` per chain Rule 4 — multi-subsystem co-existence (7 viewer types × Molecule semantic-type rendering verified across all of them simultaneously). Not a smoke (delegates standard UI to visibility-and-positioning.md), not source-matrix (no parameter enumeration over data sources), not bug-focused (no `related_bugs` triggering Rule 3).
- **Why this `priority`:** `p1` — regression coverage for molecule rendering. Not p0 because the scenario does not exercise the canonical legend UI flows (those live in the smoke), but Molecule semantic-type rendering is a high-value invariant for chemistry workflows.
- **Why this `coverage_type`:** `regression` — protects against rendering regressions in legend / marker glyphs across viewer types when Molecule semType detection changes.
- **Sibling tests consulted:**
  - `visibility-and-positioning-spec.ts` (section UI smoke) — confirms helper imports + spec shape.
  - `color-consistency-spec.ts` — same multi-viewer pattern (6 viewers); informs viewer-iteration helper.
- **Helpers reused:** `loginToDatagrok`, `softStep`, `specTestOptions`, `stepErrors` from `spec-login` (registered). No new helpers proposed.
- **Bug library consulted:** yes — 1 bug intersects sub_features_covered (GROK-19083 affects `legend.show-main-item-icons`, `legend.column`). Listed in `related_bugs`. The cross-cutting bug-focused spec `legend-grok-19083-spec.ts` is proposed by chain; this scenario verifies the positive baseline only.
- **Decision log queried:** yes — empty for `legend` feature.
- **Post-Gate-F frontmatter expansion (2026-05-07):** Critic F at section-complete returned SCOPE_REDUCTION with `F-REDUNDANCY-01` flagging this scenario as a strict-subset of `scatterplot.md` under the original 2-feature `sub_features_covered` set. Frontmatter expanded with `legend.item.marker-picker` (Scenario 2 step 4 — Color = `Series` with Marker = `Core` exercises the marker dimension of the legend; legend.item.marker-picker is the per-category marker picker which renders the molecule glyphs in this scenario's marker mode) and `legend.refresh.on-data-change` (Scenarios 3+4 layout and project save/reload trigger the legend rebuild on data restoration). These are genuine sub-feature exercises; the original migration omitted them as "setup-only", but they are tested invariants — Marker glyph rendering depends on `legend.item.marker-picker`'s code path; persistence rebuild depends on `legend.refresh.on-data-change`. Expansion breaks the strict-subset relation against `scatterplot.md`. SR application 2026-05-07 per Gate F SR routing dispatch (Critic A re-invocation expected — round 2 budget consumed).
- **UI delegation:** standard legend UI flows (visibility / position / color picker / save-dialog widgets) delegated to `visibility-and-positioning.md` per chain `ui_coverage_plan.delegated_scenarios[]`. This scenario does not duplicate those exercises.
- **Cross-cutting bug citations from chain YAML:** chain `bug_focused_candidates[]` lists GROK-19083 with span `structure-rendering.md:Step 4` (Marker=Core, Color=Core) — the migrated body's Scenario 2 step 1 maps to this span. Awareness only; spec-level downstream.

## Opt-outs (SCOPE_REDUCTION proposals)

(none)

_All 7 original numbered steps preserved as scenario steps. No technical or atlas-level dependency blocks any flow._

## Deferred items (NOT opt-outs)

(none)

_All steps preserved._

## Edge cases

- **Molecule semType detection latency at table load:** preserved as Setup step 1 ("wait for `Molecule` semantic-type detection"). Without this wait, downstream legend column assignment may render text instead of thumbnails (atlas-implicit edge — covered).
- **Mixed Color + Marker semantics on Scatter plot:** preserved as Scenario 2 steps 3–5 (Color=Series + Marker=Core combines categorical color with molecule glyphs). Reinforces bug GROK-19083 invariant context (host-viewer property change must reflect in legend without sync drift).
- **Persistence of rendered structures across layout AND project save:** two separate scenarios (3 and 4) — layout save + project save are independent serialization paths; both must round-trip structure rendering correctly. (Implicit edge from original — explicit in migrated.)

## Unresolved ambiguities

- **Scenario 2 step 5 verification depth:** original step says "legend now shows `Series` colors alongside the structure markers" — unclear whether the legend should display structures + colors as a combined item, or separate structure marker + color swatch entries. Atlas's `legend.extra-column` sub-feature (secondary column rendering) is the likely backing behavior, but this scenario does not formally cover `legend.extra-column`. Flag for downstream Automator: render shape may need atlas curation pass to disambiguate.
