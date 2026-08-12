# Stages in Diff Studio (Acid Production) — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open Diff Studio + Acid Production from Library | 20s | PASS | PASSED | `runDiffStudio` + dblclick `.diff-studio-hub-card` (hasText 'Acid Production'). The library card opens a view named **'GA-production'** (not 'Acid Production'), with 19 input hosts (`1-st-stage`, `overall`, `biomass`, `glucose`, `oxygen`, `acid`, `muM`, `alpha`, `beta`, `gamma`, `lambda`, `delta`, `phi`, `Ks`, `Ko`, `Kla`, `Cod`, `initial`, `step`) |
| 2 | Multiaxis and Facet tabs present | 5s | PASS | PASSED | Tabs: 'Multiaxis', 'Facet', 'Grid' — both Multiaxis and Facet render line charts with Aspergillus niger model curves |
| 3 | Modify `1-st stage` input; URL updates live | 15s | PASS | PASSED | `[name="input-host-1-st-stage"] input.ui-input-editor`: typed 50 (was 60), Tab. URL settles to `1-ststage=50`. Model re-simulates in real time. No clicker icons on this input — it uses a numeric textbox + range slider combo (step=0.6, min=20, max=80) |
| 4 | Tooltips on inputs | 10s | PASS | PASSED | Dispatched `mouseover` on labels: `1-st-stage` → "Duration of the 1-st stage"; `biomass` → "Aspergillus niger biomass"; `glucose` → "Glucose". All three render in `.d4-tooltip` within ~900ms |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 9s |
| grok-browser execution (scenario steps) | 40s |
| Execute via grok-browser (total) | 1m 49s |
| Spec file generation | 47s |
| Spec script execution | 24s |
| **Total scenario run (with model)** | 3m |

## Summary

Stages (Acid Production) scenario reproduces fully on dev.datagrok.ai. The library card opens a view named "GA-production" internally — scenario wording "Acid Production" refers to the card label, not the view title. All 19 input hosts are reachable via `[name="input-host-*"]`, `1-st-stage` accepts a direct numeric edit (via keyboard, no clickers on this input), URL reflects the change live, and tooltips on `1-st-stage`/`biomass`/`glucose` render descriptive help text. Playwright spec PASSED cleanly in 20.7s. **Total scenario run (with model)**: 3m.

## Retrospective

### What worked well
- `DG.Func.find({name:'runDiffStudio', package:'DiffStudio'})[0].prepare() + call() + grok.shell.addView()` is the canonical entry for any library scenario — identical pattern worked across all 8 DiffStudio scenarios
- `1-st-stage` naming (with hyphens) maps cleanly to `[name="input-host-1-st-stage"]` — no special escaping needed
- Tooltips on labeled inputs are consistently available via raw `mouseover` + ~900ms wait; descriptive text is specific enough to assert per-input

### What did not work
- View name mismatch: scenario says "Acid Production" but `grok.shell.v.name` is "GA-production" (the underlying model function name). Caused an unnecessary double-check in step 1
- Numeric inputs in this model do not expose `+`/`-` clickers — only a slider track. "Use clickers or simply change numbers" in step 3 misleads testers who look for `[name="icon-plus"]`

### Suggestions for the platform
- Make the view name match the library card label (display as "Acid Production" instead of "GA-production" internally) so scenarios and automation can assert on a single string
- Apply `[name="icon-plus"]` / `[name="icon-minus"]` clickers uniformly to integer AND float inputs — currently some have them (PK-PD `count`) and some don't (Acid Production `1-st-stage`)

### Suggestions for the scenario
- Step 1 should mention the view's internal name discrepancy — tester assertions like "view title = Acid Production" will fail; the correct title is "GA-production"
- Step 3 wording "use clickers or simply change numbers" implies clickers exist for every input — list which inputs have clickers (none in Acid Production; just slider + text) vs. keyboard-only edits
- Step 4 lists "various input fields" — pick 3–4 specific inputs with their expected tooltip text (e.g., `1-st-stage` → "Duration of the 1-st stage") so the check is concrete
