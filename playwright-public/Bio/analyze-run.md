# Bio Analyze — Run Results

**Date**: 2026-04-23
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | FASTA.csv: Bio > Analyze > Sequence Space (defaults) | 9s | PASS | PASSED | 9 cols after (Embed_X/Y, Cluster added); UMAP + Hamming |
| 2 | FASTA.csv: Bio > Analyze > Activity Cliffs (defaults) | 4s | PASS | PASSED | 12 cols after; sali and Embed_X_2/Y_2 added |
| 3 | FASTA.csv: Bio > Analyze > Composition | 2s | PASS | PASSED | WebLogo viewer opened |
| 4 | FASTA.csv: Check Composition viewer props (Gear → Context Panel) | 2s | PASS | PASSED | Gear icon absent on WebLogo title bar; used right-click → Properties... fallback. Panel shows Data / Layout / Behavior / Style sections |
| 5 | FASTA.csv: Sequence Space with t-SNE + Levenshtein | 3s | PASS | PASSED | Scatter plot rebuilt; cols=15 (second embed set) |
| 6 | HELM.csv: Sequence Space (defaults) | 12s | PASS | PASSED | 5 cols after (was 2) |
| 7 | HELM.csv: Activity Cliffs (defaults) | 6s | PASS | PASSED | 7 cols after |
| 8 | HELM.csv: Composition | 2s | PASS | PASSED | WebLogo viewer opened |
| 9 | MSA.csv: Sequence Space (defaults) | 12s | PASS | PASSED | 5 cols after |
| 10 | MSA.csv: Activity Cliffs (defaults) | 6s | PASS | PASSED | 7 cols after |
| 11 | MSA.csv: Composition | 2s | PASS | PASSED | WebLogo viewer opened |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~1m 30s |
| grok-browser execution (scenario steps) | ~1m 10s |
| Execute via grok-browser (total) | ~2m 40s |
| Spec file generation | ~20s |
| Spec script execution | 1m 21s |
| **Total scenario run (with model)** | ~4m 20s |

## Summary

All three Bio > Analyze functions (Sequence Space, Activity Cliffs, Composition) work correctly on all three dataset types (FASTA, HELM, MSA). Dialogs open with correct defaults; default and arbitrary-parameter runs both complete within a few seconds, and the resulting viewers (Scatter plot, WebLogo) are added to the view. Regenerated Playwright spec passed all 11 steps in 1m 21s. Total scenario run (with model): ~4m 20s.

## Retrospective

### What worked well
- Bio > Analyze menu items (`div-Bio---Analyze---Sequence-Space...`, `...Activity-Cliffs...`, `...Composition`) continue to work reliably via `click()` + `mouseenter` dispatch pattern.
- Dataset-loading JS path (`grok.dapi.files.readCsv` + `onSemanticTypeDetected` + 5s Bio/Chem settle) opened all three sample files consistently.
- The single-test-per-dataset Playwright structure (parameterised over FASTA/HELM/MSA) kept test isolation clean and produced 3 clear pass lines.
- Re-using the shared `spec-login.ts` helper made the spec preamble trivial and dodged the login-hang footguns called out in the skill.

### What did not work
- WebLogo viewer does **not** render title-bar icons (gear, close, help, maximize) even with `body.selenium` class — same as 2026-04-07 run. Had to reach the properties via right-click → Properties... context menu instead of clicking the gear.

### Suggestions for the platform
- Make WebLogo viewer honour the `selenium` class and render the standard title-bar icons, so automation (and real users discovering the viewer) have a consistent affordance to open settings.
- Consider adding a `gear` icon region to all viewer types by default, or explicitly documenting the WebLogo exception.

### Suggestions for the scenario
- Scenario says "sample_FASTA.csv / sample_HELM.csv / sample_MSA.csv" but actual files in `System:AppData/Bio/samples/` are `FASTA.csv` / `HELM.csv` / `MSA.csv` (no `sample_` prefix). Update the wording, or list the full System paths.
- Step "Run the function once more with arbitrary changed parameters" is ambiguous. Pin concrete alternatives per function (e.g. Sequence Space: `t-SNE + Levenshtein`; Activity Cliffs: `Similarity cutoff = 50`; Composition: `Vertical Alignment = Bottom`) so the test is reproducible.
- Clarify that the Gear check in step 4 applies only to Composition (not to Sequence Space / Activity Cliffs), and that on the current build the entry point for WebLogo properties is the right-click menu.
- Add explicit pre-conditions: file list + expected semantic type (`Macromolecule`) for the main column.
