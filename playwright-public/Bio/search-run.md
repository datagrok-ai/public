# Bio Subsequence Search — Run Results

**Date**: 2026-04-24
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open sample_FASTA.csv | 6s | PASS | PASSED | Opened `System:AppData/Bio/samples/FASTA.csv`; 64 rows, 6 columns; `Sequence` column detected as `Macromolecule`. |
| 2 | Open Bio > Search > Subsequence Search | 2s | PASS | PASSED | Menu item name is `div-Bio---Search---Subsequence-Search-...` (trailing hyphen — the label is `"Subsequence Search "` with a trailing space). Action adds a single Bio substructure filter card for the `Sequence` column; no dialog opens. |
| 3 | Set Sequence filter to MHAILRYFIRRLFYHIFYKIYSLISKKHQSLPSDVRQF | 3s | PASS | PASSED | Typed into `[name="viewer-Filters"] input[placeholder="Substructure"]`. Table filtered from 64 rows to 1 matching row. |
| 4 | Click Reset Filter | 3s | PASS | PASSED | Reset (↺) is the `icon-arrow-rotate-left` on the filter group header. No confirmation prompt — the substructure input cleared and row count returned to 64. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 45s |
| grok-browser execution (scenario steps) | 14s |
| Execute via grok-browser (total) | 1m 0s |
| Spec file generation | 25s |
| Spec script execution | 14s |
| **Total scenario run (with model)** | 1m 40s |

## Summary

All four steps passed on the dev server with both the MCP reproduction and the generated Playwright spec. The scenario exercises the Bio substructure filter path end-to-end: the Subsequence Search action adds a dedicated filter card for the Macromolecule column, typing a known sequence narrowed 64 FASTA rows to the single expected match, and the global Reset (↺) icon cleared the input and restored the full row count. Total scenario run (with model): **1m 40s**.

## Retrospective

### What worked well
- Existing `spec-login.ts` helper and `softStep`/`stepErrors` pattern made the spec a thin transcription of the MCP run — no login or env-var boilerplate to drift.
- The Bio substructure filter exposes a standard `<input>` with `placeholder="Substructure"`, so programmatic value-setting + synthetic `input`/`Enter` events reliably triggered filtering without needing the full keyboard-type path.
- `df.filter.trueCount` is a stable verification point for both the filtered (1) and reset (64) states.

### What did not work
- The `name=` attribute for the menu item is `div-Bio---Search---Subsequence-Search-...` — the trailing hyphen comes from the label `"Subsequence Search "` having a trailing space, which then gets hyphen-normalised before the `...` dialog suffix. This was not guessable from the scenario text; required one exploratory call to enumerate `[name^="div-Bio---"]` items.

### Suggestions for the platform
- Trim trailing whitespace from menu item captions before generating the `name=` attribute — `div-Bio---Search---Subsequence-Search...` (without the extra hyphen) would match the rest of the Bio menu family and the naming convention in `grok-browser/references/`.
- The Subsequence Search entry carries a `...` dialog suffix in its `name=` but no dialog is actually opened — it acts directly on the filter panel. Either drop the `...` (it breaks the "dialog suffix" convention from `grok-browser/SKILL.md`) or route the action through an explicit dialog.

### Suggestions for the scenario
- Clarify that step 1 refers to `System:AppData/Bio/samples/FASTA.csv` (the scenario says "sample_FASTA.csv" but no file by that exact name exists in the Bio samples folder).
- Add an expected-outcome line to step 3 (e.g. "1 of 64 rows should remain") so failures are easier to diagnose.
- Step 4 should specify which reset — the filter-group Reset (↺) at the top of the filter panel, not a per-card reset — to avoid ambiguity.
