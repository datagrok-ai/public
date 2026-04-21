# Pareto Front Viewer — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open cars-with-missing.csv from Demo Files | 2s | FAIL | FAILED | `grok.dapi.files.readCsv('System:DemoFiles/cars-with-missing.csv')` returns 502. Listing `System:DemoFiles/` shows only `cars.csv`, `carb_acids.csv`, `car[1-5].jpg` — the named file is absent on dev. |
| 2 | Add Pareto Front viewer, open Properties panel | 0s | SKIP | SKIPPED | Blocked by step 1 (no dataset to attach the viewer to). |
| 3 | Minimize/Maximize lists exclude empty & string cols | 0s | SKIP | SKIPPED | Blocked by step 1 — needs `cars-with-missing.csv` specifically for the empty-column case. |
| 4 | Select all in Maximize; expect warning | 0s | SKIP | SKIPPED | Blocked by step 1. |
| 5 | Open cars.csv; add Pareto Front; "model" auto-selected as Label | 5s | PASS | PASSED | `labelColumnsColumnNames = ["model"]`, `minimizeColumnNames = ["highway.mpg","price"]`, `xAxisColumnName = "highway.mpg"`, `yAxisColumnName = "price"`, `autoLabelsSelection = true`. |
| 6 | Open demog; Label empty by default OR unique-cat auto-selects | 4s | PASS | PASSED | `labelColumnsColumnNames = ["USUBJID"]`; USUBJID has 5850 unique values for 5850 rows — matches the scenario's "auto-select only if unique" branch. |
| 7 | Review all viewer properties (Description, Objectives, Axes, Labels, Legend) | 1s | PASS | PASSED | All five expected property categories present on `paretoV.props.getProperties()`; no exceptions. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 2m |
| grok-browser execution (scenario steps) | 15s |
| Execute via grok-browser (total) | 4m |
| Spec file generation | 2m |
| Spec script execution | 32s |
| **Total scenario run (with model)** | 6m 32s |

All rows are full-phase wall-clock (incl. model thinking and retries), not just tool latency. The two `scenario steps` rows sum to roughly `Execute via grok-browser (total)`.

## Summary

3 of 7 steps passed, 1 failed, 3 were skipped due to the missing prerequisite dataset. The Pareto Front viewer itself is healthy on dev: on `cars.csv` it auto-selects `model` as Label and picks `highway.mpg` / `price` for axes; on `demog.csv` it auto-selects `USUBJID` (correctly — that column has unique values), and its property model exposes the full Description / Objectives / Axes / Labels / Legend categories without errors. The blocker is environmental: `cars-with-missing.csv` is not present under `System:DemoFiles/` on dev, so the empty-column / warning branch of the scenario (steps 1-4) cannot run. **Total scenario run (with model)**: 6m 32s.

## Retrospective

### What worked well
- JS-API driven reproduction (`addViewer('Pareto Front')` + `props.getProperties()`) was fast and gave precise diagnostics.
- The 2b log transcribed cleanly into `softStep()` blocks — spec run outcomes matched the MCP run 1:1.
- `process.env.DATAGROK_URL` + explicit login worked in a fresh Chromium context without any MCP coupling.

### What did not work
- `cars-with-missing.csv` is not deployed under `System:DemoFiles/` on dev, so `readCsv` returned a bare `502`. No client-side cue distinguished "file missing" from "server down".
- Default Playwright `testMatch` does not match `-spec.ts` naming; needed a throwaway `--config` (testMatch `/-spec\.ts$/`) placed at the repo root, plus `NODE_PATH` pointing at the reddata `node_modules` so the config could resolve `@playwright/test`.

### Suggestions for the platform
- `grok.dapi.files.readCsv` on a missing path should return a 404 with the missing path in the message, not a bare `502`. Today the failure mode is indistinguishable from a real server outage.
- Surface a typed error (e.g. `FileNotFoundError` with `path`) from `grok.dapi.files.readCsv` / `list` so scripts can branch without string-matching.

### Suggestions for the scenario
- `cars-with-missing.csv` is referenced but not present under `System:DemoFiles` on dev; add the dataset to DemoFiles or replace the reference with a `cars.csv` variant (e.g. derived in-place by nulling a few cells) that actually exists on the target envs.
- Step 6 wording: the "empty by default" branch is misleading when a unique-valued category column is present (e.g. `USUBJID` on demog). Clarify that auto-select always wins over empty when a unique-value category column exists.
- Step 1: add a precondition line listing the exact file path and fail fast with a clear message if the file is absent.
