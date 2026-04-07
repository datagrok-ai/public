# Import SWAGGER — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: FAIL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Find openweathermap.yaml in Samples package | SKIP | SKIP | File exists at `packages/Samples/swaggers/openweathermap.yaml` in the repository |
| 2 | Download raw file to local machine | SKIP | SKIP | Cannot be automated — requires manual browser download to local filesystem |
| 3 | Drag-and-drop yaml file to Datagrok | SKIP | SKIP | Cannot be automated in this headless/remote browser environment — requires local file access |
| 4 | Go to Browse > Platform > Functions > OpenAPI > OpenWeatherMap | SKIP | SKIP | Skipped because step 3 could not be performed |
| 5 | Right-click connection and select Edit | SKIP | SKIP | Skipped |
| 6 | Enter the ApiKey | SKIP | SKIP | Skipped — ApiKey requires QA team |
| 7 | Run all queries | SKIP | SKIP | Skipped |

## Summary

All 7 steps skipped. This scenario requires manual interaction: downloading a YAML file to the local machine and drag-dropping it into the browser. This is not automatable in a remote browser automation context. The scenario also requires a private OpenWeatherMap API key from the QA team.

## Retrospective

### What worked well
- N/A (scenario not executed)

### What did not work
- Drag-and-drop file upload from local filesystem not automatable via JS evaluate_script in a remote browser session

### Suggestions for the platform
- Add a UI option to import a SWAGGER/OpenAPI file via URL (without requiring local download + drag-drop)
- Alternatively, provide a "File picker" button in the SWAGGER import workflow

### Suggestions for the scenario
- Add a note: "This step requires a local file — cannot be automated without file system access"
- Consider providing the openweathermap.yaml as a file already in DemoFiles so it can be imported from within Datagrok
- The ApiKey should be stored in a shared test credential store rather than requiring QA manually
