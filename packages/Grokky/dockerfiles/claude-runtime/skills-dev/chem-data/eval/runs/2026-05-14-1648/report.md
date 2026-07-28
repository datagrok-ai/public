# Eval report — chem-data

Started: 2026-05-14T13:48:32.320Z
Duration: 78.5s
Model: `sonnet`
System prompt size: 28,817 chars
Inlined skills: `datagrok-exec`, `chem-data`

## Summary

- **Total prompts:** 17
- **Passed:** 17
- **Failed:** 0
- **Pass rate:** 100.0%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 6 | 6 | 100% |
| cross-package | 3 | 3 | 100% |
| footgun | 2 | 2 | 100% |
| identifier-resolution | 4 | 4 | 100% |
| wrong-tool-trap | 2 | 2 | 100% |

### Latency

- Total: 78.5s wall
- Mean: 13.2s · Median: 11.5s · p95: 48.1s
- Slowest 3: `wrong-tool-thrombin-compounds` (48.1s), `curves-ic50-from-data` (17.5s), `inchikey-to-chembl` (16.4s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`chem-data`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| chembl-substructure-search | atomic-do | pass | pass | pass | pass | **PASS** | 7.4s |
| chembl-similarity-search | atomic-do | pass | pass | pass | pass | **PASS** | 7.9s |
| chembl-compounds-by-organism | atomic-do | pass | pass | pass | pass | **PASS** | 10.6s |
| moltrack-register-compound | atomic-do | pass | pass | pass | pass | **PASS** | 7.8s |
| moltrack-lookup-corporate-id | atomic-do | pass | pass | pass | pass | **PASS** | 6.6s |
| curves-fit | atomic-do | pass | pass | pass | pass | **PASS** | 16.3s |
| chemblid-to-smiles | identifier-resolution | pass | pass | pass | pass | **PASS** | 7.8s |
| drugname-to-smiles | identifier-resolution | pass | pass | pass | pass | **PASS** | 8.4s |
| inchikey-to-chembl | identifier-resolution | pass | pass | pass | pass | **PASS** | 16.4s |
| batch-chembl-ids-to-smiles | identifier-resolution | pass | pass | pass | pass | **PASS** | 14.7s |
| curves-ic50-from-data | cross-package | pass | pass | pass | pass | **PASS** | 17.5s |
| hittriage-start-campaign | cross-package | pass | pass | pass | pass | **PASS** | 6.2s |
| moltrack-advanced-search | cross-package | pass | pass | pass | pass | **PASS** | 12.4s |
| wrong-tool-thrombin-compounds | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 48.1s |
| wrong-tool-lookup-by-name | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 13.5s |
| footgun-chembl-namespace | footgun | pass | pass | pass | pass | **PASS** | 11.8s |
| footgun-names-to-smiles-depends-on-chembl | footgun | pass | pass | pass | pass | **PASS** | 11.5s |

## Suggested skill improvements

- No structural issues detected. Either the eval passed broadly, or the failures are one-offs requiring case-by-case review.

## Anomalies

- `wrong-tool-thrombin-compounds`: emitted 6 exec blocks — unusually verbose; check whether the skill should encourage consolidation.
