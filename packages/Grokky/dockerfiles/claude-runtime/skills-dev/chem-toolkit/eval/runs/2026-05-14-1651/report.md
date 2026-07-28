# Eval report — chem-toolkit

Started: 2026-05-14T13:51:15.636Z
Duration: 59.7s
Model: `sonnet`
System prompt size: 47,467 chars
Inlined skills: `datagrok-exec`, `chem-toolkit`

## Summary

- **Total prompts:** 14
- **Passed:** 13
- **Failed:** 1
- **Pass rate:** 92.9%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 7 | 7 | 100% |
| cross-package | 2 | 2 | 100% |
| footgun | 1 | 2 | 50% |
| identifier-resolution | 2 | 2 | 100% |
| wrong-tool-trap | 1 | 1 | 100% |

### Latency

- Total: 59.7s wall
- Mean: 13.0s · Median: 11.0s · p95: 28.7s
- Slowest 3: `footgun-chemprop-train` (28.7s), `drugname-to-smiles-then-mw` (20.5s), `chembl-id-then-properties` (17.2s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`chem-toolkit`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| lipinski-properties | atomic-do | pass | pass | pass | pass | **PASS** | 9.3s |
| smiles-to-inchi | atomic-do | pass | pass | pass | pass | **PASS** | 10.6s |
| convert-notation-molblock | atomic-do | pass | pass | pass | pass | **PASS** | 9.3s |
| chemical-space-umap | atomic-do | pass | pass | pass | pass | **PASS** | 10.6s |
| similarity-search | atomic-do | pass | pass | pass | pass | **PASS** | 13.3s |
| substructure-search-benzene | atomic-do | pass | pass | pass | pass | **PASS** | 12.0s |
| rgroup-decomposition | atomic-do | pass | pass | pass | pass | **PASS** | 11.0s |
| admet-aspirin | cross-package | pass | pass | pass | pass | **PASS** | 7.2s |
| admet-table-full-panel | cross-package | pass | pass | pass | pass | **PASS** | 10.1s |
| chembl-id-then-properties | identifier-resolution | pass | pass | pass | pass | **PASS** | 17.2s |
| drugname-to-smiles-then-mw | identifier-resolution | pass | pass | pass | pass | **PASS** | 20.5s |
| wrong-tool-admet-as-chem | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 8.3s |
| footgun-logs-property | footgun | pass | pass | pass | pass | **PASS** | 13.7s |
| footgun-chemprop-train | footgun | pass | pass | pass | fail | **FAIL** | 28.7s |

## Failed prompts (deep)

### footgun-chemprop-train — footgun

**Prompt:** Train a Chemprop regression model on this dataset to predict the 'activity' column.
**Intent:** Chemprop training has a 25-arg snake_case signature. predictColumn + dataset_type='regression' + a sensible epochs value should appear. Must NOT confuse with applyChemprop (inference) or isApplicableNN (framework hook).
**Trap:** trainChemprop's signature is unusually long and snake_case; Claude often forgets dataset_type or mistakes apply for train.

**Rubric:**

- path: PASS
- helpers: PASS
- symbols: PASS
- forbidden: FAIL (`Chem:applyChemprop\b`)

**Hypothesis:** Forbidden pattern(s) matched (`Chem:applyChemprop\b`). Skill's anti-pattern callout isn't strong enough — consider making the wrong path throw, or surface a more obvious "DO NOT" callout near the relevant example.

**Emitted exec blocks:**

```javascript
// Find the molecule column and train Chemprop
const molCol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);
if (!molCol) throw new Error('No molecule column found in the table.');

const activityCol = t.col('activity');
if (!activityCol) throw new Error("No 'activity' column found in the table.");

grok.shell.info(`Training Chemprop model: molecules from '${molCol.name}', predicting '${activityCol.name}'...`);

const model = await grok.functions.call('Chem:trainChemprop', {
  df: t,
  predictColumn: activityCol,
  dataset_type: 'regression',
  epochs: 50,
});

grok.shell.info('Chemprop model trained successfully!');

// Apply the model back to the table
const resultDf = await grok.functions.call('Chem:applyChemprop', {
  df: t,
  model: model,
});

return DG.Viewer.grid(resultDf).root;
```

<details><summary>Full response text</summary>

```
I'll train a Chemprop regression model on the current table, using the molecule column to predict `activity`.

```datagrok-exec
// Find the molecule column and train Chemprop
const molCol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);
if (!molCol) throw new Error('No molecule column found in the table.');

const activityCol = t.col('activity');
if (!activityCol) throw new Error("No 'activity' column found in the table.");

grok.shell.info(`Training Chemprop model: molecules from '${molCol.name}', predicting '${activityCol.name}'...`);

const model = await grok.functions.call('Chem:trainChemprop', {
  df: t,
  predictColumn: activityCol,
  dataset_type: 'regression',
  epochs: 50,
});

grok.shell.info('Chemprop model trained successfully!');

// Apply the model back to the table
const resultDf = await grok.functions.call('Chem:applyChemprop', {
  df: t,
  model: model,
});

return DG.Viewer.grid(resultDf).root;
```

Here's what this does:

1. **Finds the molecule column** (by semantic type `Molecule`) and the `activity` column in the current table.
2. **Trains a Chemprop regression model** (`dataset_type: 'regression'`) for 50 epochs using the Chemprop message-passing neural network (MPNN).
3. **Applies the trained model** back to the dataset and displays the predictions.

> **Note:** Training runs in the `chem-chemprop` Docker container — it may take a minute or two depending on dataset size. You'll see a notification when it completes. If you'd like to adjust hyperparameters (e.g., more `epochs`, custom `dropout`, `depth`, etc.), just let me know!
```

</details>

## Suggested skill improvements

- Forbidden pattern `Chem:applyChemprop\b` slipped through for: footgun-chemprop-train. Anti-pattern callout needs to be more prominent — consider a "WRONG / RIGHT" pair right before the relevant example.

## Anomalies

- None.
