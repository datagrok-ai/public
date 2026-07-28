# Eval report — chem-toolkit

Started: 2026-05-14T13:41:42.105Z
Duration: 117.6s
Model: `sonnet`
System prompt size: 47,817 chars
Inlined skills: `datagrok-exec`, `chem-toolkit`

## Summary

- **Total prompts:** 19
- **Passed:** 16
- **Failed:** 3
- **Pass rate:** 84.2%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 7 | 7 | 100% |
| cross-package | 4 | 5 | 80% |
| footgun | 2 | 2 | 100% |
| identifier-resolution | 1 | 2 | 50% |
| wrong-tool-trap | 2 | 3 | 67% |

### Latency

- Total: 117.6s wall
- Mean: 17.5s · Median: 12.5s · p95: 73.2s
- Slowest 3: `wrong-tool-dock-as-chem` (73.2s), `footgun-chemprop-train` (60.1s), `chembl-id-then-properties` (18.6s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`chem-toolkit`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| lipinski-properties | atomic-do | pass | pass | pass | pass | **PASS** | 8.8s |
| smiles-to-inchi | atomic-do | pass | pass | pass | pass | **PASS** | 9.9s |
| convert-notation-molblock | atomic-do | pass | pass | pass | pass | **PASS** | 13.3s |
| chemical-space-umap | atomic-do | pass | pass | pass | pass | **PASS** | 15.9s |
| similarity-search | atomic-do | pass | pass | pass | pass | **PASS** | 9.2s |
| substructure-search-benzene | atomic-do | pass | pass | pass | pass | **PASS** | 11.1s |
| rgroup-decomposition | atomic-do | pass | pass | pass | pass | **PASS** | 11.7s |
| admet-aspirin | cross-package | pass | pass | pass | pass | **PASS** | 9.7s |
| admet-table-full-panel | cross-package | pass | pass | pass | pass | **PASS** | 6.6s |
| reinvent-analogs | cross-package | pass | fail | fail | pass | **FAIL** | 15.0s |
| docking-autodock | cross-package | pass | pass | pass | pass | **PASS** | 13.5s |
| boltz1-docking | cross-package | pass | pass | pass | pass | **PASS** | 12.5s |
| chembl-id-then-properties | identifier-resolution | pass | pass | pass | pass | **PASS** | 18.6s |
| drugname-to-smiles-then-mw | identifier-resolution | pass | fail | pass | pass | **FAIL** | 15.2s |
| wrong-tool-admet-as-chem | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 8.0s |
| wrong-tool-generate-as-chem | wrong-tool-trap | pass | fail | fail | pass | **FAIL** | 6.5s |
| wrong-tool-dock-as-chem | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 73.2s |
| footgun-logs-property | footgun | pass | pass | pass | pass | **PASS** | 13.4s |
| footgun-chemprop-train | footgun | pass | pass | pass | pass | **PASS** | 60.1s |

## Failed prompts (deep)

### reinvent-analogs — cross-package

**Prompt:** Generate 10 variations of the seed SMILES CCN(CC)CC for the BACE1 target.
**Intent:** Analog generation → Reinvent4. Must NOT confuse with MMP analysis or reaction enumeration (which derive analogs through different mechanisms).
**Trap:** 'Generate variations / analogs / optimize seed' = Reinvent4. Reactions in Chem are for explicit SMARTS transformations, not generative.

**Rubric:**

- path: PASS
- helpers: FAIL (`Reinvent4:reinvent|Reinvent4:reinventTopMenu|call_function.*Reinvent4:(reinvent|reinventTopMenu)`)
- symbols: FAIL (`ligand`, `optimize`, `CCN\(CC\)CC`)
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`Reinvent4:reinvent|Reinvent4:reinventTopMenu|call_function.*Reinvent4:(reinvent|reinventTopMenu)`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`ligand`, `optimize`, `CCN\(CC\)CC`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
const folders = await grok.functions.call('Reinvent4:getFolders');
return ui.divV(folders.map(f => ui.divText(f)));
```

<details><summary>Full response text</summary>

```
I'll first get the available optimization targets to find the BACE1 folder, then generate the molecules.

```datagrok-exec
const folders = await grok.functions.call('Reinvent4:getFolders');
return ui.divV(folders.map(f => ui.divText(f)));
```
```

</details>

### drugname-to-smiles-then-mw — identifier-resolution

**Prompt:** Compute MW and LogP for ibuprofen.
**Intent:** Drug name → SMILES via Chem:namesToSmiles (ChEMBL-backed), then properties. namesToSmiles is preferred over freeTextToSmiles (LLM fallback) for clean named compounds.

**Rubric:**

- path: PASS
- helpers: FAIL (`Chem:namesToSmiles|call_function.*Chem:namesToSmiles`)
- symbols: PASS
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`Chem:namesToSmiles|call_function.*Chem:namesToSmiles`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper.

**Emitted exec blocks:**

```javascript
// Ibuprofen SMILES (well-known structure)
const ibuprofenSmiles = 'CC(C)Cc1ccc(cc1)C(C)C(=O)O';

// Create a small DataFrame with one molecule
const df = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('molecule', [ibuprofenSmiles])
]);
df.columns.byName('molecule').semType = DG.SEMTYPE.MOLECULE;

// Get MW and LogP
const mol = df.col('molecule');
const propsDf = await grok.functions.call('Chem:getProperties', {
  molecules: mol,
  selected: ['MW', 'LogP']
});

// Build result map
const mw = propsDf.col('MW').get(0);
const logp = propsDf.col('LogP').get(0);

return ui.tableFromMap({
  'Compound': 'Ibuprofen',
  'SMILES': ibuprofenSmiles,
  'MW (Da)': mw?.toFixed(4),
  'LogP': logp?.toFixed(4)
});
```

<details><summary>Full response text</summary>

```
I'll use the known Chem package functions to:
1. Convert ibuprofen's name to SMILES
2. Compute MW and LogP from that SMILES

```datagrok-exec
// Ibuprofen SMILES (well-known structure)
const ibuprofenSmiles = 'CC(C)Cc1ccc(cc1)C(C)C(=O)O';

// Create a small DataFrame with one molecule
const df = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('molecule', [ibuprofenSmiles])
]);
df.columns.byName('molecule').semType = DG.SEMTYPE.MOLECULE;

// Get MW and LogP
const mol = df.col('molecule');
const propsDf = await grok.functions.call('Chem:getProperties', {
  molecules: mol,
  selected: ['MW', 'LogP']
});

// Build result map
const mw = propsDf.col('MW').get(0);
const logp = propsDf.col('LogP').get(0);

return ui.tableFromMap({
  'Compound': 'Ibuprofen',
  'SMILES': ibuprofenSmiles,
  'MW (Da)': mw?.toFixed(4),
  'LogP': logp?.toFixed(4)
});
```

This computes **MW** and **LogP** for ibuprofen using RDKit via `Chem:getProperties`. The expected values are approximately:

| Property | Value |
|---|---|
| **MW** | ~206.28 Da |
| **LogP** | ~3.5 |
```

</details>

### wrong-tool-generate-as-chem — wrong-tool-trap

**Prompt:** Generate 20 new analogs of this lead molecule.
**Intent:** Trap: 'generate analogs' is generative chemistry (Reinvent4), not reaction enumeration or MMP / R-group analysis (which describe existing series, not invent new ones).
**Trap:** Generation = Reinvent4. Reaction enumeration / MMP / R-group operate on existing molecules.

**Rubric:**

- path: PASS
- helpers: FAIL (`Reinvent4:reinvent|Reinvent4:reinventTopMenu|call_function.*Reinvent4:(reinvent|reinventTopMenu)`)
- symbols: FAIL (`ligand`)
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`Reinvent4:reinvent|Reinvent4:reinventTopMenu|call_function.*Reinvent4:(reinvent|reinventTopMenu)`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`ligand`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
const mol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);
if (!mol) throw new Error('No molecule column found in the current table.');
const seed = mol.get(t.currentRowIdx >= 0 ? t.currentRowIdx : 0);
return ui.divText(`Seed molecule (row ${t.currentRowIdx >= 0 ? t.currentRowIdx : 0}): ${seed}`);
```

<details><summary>Full response text</summary>

```
I need to find a molecule to use as the seed. Let me check what's in the current table.

```datagrok-exec
const mol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);
if (!mol) throw new Error('No molecule column found in the current table.');
const seed = mol.get(t.currentRowIdx >= 0 ? t.currentRowIdx : 0);
return ui.divText(`Seed molecule (row ${t.currentRowIdx >= 0 ? t.currentRowIdx : 0}): ${seed}`);
```
```

</details>

## Suggested skill improvements

- Expected helper `Reinvent4:reinvent|Reinvent4:reinventTopMenu|call_function.*Reinvent4:(reinvent|reinventTopMenu)` not reached for: reinvent-analogs, wrong-tool-generate-as-chem. Consider tightening the skill's routing table or example density for this case.
- Expected helper `Chem:namesToSmiles|call_function.*Chem:namesToSmiles` not reached for: drugname-to-smiles-then-mw. Consider tightening the skill's routing table or example density for this case.

## Anomalies

- None.
