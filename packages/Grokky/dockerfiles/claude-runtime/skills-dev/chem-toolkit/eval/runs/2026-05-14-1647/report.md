# Eval report — chem-toolkit

Started: 2026-05-14T13:47:15.099Z
Duration: 76.6s
Model: `sonnet`
System prompt size: 47,817 chars
Inlined skills: `datagrok-exec`, `chem-toolkit`

## Summary

- **Total prompts:** 19
- **Passed:** 14
- **Failed:** 5
- **Pass rate:** 73.7%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 7 | 7 | 100% |
| cross-package | 4 | 5 | 80% |
| footgun | 1 | 2 | 50% |
| identifier-resolution | 1 | 2 | 50% |
| wrong-tool-trap | 1 | 3 | 33% |

### Latency

- Total: 76.6s wall
- Mean: 13.7s · Median: 13.1s · p95: 24.0s
- Slowest 3: `footgun-chemprop-train` (24.0s), `reinvent-analogs` (23.6s), `smiles-to-inchi` (22.8s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`chem-toolkit`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| lipinski-properties | atomic-do | pass | pass | pass | pass | **PASS** | 9.3s |
| smiles-to-inchi | atomic-do | pass | pass | pass | pass | **PASS** | 22.8s |
| convert-notation-molblock | atomic-do | pass | pass | pass | pass | **PASS** | 11.6s |
| chemical-space-umap | atomic-do | pass | pass | pass | pass | **PASS** | 13.1s |
| similarity-search | atomic-do | pass | pass | pass | pass | **PASS** | 14.0s |
| substructure-search-benzene | atomic-do | pass | pass | pass | pass | **PASS** | 15.5s |
| rgroup-decomposition | atomic-do | pass | pass | pass | pass | **PASS** | 13.1s |
| admet-aspirin | cross-package | pass | pass | pass | pass | **PASS** | 9.9s |
| admet-table-full-panel | cross-package | pass | pass | pass | pass | **PASS** | 9.0s |
| reinvent-analogs | cross-package | pass | pass | pass | pass | **PASS** | 23.6s |
| docking-autodock | cross-package | pass | pass | pass | pass | **PASS** | 15.8s |
| boltz1-docking | cross-package | pass | fail | fail | pass | **FAIL** | 6.3s |
| chembl-id-then-properties | identifier-resolution | pass | pass | pass | pass | **PASS** | 17.6s |
| drugname-to-smiles-then-mw | identifier-resolution | pass | fail | fail | pass | **FAIL** | 13.9s |
| wrong-tool-admet-as-chem | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 7.4s |
| wrong-tool-generate-as-chem | wrong-tool-trap | pass | fail | fail | pass | **FAIL** | 11.8s |
| wrong-tool-dock-as-chem | wrong-tool-trap | pass | fail | fail | pass | **FAIL** | 6.4s |
| footgun-logs-property | footgun | pass | pass | fail | pass | **FAIL** | 15.3s |
| footgun-chemprop-train | footgun | pass | pass | pass | pass | **PASS** | 24.0s |

## Failed prompts (deep)

### boltz1-docking — cross-package

**Prompt:** Predict a binding pose for this ligand using Boltz-1.
**Intent:** Boltz-1 routing. The user asked for Boltz specifically; must pick Boltz1:docking (small-mol) NOT Boltz1:folding (proteins) and NOT AutoDock (different tool).

**Rubric:**

- path: PASS
- helpers: FAIL (`Boltz1:docking|call_function.*Boltz1:docking`)
- symbols: FAIL (`ligands|ligand`, `config`)
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`Boltz1:docking|call_function.*Boltz1:docking`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`ligands|ligand`, `config`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
// Check what columns are available
const cols = [...t.columns].map(c => `${c.name} (${c.type}, semType: ${c.semType})`);
const molCol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);
return ui.tableFromMap({
  'Columns': cols.join('\n'),
  'Molecule column': molCol ? molCol.name : 'None found'
});
```

<details><summary>Full response text</summary>

```
I need to find the ligand to dock. Let me check the current table for a molecule column.

```datagrok-exec
// Check what columns are available
const cols = [...t.columns].map(c => `${c.name} (${c.type}, semType: ${c.semType})`);
const molCol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);
return ui.tableFromMap({
  'Columns': cols.join('\n'),
  'Molecule column': molCol ? molCol.name : 'None found'
});
```
```

</details>

### drugname-to-smiles-then-mw — identifier-resolution

**Prompt:** Compute MW and LogP for ibuprofen.
**Intent:** Drug name → SMILES via Chem:namesToSmiles (ChEMBL-backed), then properties. namesToSmiles is preferred over freeTextToSmiles (LLM fallback) for clean named compounds.

**Rubric:**

- path: PASS
- helpers: FAIL (`(Chem|Chembl):namesToSmiles|call_function.*(Chem|Chembl):namesToSmiles`)
- symbols: FAIL (`ibuprofen`)
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`(Chem|Chembl):namesToSmiles|call_function.*(Chem|Chembl):namesToSmiles`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`ibuprofen`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
// Ibuprofen canonical SMILES (well-known structure)
const smiles = 'CC(C)Cc1ccc(cc1)C(C)C(=O)O';

// Build a one-row DataFrame
const df = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('molecule', [smiles])
]);
df.col('molecule').semType = DG.SEMTYPE.MOLECULE;

// Compute MW and LogP
const props = await grok.functions.call('Chem:getProperties', {
  molecules: df.col('molecule'),
  selected: ['MW', 'LogP']
});

// Build a result table
const resultDf = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('Compound', ['Ibuprofen']),
  DG.Column.fromStrings('SMILES', [smiles]),
  DG.Column.fromFloat32Array('MW', new Float32Array([props.col('MW').get(0)])),
  DG.Column.fromFloat32Array('LogP', new Float32Array([props.col('LogP').get(0)])),
]);

return DG.Viewer.grid(resultDf).root;
```

<details><summary>Full response text</summary>

```
I'll look up ibuprofen's SMILES and compute MW + LogP in one block — first converting the name to SMILES, then computing properties.

```datagrok-exec
// Ibuprofen canonical SMILES (well-known structure)
const smiles = 'CC(C)Cc1ccc(cc1)C(C)C(=O)O';

// Build a one-row DataFrame
const df = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('molecule', [smiles])
]);
df.col('molecule').semType = DG.SEMTYPE.MOLECULE;

// Compute MW and LogP
const props = await grok.functions.call('Chem:getProperties', {
  molecules: df.col('molecule'),
  selected: ['MW', 'LogP']
});

// Build a result table
const resultDf = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('Compound', ['Ibuprofen']),
  DG.Column.fromStrings('SMILES', [smiles]),
  DG.Column.fromFloat32Array('MW', new Float32Array([props.col('MW').get(0)])),
  DG.Column.fromFloat32Array('LogP', new Float32Array([props.col('LogP').get(0)])),
]);

return DG.Viewer.grid(resultDf).root;
```

This computes:
- **MW** — molecular weight (Da) via OpenChemLib
- **LogP** — Crippen/OCL partition coefficient

The result is displayed as an interactive grid. Ibuprofen (MW ≈ 206.3, LogP ≈ 3.5) is a classic lipophilic NSAID, well within drug-like space.
```

</details>

### wrong-tool-generate-as-chem — wrong-tool-trap

**Prompt:** Generate 20 new analogs of the lead molecule CCN(CC)CC. Use Reinvent4 with the 'optimize' run type.
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
// Step 1: List available Reinvent4 folders/targets
const folders = await grok.functions.call('Reinvent4:getFolders', {});
return ui.tableFromMap({ 'Available targets': folders.join(', ') });
```

<details><summary>Full response text</summary>

```
I'll first check the available Reinvent4 optimization targets, then generate the analogs.

```datagrok-exec
// Step 1: List available Reinvent4 folders/targets
const folders = await grok.functions.call('Reinvent4:getFolders', {});
return ui.tableFromMap({ 'Available targets': folders.join(', ') });
```
```

</details>

### wrong-tool-dock-as-chem — wrong-tool-trap

**Prompt:** Dock these ligands into the protein target.
**Intent:** Trap: 'docking' must route to Docking or Boltz1, not to Chem search/property functions.
**Trap:** Docking keyword → Docking package (AutoDock) or Boltz1 (deep learning), never Chem.

**Rubric:**

- path: PASS
- helpers: FAIL (`Docking:getAutodockResults|Docking:runAutodock|Boltz1:docking|call_function.*(Docking:(getAutodockResults|runAutodock)|Boltz1:docking)`)
- symbols: FAIL (`ligands`)
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`Docking:getAutodockResults|Docking:runAutodock|Boltz1:docking|call_function.*(Docking:(getAutodockResults|runAutodock)|Boltz1:docking)`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`ligands`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
// Find the molecule column and list docking targets in parallel
const molCol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);
const targets = await grok.functions.call('Docking:getConfigFiles', {});
return ui.tableFromMap({
  'Molecule column': molCol ? molCol.name : 'NOT FOUND',
  'Available targets': targets.length > 0 ? targets.join(', ') : 'NONE'
});
```

<details><summary>Full response text</summary>

```
I need to find the available docking targets and the molecule column in the current table. Let me do both at once.

```datagrok-exec
// Find the molecule column and list docking targets in parallel
const molCol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);
const targets = await grok.functions.call('Docking:getConfigFiles', {});
return ui.tableFromMap({
  'Molecule column': molCol ? molCol.name : 'NOT FOUND',
  'Available targets': targets.length > 0 ? targets.join(', ') : 'NONE'
});
```
```

</details>

### footgun-logs-property — footgun

**Prompt:** Compute logS (solubility) for these molecules.
**Intent:** logS is a parameter of addChemPropertiesColumns / getProperties, not a standalone function. Must NOT invent Chem:CalculateLogS (the catalog only exposes the property flag) and must NOT confuse with cLogP.
**Trap:** logS lives inside addChemPropertiesColumns({logS:true}) / getProperties(['logS']); there is no Chem:CalculateLogS.

**Rubric:**

- path: PASS
- helpers: PASS
- symbols: FAIL (`logS`)
- forbidden: PASS

**Hypothesis:** Expected symbol(s) absent (`logS`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
const mol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);
if (!mol) throw new Error('No Molecule column found in the current table.');

const res = await grok.functions.call('Chem:getProperties', {
  molecules: mol,
  selected: ['LogS']
});

// Append the LogS column to the table
for (let i = 0; i < res.columns.length; i++) {
  const col = res.columns.byIndex(i);
  if (!t.columns.contains(col.name))
    t.columns.add(col);
}
```

<details><summary>Full response text</summary>

```
I need to find the right function for computing logS. Let me check FUNCTIONS.md.Based on the chem-toolkit skill, **logS** is one of the chemical properties available via `Chem:addChemPropertiesColumns` and `Chem:getProperties`. I'll compute it for the molecules in your current table.

```datagrok-exec
const mol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);
if (!mol) throw new Error('No Molecule column found in the current table.');

const res = await grok.functions.call('Chem:getProperties', {
  molecules: mol,
  selected: ['LogS']
});

// Append the LogS column to the table
for (let i = 0; i < res.columns.length; i++) {
  const col = res.columns.byIndex(i);
  if (!t.columns.contains(col.name))
    t.columns.add(col);
}
```

This computes **LogS (aqueous solubility)** using the OpenChemLib model and appends the result as a new column to your table. The values represent log₁₀ of the molar solubility — higher (less negative) values indicate better solubility.
```

</details>

## Suggested skill improvements

- Expected helper `Boltz1:docking|call_function.*Boltz1:docking` not reached for: boltz1-docking. Consider tightening the skill's routing table or example density for this case.
- Expected helper `(Chem|Chembl):namesToSmiles|call_function.*(Chem|Chembl):namesToSmiles` not reached for: drugname-to-smiles-then-mw. Consider tightening the skill's routing table or example density for this case.
- Expected helper `Reinvent4:reinvent|Reinvent4:reinventTopMenu|call_function.*Reinvent4:(reinvent|reinventTopMenu)` not reached for: wrong-tool-generate-as-chem. Consider tightening the skill's routing table or example density for this case.
- Expected helper `Docking:getAutodockResults|Docking:runAutodock|Boltz1:docking|call_function.*(Docking:(getAutodockResults|runAutodock)|Boltz1:docking)` not reached for: wrong-tool-dock-as-chem. Consider tightening the skill's routing table or example density for this case.

## Anomalies

- None.
