# Eval report — chem-toolkit

Started: 2026-05-14T13:25:58.064Z
Duration: 96.1s
Model: `sonnet`
System prompt size: 37,251 chars
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

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`chem-toolkit`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| lipinski-properties | atomic-do | pass | pass | pass | pass | **PASS** | 12.7s |
| smiles-to-inchi | atomic-do | pass | pass | pass | pass | **PASS** | 10.2s |
| convert-notation-molblock | atomic-do | pass | pass | pass | pass | **PASS** | 15.6s |
| chemical-space-umap | atomic-do | pass | pass | pass | pass | **PASS** | 10.5s |
| similarity-search | atomic-do | pass | pass | pass | pass | **PASS** | 14.3s |
| substructure-search-benzene | atomic-do | pass | pass | pass | pass | **PASS** | 12.0s |
| rgroup-decomposition | atomic-do | pass | pass | pass | pass | **PASS** | 14.4s |
| admet-aspirin | cross-package | pass | pass | pass | pass | **PASS** | 16.7s |
| admet-table-full-panel | cross-package | pass | pass | pass | pass | **PASS** | 9.9s |
| reinvent-analogs | cross-package | pass | pass | pass | pass | **PASS** | 19.3s |
| docking-autodock | cross-package | pass | pass | pass | pass | **PASS** | 15.3s |
| boltz1-docking | cross-package | pass | fail | fail | pass | **FAIL** | 18.2s |
| chembl-id-then-properties | identifier-resolution | pass | pass | pass | pass | **PASS** | 22.3s |
| drugname-to-smiles-then-mw | identifier-resolution | pass | fail | fail | pass | **FAIL** | 13.4s |
| wrong-tool-admet-as-chem | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 24.5s |
| wrong-tool-generate-as-chem | wrong-tool-trap | fail | fail | fail | pass | **FAIL** | 7.5s |
| wrong-tool-dock-as-chem | wrong-tool-trap | pass | fail | fail | pass | **FAIL** | 13.5s |
| footgun-logs-property | footgun | pass | pass | pass | pass | **PASS** | 15.1s |
| footgun-chemprop-train | footgun | pass | pass | fail | pass | **FAIL** | 38.0s |

## Failed prompts (deep)

### boltz1-docking — cross-package

**Prompt:** Predict a binding pose for this ligand using Boltz-1.
**Intent:** Boltz-1 routing. The user asked for Boltz specifically; must pick Boltz1:docking (small-mol) NOT Boltz1:folding (proteins) and NOT AutoDock (different tool).

**Rubric:**

- path: PASS
- helpers: FAIL (`Boltz1:docking|call_function.*Boltz1:docking`)
- symbols: FAIL (`ligands|ligand`)
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`Boltz1:docking|call_function.*Boltz1:docking`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`ligands|ligand`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
// Get available Boltz configs and current table info
const configs = await grok.functions.call('Boltz1:getBoltzConfigFolders', {});
const hasMolCol = (view.type === 'TableView') 
  ? t.columns.toList().filter(c => c.semType === DG.SEMTYPE.MOLECULE) 
  : [];

return ui.tableFromMap({
  'Available Boltz configs': configs?.toString() ?? '(none)',
  'Molecule columns in table': hasMolCol.map(c => c.name).join(', ') || '(none / no table open)',
});
```

<details><summary>Full response text</summary>

```
I need to know which ligand you'd like to dock. Could you clarify?Let me first check what Boltz-1 config folders are available and what's in the current table.```datagrok-exec
// Get available Boltz configs and current table info
const configs = await grok.functions.call('Boltz1:getBoltzConfigFolders', {});
const hasMolCol = (view.type === 'TableView') 
  ? t.columns.toList().filter(c => c.semType === DG.SEMTYPE.MOLECULE) 
  : [];

return ui.tableFromMap({
  'Available Boltz configs': configs?.toString() ?? '(none)',
  'Molecule columns in table': hasMolCol.map(c => c.name).join(', ') || '(none / no table open)',
});
```

Once I see what configs and molecule columns are available, I can proceed with the docking. Could you also let me know:

1. **Which ligand** to use — the molecule column in the current table, or a specific SMILES string?
2. **Which Boltz-1 config** (target) to dock against (I'll show the available ones above)?
```

</details>

### drugname-to-smiles-then-mw — identifier-resolution

**Prompt:** Compute MW and LogP for ibuprofen.
**Intent:** Drug name → SMILES via Chem:namesToSmiles (ChEMBL-backed), then properties. namesToSmiles is preferred over freeTextToSmiles (LLM fallback) for clean named compounds.

**Rubric:**

- path: PASS
- helpers: FAIL (`Chem:namesToSmiles|call_function.*Chem:namesToSmiles`)
- symbols: FAIL (`ibuprofen`)
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`Chem:namesToSmiles|call_function.*Chem:namesToSmiles`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`ibuprofen`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
const smiles = 'CC(C)Cc1ccc(cc1)C(C)C(=O)O'; // Ibuprofen

// Create a small dataframe with one molecule
const df = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('molecule', [smiles])
]);
df.col('molecule').semType = DG.SEMTYPE.MOLECULE;

// Calculate properties
const resDf = await grok.functions.call('Chem:getProperties', {
  molecules: df.col('molecule'),
  selected: ['MW', 'LogP']
});

// Build a display map
const map = { 'Molecule': smiles };
for (const col of resDf.columns.toList()) {
  map[col.name] = col.get(0);
}

return ui.tableFromMap(map);
```

<details><summary>Full response text</summary>

```
I'll use the Chem package to compute MW and LogP for ibuprofen (SMILES: `CC(C)Cc1ccc(cc1)C(C)C(=O)O`) directly:

```datagrok-exec
const smiles = 'CC(C)Cc1ccc(cc1)C(C)C(=O)O'; // Ibuprofen

// Create a small dataframe with one molecule
const df = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('molecule', [smiles])
]);
df.col('molecule').semType = DG.SEMTYPE.MOLECULE;

// Calculate properties
const resDf = await grok.functions.call('Chem:getProperties', {
  molecules: df.col('molecule'),
  selected: ['MW', 'LogP']
});

// Build a display map
const map = { 'Molecule': smiles };
for (const col of resDf.columns.toList()) {
  map[col.name] = col.get(0);
}

return ui.tableFromMap(map);
```

This will compute **MW** (molecular weight) and **LogP** for ibuprofen using the RDKit/OCL-backed `Chem:getProperties` function and display the results inline.
```

</details>

### wrong-tool-generate-as-chem — wrong-tool-trap

**Prompt:** Generate 20 new analogs of this lead molecule.
**Intent:** Trap: 'generate analogs' is generative chemistry (Reinvent4), not reaction enumeration or MMP / R-group analysis (which describe existing series, not invent new ones).
**Trap:** Generation = Reinvent4. Reaction enumeration / MMP / R-group operate on existing molecules.

**Rubric:**

- path: FAIL (`>=1 datagrok-exec block`)
- helpers: FAIL (`Reinvent4:reinvent|Reinvent4:reinventTopMenu|call_function.*Reinvent4:(reinvent|reinventTopMenu)`)
- symbols: FAIL (`ligand`)
- forbidden: PASS

**Hypothesis:** Claude did not emit a `datagrok-exec` block at all — likely treating the request as conversational. Skill may need a stronger "always emit a block" directive for this intent class. Claude bypassed the expected helper(s) (`Reinvent4:reinvent|Reinvent4:reinventTopMenu|call_function.*Reinvent4:(reinvent|reinventTopMenu)`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`ligand`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:** (none)

<details><summary>Full response text</summary>

```
I need to know which molecule you'd like to generate analogs for. Could you please provide the lead molecule?Please provide the SMILES string for your lead molecule, or select the molecule in the grid and let me know — I'll then use Reinvent4 to generate 20 analogs optimized for your target.
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
const targets = await grok.functions.call('Docking:getConfigFiles');
return ui.tableFromMap({ 'Available targets': targets.join('\n') });
```

<details><summary>Full response text</summary>

```
I'd be happy to help dock ligands into a protein target! To proceed, I need a few clarifications.I need a couple of details to set up the docking run:

1. **Docking engine** — AutoDock Vina (classic, fast) or Boltz-1 (deep learning)?
2. **Protein target** — Which target should the ligands be docked into? I can list the available configured targets if you'd like.
3. **Ligand column** — Which column in the current table contains the ligand molecules?
4. **Number of poses** — How many docking poses per ligand would you like (default is usually 1–3)?

Let me fetch the available targets so you can choose:

```datagrok-exec
const targets = await grok.functions.call('Docking:getConfigFiles');
return ui.tableFromMap({ 'Available targets': targets.join('\n') });
```

Please share these details and I'll kick off the docking run right away!
```

</details>

### footgun-chemprop-train — footgun

**Prompt:** Train a Chemprop regression model on this dataset to predict the 'activity' column.
**Intent:** Chemprop training has a 25-arg snake_case signature. predictColumn + dataset_type='regression' + a sensible epochs value should appear. Must NOT confuse with applyChemprop (inference) or isApplicableNN (framework hook).
**Trap:** trainChemprop's signature is unusually long and snake_case; Claude often forgets dataset_type or mistakes apply for train.

**Rubric:**

- path: PASS
- helpers: PASS
- symbols: FAIL (`epochs`)
- forbidden: PASS

**Hypothesis:** Expected symbol(s) absent (`epochs`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
// Find the molecule column
const molCol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);
if (!molCol) throw new Error('No molecule column found in the current table.');

const activityCol = t.col('activity');
if (!activityCol) throw new Error('No "activity" column found in the current table.');

grok.shell.info('Training Chemprop regression model — this may take a few minutes...');

const model = await grok.functions.call('Chem:trainChemprop', {
  df: t,
  predictColumn: activityCol,
  dataset_type: 'regression',
});

grok.shell.info('Chemprop model trained successfully!');
return ui.tableFromMap({ 'Status': 'Training complete', 'Model type': 'Chemprop regression', 'Target column': 'activity', 'Molecule column': molCol.name });
```

<details><summary>Full response text</summary>

```
I'll train a Chemprop regression model using the current table's molecule column to predict `activity`. Let me run this now:

```datagrok-exec
// Find the molecule column
const molCol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);
if (!molCol) throw new Error('No molecule column found in the current table.');

const activityCol = t.col('activity');
if (!activityCol) throw new Error('No "activity" column found in the current table.');

grok.shell.info('Training Chemprop regression model — this may take a few minutes...');

const model = await grok.functions.call('Chem:trainChemprop', {
  df: t,
  predictColumn: activityCol,
  dataset_type: 'regression',
});

grok.shell.info('Chemprop model trained successfully!');
return ui.tableFromMap({ 'Status': 'Training complete', 'Model type': 'Chemprop regression', 'Target column': 'activity', 'Molecule column': molCol.name });
```

This will:
1. **Detect** the molecule column automatically (by semantic type `Molecule`)
2. **Train** a Chemprop Message Passing Neural Network (MPNN) regression model on `activity`
3. **Return** the trained model blob (usable later with `Chem:applyChemprop`)

Training runs inside the `chem-chemprop` Docker container and may take a few minutes depending on dataset size. Let me know if you'd like to adjust any hyperparameters (epochs, depth, dropout, learning rate, etc.) or apply the model to new molecules afterwards!
```

</details>

## Suggested skill improvements

- Expected helper `Boltz1:docking|call_function.*Boltz1:docking` not reached for: boltz1-docking. Consider tightening the skill's routing table or example density for this case.
- Expected helper `Chem:namesToSmiles|call_function.*Chem:namesToSmiles` not reached for: drugname-to-smiles-then-mw. Consider tightening the skill's routing table or example density for this case.
- Expected helper `Reinvent4:reinvent|Reinvent4:reinventTopMenu|call_function.*Reinvent4:(reinvent|reinventTopMenu)` not reached for: wrong-tool-generate-as-chem. Consider tightening the skill's routing table or example density for this case.
- Expected helper `Docking:getAutodockResults|Docking:runAutodock|Boltz1:docking|call_function.*(Docking:(getAutodockResults|runAutodock)|Boltz1:docking)` not reached for: wrong-tool-dock-as-chem. Consider tightening the skill's routing table or example density for this case.

## Anomalies

- None.
