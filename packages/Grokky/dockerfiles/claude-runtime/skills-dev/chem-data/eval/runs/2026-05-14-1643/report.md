# Eval report — chem-data

Started: 2026-05-14T13:43:40.240Z
Duration: 49.3s
Model: `sonnet`
System prompt size: 28,817 chars
Inlined skills: `datagrok-exec`, `chem-data`

## Summary

- **Total prompts:** 17
- **Passed:** 11
- **Failed:** 6
- **Pass rate:** 64.7%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 5 | 6 | 83% |
| cross-package | 2 | 3 | 67% |
| footgun | 1 | 2 | 50% |
| identifier-resolution | 2 | 4 | 50% |
| wrong-tool-trap | 1 | 2 | 50% |

### Latency

- Total: 49.3s wall
- Mean: 11.0s · Median: 10.5s · p95: 23.6s
- Slowest 3: `inchikey-to-chembl` (23.6s), `batch-chembl-ids-to-smiles` (16.5s), `moltrack-register-compound` (14.0s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`chem-data`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| chembl-substructure-search | atomic-do | pass | pass | pass | pass | **PASS** | 9.5s |
| chembl-similarity-search | atomic-do | pass | pass | pass | pass | **PASS** | 7.9s |
| chembl-compounds-by-organism | atomic-do | pass | pass | pass | pass | **PASS** | 6.0s |
| moltrack-register-compound | atomic-do | fail | fail | fail | pass | **FAIL** | 14.0s |
| moltrack-lookup-corporate-id | atomic-do | pass | pass | pass | pass | **PASS** | 6.1s |
| curves-fit | atomic-do | pass | pass | pass | pass | **PASS** | 9.6s |
| chemblid-to-smiles | identifier-resolution | pass | pass | pass | pass | **PASS** | 13.3s |
| drugname-to-smiles | identifier-resolution | pass | fail | pass | pass | **FAIL** | 9.1s |
| inchikey-to-chembl | identifier-resolution | pass | pass | fail | pass | **FAIL** | 23.6s |
| batch-chembl-ids-to-smiles | identifier-resolution | pass | pass | pass | pass | **PASS** | 16.5s |
| curves-ic50-from-data | cross-package | pass | fail | fail | pass | **FAIL** | 4.8s |
| hittriage-start-campaign | cross-package | pass | pass | pass | pass | **PASS** | 6.5s |
| moltrack-advanced-search | cross-package | pass | pass | pass | pass | **PASS** | 11.7s |
| wrong-tool-thrombin-compounds | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 13.5s |
| wrong-tool-lookup-by-name | wrong-tool-trap | pass | fail | pass | pass | **FAIL** | 12.4s |
| footgun-chembl-namespace | footgun | pass | pass | pass | pass | **PASS** | 10.5s |
| footgun-names-to-smiles-depends-on-chembl | footgun | pass | fail | pass | pass | **FAIL** | 11.6s |

## Failed prompts (deep)

### moltrack-register-compound — atomic-do

**Prompt:** Register this compound into MolTrack with its properties payload.
**Intent:** Registration entry points: registerMolTrackProperties for property defs, registerBulk for CSV. Must NOT confuse with lookup or search functions.

**Rubric:**

- path: FAIL (`>=1 datagrok-exec block`)
- helpers: FAIL (`MolTrack:registerMolTrackProperties|MolTrack:registerBulk|call_function.*MolTrack:(registerMolTrackProperties|registerBulk)`)
- symbols: FAIL (`jsonPayload|payload|csvFile`)
- forbidden: PASS

**Hypothesis:** Claude did not emit a `datagrok-exec` block at all — likely treating the request as conversational. Skill may need a stronger "always emit a block" directive for this intent class. Claude bypassed the expected helper(s) (`MolTrack:registerMolTrackProperties|MolTrack:registerBulk|call_function.*MolTrack:(registerMolTrackProperties|registerBulk)`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`jsonPayload|payload|csvFile`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:** (none)

<details><summary>Full response text</summary>

```
I need a bit more information to register the compound correctly.It looks like the question prompt wasn't answered. To register a compound into MolTrack, I need:

1. **Structure** — the SMILES (or InChI / molblock) of the compound you want to register.
2. **Properties** *(optional)* — any additional fields you want to store, for example:
   - `name` / `corporate_id`
   - `purity`, `molecular_weight`
   - `project`, `scientist`, etc.

Please share those details and I'll emit the registration call right away. For example:

> *"Register aspirin (`CC(=O)Oc1ccccc1C(=O)O`) with project = 'Analgesics' and purity = 99.5"*
```

</details>

### drugname-to-smiles — identifier-resolution

**Prompt:** What is the SMILES of aspirin?
**Intent:** Free-text drug name → SMILES. Chem:namesToSmiles is ChEMBL-backed lookup for clean names. Must NOT use chemblIdToSmilesTs (only accepts CHEMBL[0-9]+).
**Trap:** namesToSmiles is a Chem function but backed by ChEMBL DB — depends on Chembl install.

**Rubric:**

- path: PASS
- helpers: FAIL (`Chem:namesToSmiles|call_function.*Chem:namesToSmiles`)
- symbols: PASS
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`Chem:namesToSmiles|call_function.*Chem:namesToSmiles`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper.

**Emitted exec blocks:**

```javascript
const result = await grok.functions.call('Chembl:namesToSmiles', { names: ['aspirin'] });
const smiles = result.columns.byName('canonical_smiles').get(0);
return ui.tableFromMap({ 'Compound': 'Aspirin', 'SMILES': smiles });
```

<details><summary>Full response text</summary>

```
```datagrok-exec
const result = await grok.functions.call('Chembl:namesToSmiles', { names: ['aspirin'] });
const smiles = result.columns.byName('canonical_smiles').get(0);
return ui.tableFromMap({ 'Compound': 'Aspirin', 'SMILES': smiles });
```

The SMILES for **aspirin** will be retrieved directly from the ChEMBL database and displayed above. The canonical SMILES for aspirin is typically:

**`CC(=O)Oc1ccccc1C(=O)O`**

This represents acetylsalicylic acid — an acetyl group (`CC(=O)O-`) attached to a benzene ring that also carries a carboxylic acid group (`-C(=O)O`).
```

</details>

### inchikey-to-chembl — identifier-resolution

**Prompt:** Get the ChEMBL ID for InChI key BSYNRYMUTXBXSQ-UHFFFAOYSA-N.
**Intent:** Cross-ID mapping → Chem:getMapIdentifiers / mapIdentifiersTransform (the catalog's ID-mapping entry). chemblIdToSmilesTs only goes ID→SMILES, not the reverse, so it's wrong here.

**Rubric:**

- path: PASS
- helpers: PASS
- symbols: FAIL (`toSource|to_source|CHEMBL|chembl`)
- forbidden: PASS

**Hypothesis:** Expected symbol(s) absent (`toSource|to_source|CHEMBL|chembl`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
// Use the UniChem/ChEMBL identifier mapping to resolve InChI key → ChEMBL ID
const df = DG.DataFrame.fromCsv('inchi_key\nBSYNRYMUTXBXSQ-UHFFFAOYSA-N');

const result = await grok.functions.call('Chem:mapIdentifiersTransform', {
  table: df,
  fromIdentifier: 'InChIKey',
  toIdentifier: 'ChEMBL ID',
  keyColumn: df.col('inchi_key'),
});

// Show the result table
return DG.Viewer.grid(df).root;
```

<details><summary>Full response text</summary>

```
The `mapIdentifiers` function from the `Chem` package handles InChI key → ChEMBL ID conversion. Let me look it up and run it:Based on the `datagrok-chem-data` skill, the route for **InChI key → ChEMBL ID** is `Chem:mapIdentifiersTransform`. I'll call it now:

```datagrok-exec
// Use the UniChem/ChEMBL identifier mapping to resolve InChI key → ChEMBL ID
const df = DG.DataFrame.fromCsv('inchi_key\nBSYNRYMUTXBXSQ-UHFFFAOYSA-N');

const result = await grok.functions.call('Chem:mapIdentifiersTransform', {
  table: df,
  fromIdentifier: 'InChIKey',
  toIdentifier: 'ChEMBL ID',
  keyColumn: df.col('inchi_key'),
});

// Show the result table
return DG.Viewer.grid(df).root;
```

This runs `Chem:mapIdentifiersTransform` to look up the ChEMBL ID for InChI key **BSYNRYMUTXBXSQ-UHFFFAOYSA-N** (which corresponds to **Aspirin / Acetylsalicylic acid**, expected ChEMBL ID: `CHEMBL25`). The result will be appended as a new column in the table above.

> **Note:** The InChI key `BSYNRYMUTXBXSQ-UHFFFAOYSA-N` is the standard key for **aspirin (acetylsalicylic acid)**. Its ChEMBL ID is **CHEMBL25**. The exec block above confirms this via the live platform lookup.
```

</details>

### curves-ic50-from-data — cross-package

**Prompt:** Compute IC50 values from this dose-response dataset.
**Intent:** IC50 lives in Curves. Two-step: dataToCurves to fit, then addStatisticsColumn(table, 'IC50', 'IC50', 0). Pass if Claude reaches either function with IC50 mentioned.

**Rubric:**

- path: PASS
- helpers: FAIL (`Curves:dataToCurves|Curves:addStatisticsColumn|call_function.*Curves:(dataToCurves|addStatisticsColumn)`)
- symbols: FAIL (`IC50|ic50`)
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`Curves:dataToCurves|Curves:addStatisticsColumn|call_function.*Curves:(dataToCurves|addStatisticsColumn)`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`IC50|ic50`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
// Inspect available columns
return ui.tableFromMap(
  Object.fromEntries(t.columns.toList().map(c => [c.name, c.type]))
);
```

<details><summary>Full response text</summary>

```
I'll help you compute IC50 values from your dose-response dataset. First, let me check the current table columns to identify the right columns for concentration, readout, and compound IDs.

```datagrok-exec
// Inspect available columns
return ui.tableFromMap(
  Object.fromEntries(t.columns.toList().map(c => [c.name, c.type]))
);
```
```

</details>

### wrong-tool-lookup-by-name — wrong-tool-trap

**Prompt:** Look up acetylsalicylic acid in the database.
**Intent:** Drug-name DB lookup → Chem:namesToSmiles (ChEMBL-backed). Must NOT freeform LLM-parse or treat as a SMILES validation problem.
**Trap:** 'acetylsalicylic acid' is a chemical name — needs name lookup, not validation/parsing.

**Rubric:**

- path: PASS
- helpers: FAIL (`Chem:namesToSmiles|Chembl:chemblSimilaritySearch|call_function.*(Chem:namesToSmiles|Chembl:)`)
- symbols: PASS
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`Chem:namesToSmiles|Chembl:chemblSimilaritySearch|call_function.*(Chem:namesToSmiles|Chembl:)`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper.

**Emitted exec blocks:**

```javascript
const result = await grok.functions.call('Chembl:namesToSmiles', { names: ['acetylsalicylic acid'] });
return DG.Viewer.grid(result).root;
```

<details><summary>Full response text</summary>

```
I'll look up acetylsalicylic acid (aspirin) by resolving its name to a SMILES string via the ChEMBL name-to-SMILES lookup:

```datagrok-exec
const result = await grok.functions.call('Chembl:namesToSmiles', { names: ['acetylsalicylic acid'] });
return DG.Viewer.grid(result).root;
```

This queries the local ChEMBL database to resolve **acetylsalicylic acid** (aspirin) to its canonical SMILES. The result table will show the matched SMILES string from ChEMBL. If you'd like to go further — such as running a similarity/substructure search or pulling bioactivity data — just let me know!
```

</details>

### footgun-names-to-smiles-depends-on-chembl — footgun

**Prompt:** Convert this list of drug names ['aspirin', 'ibuprofen', 'caffeine'] to SMILES.
**Intent:** namesToSmiles needs the Chembl DB installed — Claude should pick it (catalog says so) without surfacing infra concerns unprompted, but must NOT fall back to freeTextToSmiles which is the LLM path.
**Trap:** namesToSmiles is exposed via Chem but backed by ChEMBL. Catalog doesn't make the dependency explicit — borderline ambiguous wording.

**Rubric:**

- path: PASS
- helpers: FAIL (`Chem:namesToSmiles|call_function.*Chem:namesToSmiles`)
- symbols: PASS
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`Chem:namesToSmiles|call_function.*Chem:namesToSmiles`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper.

**Emitted exec blocks:**

```javascript
const df = await grok.functions.call('Chembl:namesToSmiles', {
  names: ['aspirin', 'ibuprofen', 'caffeine']
});
return DG.Viewer.grid(df).root;
```

<details><summary>Full response text</summary>

```
Based on the skill documentation, I can use `Chembl:namesToSmiles` — a SQL query that converts drug names to SMILES. Let me call it now:

```datagrok-exec
const df = await grok.functions.call('Chembl:namesToSmiles', {
  names: ['aspirin', 'ibuprofen', 'caffeine']
});
return DG.Viewer.grid(df).root;
```

This calls the `Chembl:namesToSmiles` query, which looks up each name in the local ChEMBL database and returns a DataFrame with the corresponding `canonical_smiles`. The result is displayed as a grid above.

If any name isn't found, it will appear with a null/empty SMILES — in that case, double-check the spelling or try alternate names (e.g. `"acetylsalicylic acid"` for aspirin).
```

</details>

## Suggested skill improvements

- Expected helper `MolTrack:registerMolTrackProperties|MolTrack:registerBulk|call_function.*MolTrack:(registerMolTrackProperties|registerBulk)` not reached for: moltrack-register-compound. Consider tightening the skill's routing table or example density for this case.
- Expected helper `Chem:namesToSmiles|call_function.*Chem:namesToSmiles` not reached for: drugname-to-smiles, footgun-names-to-smiles-depends-on-chembl. Consider tightening the skill's routing table or example density for this case.
- Expected helper `Curves:dataToCurves|Curves:addStatisticsColumn|call_function.*Curves:(dataToCurves|addStatisticsColumn)` not reached for: curves-ic50-from-data. Consider tightening the skill's routing table or example density for this case.
- Expected helper `Chem:namesToSmiles|Chembl:chemblSimilaritySearch|call_function.*(Chem:namesToSmiles|Chembl:)` not reached for: wrong-tool-lookup-by-name. Consider tightening the skill's routing table or example density for this case.

## Anomalies

- None.
