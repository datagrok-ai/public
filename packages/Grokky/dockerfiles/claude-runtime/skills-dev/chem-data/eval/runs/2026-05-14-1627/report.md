# Eval report — chem-data

Started: 2026-05-14T13:27:34.705Z
Duration: 67.9s
Model: `sonnet`
System prompt size: 24,712 chars
Inlined skills: `datagrok-exec`, `chem-data`

## Summary

- **Total prompts:** 17
- **Passed:** 12
- **Failed:** 5
- **Pass rate:** 70.6%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 5 | 6 | 83% |
| cross-package | 2 | 3 | 67% |
| footgun | 1 | 2 | 50% |
| identifier-resolution | 2 | 4 | 50% |
| wrong-tool-trap | 2 | 2 | 100% |

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`chem-data`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| chembl-substructure-search | atomic-do | pass | pass | pass | pass | **PASS** | 13.3s |
| chembl-similarity-search | atomic-do | pass | pass | pass | pass | **PASS** | 12.8s |
| chembl-compounds-by-organism | atomic-do | pass | pass | pass | pass | **PASS** | 13.6s |
| moltrack-register-compound | atomic-do | fail | fail | fail | pass | **FAIL** | 13.0s |
| moltrack-lookup-corporate-id | atomic-do | pass | pass | pass | pass | **PASS** | 11.9s |
| curves-fit | atomic-do | pass | pass | pass | pass | **PASS** | 14.0s |
| chemblid-to-smiles | identifier-resolution | pass | pass | pass | pass | **PASS** | 7.3s |
| drugname-to-smiles | identifier-resolution | fail | fail | fail | pass | **FAIL** | 4.1s |
| inchikey-to-chembl | identifier-resolution | pass | fail | pass | pass | **FAIL** | 31.4s |
| batch-chembl-ids-to-smiles | identifier-resolution | pass | pass | pass | pass | **PASS** | 14.7s |
| curves-ic50-from-data | cross-package | fail | fail | fail | pass | **FAIL** | 5.2s |
| hittriage-start-campaign | cross-package | pass | pass | pass | pass | **PASS** | 5.9s |
| moltrack-advanced-search | cross-package | pass | pass | pass | pass | **PASS** | 14.2s |
| wrong-tool-thrombin-compounds | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 35.1s |
| wrong-tool-lookup-by-name | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 15.7s |
| footgun-chembl-namespace | footgun | pass | pass | pass | pass | **PASS** | 10.5s |
| footgun-names-to-smiles-depends-on-chembl | footgun | pass | fail | fail | pass | **FAIL** | 15.1s |

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
I need to ask you for the details before registering the compound.It looks like the question prompt was dismissed. To register your compound into MolTrack, I need a few details:

1. **Structure** — Please provide the **SMILES string** (or molblock) of the compound you want to register.
2. **Properties payload** — Provide the JSON properties you want to attach. For example:
   ```json
   {
     "canonical_smiles": "CCO",
     "name": "Ethanol",
     "purity": 99.5,
     "source": "Sigma-Aldrich"
   }
   ```
3. **Scope** — Should this be registered as a **compound**, a **batch**, or both?

Once you share these details, I'll call `MolTrack:registerMolTrackProperties` (for property definitions) or `MolTrack:registerBulk` (for bulk/single-row CSV registration) right away.
```

</details>

### drugname-to-smiles — identifier-resolution

**Prompt:** What is the SMILES of aspirin?
**Intent:** Free-text drug name → SMILES. Chem:namesToSmiles is ChEMBL-backed lookup for clean names. Must NOT use chemblIdToSmilesTs (only accepts CHEMBL[0-9]+).
**Trap:** namesToSmiles is a Chem function but backed by ChEMBL DB — depends on Chembl install.

**Rubric:**

- path: FAIL (`>=1 datagrok-exec block`)
- helpers: FAIL (`Chem:namesToSmiles|call_function.*Chem:namesToSmiles`)
- symbols: FAIL (`aspirin`, `names`)
- forbidden: PASS

**Hypothesis:** Claude did not emit a `datagrok-exec` block at all — likely treating the request as conversational. Skill may need a stronger "always emit a block" directive for this intent class. Claude bypassed the expected helper(s) (`Chem:namesToSmiles|call_function.*Chem:namesToSmiles`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`aspirin`, `names`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:** (none)

<details><summary>Full response text</summary>

```
The SMILES of aspirin (acetylsalicylic acid) is:

**`CC(=O)Oc1ccccc1C(=O)O`**

This represents:
- `CC(=O)O` — the acetyl ester group
- `c1ccccc1` — the benzene ring
- `C(=O)O` — the carboxylic acid group
```

</details>

### inchikey-to-chembl — identifier-resolution

**Prompt:** Get the ChEMBL ID for InChI key BSYNRYMUTXBXSQ-UHFFFAOYSA-N.
**Intent:** Cross-ID mapping → Chem:getMapIdentifiers / mapIdentifiersTransform (the catalog's ID-mapping entry). chemblIdToSmilesTs only goes ID→SMILES, not the reverse, so it's wrong here.

**Rubric:**

- path: PASS
- helpers: FAIL (`Chem:getMapIdentifiers|Chem:mapIdentifiersTransform|call_function.*Chem:(getMapIdentifiers|mapIdentifiersTransform)`)
- symbols: PASS
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`Chem:getMapIdentifiers|Chem:mapIdentifiersTransform|call_function.*Chem:(getMapIdentifiers|mapIdentifiersTransform)`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper.

**Emitted exec blocks:**

```javascript
const df = await grok.data.query('Chembl:patternSubstructureSearch', { pattern: 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N' });
return DG.Viewer.grid(df).root;
```

```javascript
const df = await grok.data.query('Chembl:patternSubstructureSearch', { pattern: 'CC(=O)Oc1ccccc1C(=O)O' });
return DG.Viewer.grid(df).root;
```

```javascript
const df = await grok.data.query('System:ExecuteQuery', {
  connection: 'Chembl:Chembl',
  query: `SELECT cs.chembl_id, md.standard_inchi_key 
          FROM compound_structures cs 
          JOIN molecule_dictionary md ON cs.molregno = md.molregno 
          WHERE cs.standard_inchi_key = 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'`
});
return DG.Viewer.grid(df).root;
```

<details><summary>Full response text</summary>

```
The InChI key `BSYNRYMUTXBXSQ-UHFFFAOYSA-N` corresponds to **Aspirin (acetylsalicylic acid)**. Let me look it up via the ChEMBL database using the UniChem/ChEMBL web service.Let me query this directly within the Datagrok platform using the ChEMBL database:

```datagrok-exec
const df = await grok.data.query('Chembl:patternSubstructureSearch', { pattern: 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N' });
return DG.Viewer.grid(df).root;
```

Let me try a more appropriate query using the ChEMBL SQL queries available:The Datagrok workspace isn't accessible from this evaluation harness directory. However, based on well-known chemical databases:

**The ChEMBL ID for InChI key `BSYNRYMUTXBXSQ-UHFFFAOYSA-N` is `CHEMBL25`.**

This corresponds to **Aspirin (acetylsalicylic acid)**, one of the most well-known drugs, which is compound #25 in ChEMBL.

To look this up programmatically within Datagrok, you can run the following against the local ChEMBL PostgreSQL connection:

```datagrok-exec
const df = await grok.data.query('Chembl:patternSubstructureSearch', { pattern: 'CC(=O)Oc1ccccc1C(=O)O' });
return DG.Viewer.grid(df).root;
```

Or, for a direct InChI key lookup via SQL:

```datagrok-exec
const df = await grok.data.query('System:ExecuteQuery', {
  connection: 'Chembl:Chembl',
  query: `SELECT cs.chembl_id, md.standard_inchi_key 
          FROM compound_structures cs 
          JOIN molecule_dictionary md ON cs.molregno = md.molregno 
          WHERE cs.standard_inchi_key = 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'`
});
return DG.Viewer.grid(df).root;
```

**Answer: `CHEMBL25` (Aspirin / Acetylsalicylic acid)**
```

</details>

### curves-ic50-from-data — cross-package

**Prompt:** Compute IC50 values from this dose-response dataset.
**Intent:** IC50 lives in Curves. Two-step: dataToCurves to fit, then addStatisticsColumn(table, 'IC50', 'IC50', 0). Pass if Claude reaches either function with IC50 mentioned.

**Rubric:**

- path: FAIL (`>=1 datagrok-exec block`)
- helpers: FAIL (`Curves:dataToCurves|Curves:addStatisticsColumn|call_function.*Curves:(dataToCurves|addStatisticsColumn)`)
- symbols: FAIL (`IC50|ic50`)
- forbidden: PASS

**Hypothesis:** Claude did not emit a `datagrok-exec` block at all — likely treating the request as conversational. Skill may need a stronger "always emit a block" directive for this intent class. Claude bypassed the expected helper(s) (`Curves:dataToCurves|Curves:addStatisticsColumn|call_function.*Curves:(dataToCurves|addStatisticsColumn)`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`IC50|ic50`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:** (none)

<details><summary>Full response text</summary>

```
I'd be happy to help compute IC50 values from your dose-response dataset! However, I don't see any dataset attached to your message.

Could you please share your data? You can:

1. **Upload a file** (CSV, Excel, etc.) with your dose-response data
2. **Paste the data** directly in the chat

Also, to set up the analysis correctly, please let me know:
- Which column contains the **concentration** values?
- Which column contains the **readout/response** values?
- Which columns identify **compound**, **assay**, and **run** (if applicable)?
```

</details>

### footgun-names-to-smiles-depends-on-chembl — footgun

**Prompt:** Convert this list of drug names ['aspirin', 'ibuprofen', 'caffeine'] to SMILES.
**Intent:** namesToSmiles needs the Chembl DB installed — Claude should pick it (catalog says so) without surfacing infra concerns unprompted, but must NOT fall back to freeTextToSmiles which is the LLM path.
**Trap:** namesToSmiles is exposed via Chem but backed by ChEMBL. Catalog doesn't make the dependency explicit — borderline ambiguous wording.

**Rubric:**

- path: PASS
- helpers: FAIL (`Chem:namesToSmiles|call_function.*Chem:namesToSmiles`)
- symbols: FAIL (`names`)
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`Chem:namesToSmiles|call_function.*Chem:namesToSmiles`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`names`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
const drugs = [
  { name: 'aspirin',    smiles: 'CC(=O)Oc1ccccc1C(=O)O' },
  { name: 'ibuprofen',  smiles: 'CC(C)Cc1ccc(cc1)C(C)C(=O)O' },
  { name: 'caffeine',   smiles: 'Cn1cnc2c1c(=O)n(c(=O)n2C)C' },
];

const df = DG.DataFrame.fromObjects(drugs);
df.name = 'Drug SMILES';

// Mark the smiles column with the correct semantic type so structures render
df.col('smiles').semType = DG.SEMTYPE.MOLECULE;

grok.shell.addTableView(df);
```

<details><summary>Full response text</summary>

```
No registered function is available for drug name → SMILES conversion, but this is a well-known mapping I can resolve directly using standard SMILES notation. I'll create a DataFrame with the drug names and their corresponding SMILES, then display it.

```datagrok-exec
const drugs = [
  { name: 'aspirin',    smiles: 'CC(=O)Oc1ccccc1C(=O)O' },
  { name: 'ibuprofen',  smiles: 'CC(C)Cc1ccc(cc1)C(C)C(=O)O' },
  { name: 'caffeine',   smiles: 'Cn1cnc2c1c(=O)n(c(=O)n2C)C' },
];

const df = DG.DataFrame.fromObjects(drugs);
df.name = 'Drug SMILES';

// Mark the smiles column with the correct semantic type so structures render
df.col('smiles').semType = DG.SEMTYPE.MOLECULE;

grok.shell.addTableView(df);
```

Here are the SMILES strings for your three compounds:

| Drug | SMILES |
|---|---|
| **Aspirin** | `CC(=O)Oc1ccccc1C(=O)O` |
| **Ibuprofen** | `CC(C)Cc1ccc(cc1)C(C)C(=O)O` |
| **Caffeine** | `Cn1cnc2c1c(=O)n(c(=O)n2C)C` |

The table has been opened in Datagrok with the `smiles` column tagged as `Molecule` semantic type, so structures will render as 2D depictions in the grid. Let me know if you'd like to run any cheminformatics analyses on these compounds!
```

</details>

## Suggested skill improvements

- Expected helper `MolTrack:registerMolTrackProperties|MolTrack:registerBulk|call_function.*MolTrack:(registerMolTrackProperties|registerBulk)` not reached for: moltrack-register-compound. Consider tightening the skill's routing table or example density for this case.
- Expected helper `Chem:namesToSmiles|call_function.*Chem:namesToSmiles` not reached for: drugname-to-smiles, footgun-names-to-smiles-depends-on-chembl. Consider tightening the skill's routing table or example density for this case.
- Expected helper `Chem:getMapIdentifiers|Chem:mapIdentifiersTransform|call_function.*Chem:(getMapIdentifiers|mapIdentifiersTransform)` not reached for: inchikey-to-chembl. Consider tightening the skill's routing table or example density for this case.
- Expected helper `Curves:dataToCurves|Curves:addStatisticsColumn|call_function.*Curves:(dataToCurves|addStatisticsColumn)` not reached for: curves-ic50-from-data. Consider tightening the skill's routing table or example density for this case.

## Anomalies

- `wrong-tool-thrombin-compounds`: emitted 4 exec blocks — unusually verbose; check whether the skill should encourage consolidation.
