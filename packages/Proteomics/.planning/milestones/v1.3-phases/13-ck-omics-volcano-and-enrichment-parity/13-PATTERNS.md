# Phase 13: CK-omics Volcano and Enrichment Parity - Pattern Map

**Mapped:** 2026-05-16
**Files analyzed:** 9 (1 new, 8 modified)
**Analogs found:** 9 / 9 (all analogs are the integration-surface files themselves — extend in place, do not rebuild)

> The integration surface for this phase is its own analog set: every file to
> touch already exists and is test-covered. This phase is *parameterizing
> existing functions* and *porting two Python algorithms verbatim*, not building
> new infrastructure. Excerpts below are the exact patterns to mirror.

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|----------------|---------------|
| `src/analysis/subcellular-location.ts` (NEW) | service | request-response + transform + cache | `src/analysis/enrichment.ts` (fetchProxy client) + `src/panels/uniprot-panel.ts` (`fetchUniProtData`, `parseAccession`) | role-match (compose two analogs) |
| `src/viewers/volcano.ts` (MODIFY) | viewer | transform | itself (`ensureNegLog10Column`, `ensureDirectionColumn`, `createVolcanoPlot`) | exact (in-place parameterize) |
| `src/analysis/enrichment.ts` (MODIFY) | service | request-response | itself (`gGOSt`, `runEnrichmentPipeline`, `showEnrichmentDialog`) | exact |
| `src/viewers/enrichment-viewers.ts` (MODIFY) | viewer | event-driven | itself (`createEnrichmentDotPlot`, `wireEnrichmentToVolcano`, `openEnrichmentVisualization`) | exact |
| `src/parsers/spectronaut-candidates-parser.ts` (MODIFY) | parser | transform (pure) | itself (`parseSpectronautCandidatesText` `getRawData` block) | exact |
| `src/analysis/differential-expression.ts` (MODIFY) | service | transform | itself (`showDEDialog` `pairs`/`comparisonInput`) | exact |
| `src/panels/uniprot-panel.ts` (MODIFY) | panel | request-response | itself (`fetchUniProtData`) — refactor to delegate to new shared module | exact |
| `src/utils/proteomics-types.ts` (MODIFY) | config | n/a | itself (`SEMTYPE` const map) | exact |
| `detectors.js` (MODIFY) | config | n/a | itself (`detectGeneSymbol` shape) | exact |
| `src/package.ts` (MODIFY) | controller | request-response | itself (`importSpectronautCandidates`, `showVolcanoPlot`) | exact |
| `src/tests/subcellular-location.ts` (NEW), `src/tests/volcano.ts` (NEW) | test | n/a | `src/tests/spectronaut-candidates-parser.ts` | role-match |

---

## Pattern Assignments

### `src/analysis/subcellular-location.ts` (NEW — service: fetch + verbatim classifier + cache)

This file composes three existing patterns. It is consumed by both `volcano.ts`
(R1/D-05 color column) and the refactored `uniprot-panel.ts` (folded cache todo).

**Pattern A — fetchProxy external HTTP (CORS-safe)**
**Analog:** `src/panels/uniprot-panel.ts` lines 77-91 and `src/analysis/enrichment.ts` lines 55-69.

```typescript
// src/panels/uniprot-panel.ts:77-91 — the exact try/warn/continue shape to mirror
async function fetchUniProtData(accession: string): Promise<UniProtEntry | null> {
  const url = `https://rest.uniprot.org/uniprotkb/${encodeURIComponent(accession)}.json` +
    '?fields=accession,protein_name,gene_names,organism_name,cc_function,go';
  try {
    const resp = await grok.dapi.fetchProxy(url);
    if (!resp.ok) {
      console.warn(`UniProt fetch for ${accession} returned status ${resp.status}`);
      return null;
    }
    return await resp.json() as UniProtEntry;
  } catch (e: any) {
    console.warn(`UniProt fetch for ${accession} failed:`, e?.message ?? e);
    return null;
  }
}
```

For the **stream** endpoint (D-01), keep the same `grok.dapi.fetchProxy` + `!resp.ok` warn-and-continue
discipline but build the chunked OR-query and parse TSV positionally (RESEARCH.md
"Verified UniProt stream response shape" — header is `Entry`, NOT `accession`; index `parts[0..4]`,
skip `lines[0]`). RESEARCH.md §"Code Examples / fetchProxy to UniProt stream" gives the literal shape.
`enrichment.ts:55-69` `fetchWithTimeout` is the analog if a timeout race is wanted (g:Profiler precedent).

**Pattern B — `parseAccession` for protein-group ID extraction (reuse, do not re-implement)**
**Analog:** `src/panels/uniprot-panel.ts` lines 40-71 (already exported, already imported by `enrichment.ts:4`).
Feed `parseAccession(rawProteinId)` per row to get the primary accession; the parser handles
`sp|…|…`, `;`-delimited groups (first wins), and `ups` suffix. CK-omics' "first UID in map" rule
maps onto this directly.

**Pattern C — persistent cache via platform KV**
**Analog (verified API):** `js-api/src/dapi.ts:735-760` `UserDataStorage`. RESEARCH.md §"Persistent
cache via platform KV":

```typescript
const STORE = 'proteomics-subcell-loc';
const cached = await grok.dapi.userDataStorage.get(STORE) ?? {}; // {acc: category}
// ... compute misses, fetch ...
await grok.dapi.userDataStorage.put(STORE, {...cached, ...fetched}); // write-through
```
Discretion: include a `__v` schema-version key (D-04 keyword map is a frozen contract → fixed version
is safe). Assumption A4: validate ~8329-entry map fits as one value in a Wave-0 spike; fall back to
prefix-sharding or IndexedDB if it bites (D-02 grants the discretion).

**Pattern D — verbatim classifier port (no I/O, pure)**
The 11-category keyword map (insertion order load-bearing for tie-breaks), exact hex palette, and
`parseSubcellularLocation(subcell, go)` algorithm are reproduced byte-for-byte in RESEARCH.md
§"CK-omics Algorithm Port: Subcellular Classification". Port verbatim — `\b`-bounded
case-insensitive substring scan over the **lowercased raw** string, min `match.index` across all
categories, subcell field fully before GO, `'Unknown'` default. **Do not strip `{ECO:…}` /
`SUBCELLULAR LOCATION:` first** (Pitfall 2).

---

### `src/viewers/volcano.ts` (MODIFY — viewer: parameterize metric + add color dimension)

**Analog: itself.** D-05/D-06 land here. The load-bearing requirement (D-06): "dots and colouring
never disagree" → one synchronized recompute function.

**Current significance-hardwire to fix (lines 13-29, parameterize the p-column):**
```typescript
export function ensureNegLog10Column(df: DG.DataFrame): string {
  const colName = 'negLog10padj';
  if (df.columns.contains(colName))
    return colName;                       // ← Pitfall 5: early-return defeats live toggle.
  const adjPCol = df.col('adj.p-value');  // ← hardwired metric — make this a parameter
  ...
  const col = df.columns.addNewFloat(colName);
  col.init((i) => {
    if (adjPCol.isNone(i)) return DG.FLOAT_NULL;
    const p = adjPCol.get(i) as number;
    return p > 0 ? -Math.log10(p) : UNDERFLOW_NEGLOG10;
  });
```
Change: take `metric: 'adj.p-value' | 'p-value'`, keep a **fixed generic column name**
(e.g. `negLog10p`), and **always re-`init` in place** (drop the early-return) so the Y-binding name
is stable but the values change on toggle (Pitfall 5; same in-place discipline as
`ensureDirectionColumn` at line 50).

**Direction column in-place update pattern to preserve (lines 48-62):**
```typescript
// Update existing column in place rather than remove + re-add, which would
// invalidate viewer/grid bindings to the column reference.
const col = df.col(colName) ?? df.columns.addNewString(colName);
col.init((i) => {
  if (fcCol.isNone(i) || adjPCol.isNone(i)) return 'not significant';
  const fc = fcCol.get(i) as number;
  const adjP = adjPCol.get(i) as number;
  if (adjP <= pThreshold && fc > fcThreshold) return 'up';
  if (adjP <= pThreshold && fc < -fcThreshold) return 'down';
  return 'not significant';
});
col.meta.colors.setCategorical(
  {'up': 0xFFCC0000, 'down': 0xFF0000CC, 'not significant': 0xFFAAAAAA},
);
```
Parameterize the p-column the same way; **do NOT remove/re-add on toggle** (the comment is the
contract). Mirror this exact `meta.colors.setCategorical(ARGB)` idiom for the 11-color location map
(D-05) — ARGB ints from RESEARCH.md §"Categorical color map on the location column".

**Threshold-line stale-cleanup pattern to reuse on every toggle (lines 94-118):**
```typescript
const hLine = -Math.log10(pThreshold);
const yFormulaPrefix = `\${${yColName}}`;
const fcFormulaPrefix = `\${log2FC}`;
df.meta.formulaLines.items = df.meta.formulaLines.items.filter((line) => {
  const f = line.formula ?? '';
  return !(typeof f === 'string' && (f.startsWith(yFormulaPrefix) || f.startsWith(fcFormulaPrefix)));
});
df.meta.formulaLines.addLine({formula: `${yFormulaPrefix} = ${hLine}`, color: '#888888', width: 1, visible: true});
```
On a Q↔P toggle, re-run this filter+addLine block so lines replace, not stack. RESEARCH.md
§"Volcano Q/P Live Toggle" gives the target `recomputeVolcano(df, sp, metric, colorDim, fcThr, pThr)`
signature. Guard `'p-value'` when `df.col('p-value') == null` (Candidates only adds it `if (pValCol)` —
see candidates parser line 124).

---

### `src/analysis/enrichment.ts` (MODIFY — service: up/down split + WP source)

**Analog: itself.** D-10/D-11 (caller change only — RESEARCH.md is explicit: keep `gGOSt`'s
request body unchanged, do **not** adopt CK-omics highlight-filtering, that is Phase 14).

**Existing single-query path to split by sign (lines 295-321):**
```typescript
const significantGenes: Set<string> = new Set();
const backgroundGenes: Set<string> = new Set();
const fcRaw = cols.log2fc.getRawData() as Float32Array | Float64Array;
const pRaw = cols.pValue.getRawData() as Float32Array | Float64Array;
for (const [row, gene] of geneForRow) {
  backgroundGenes.add(gene);                       // ← background stays ALL detected (D-10)
  const fc = fcRaw[row];
  const adjP = pRaw[row];
  if (fc !== DG.FLOAT_NULL && adjP !== DG.FLOAT_NULL &&
      adjP <= pThreshold && Math.abs(fc) >= fcThreshold)
    significantGenes.add(gene);
}
...
const results = await gGOSt(queryArray, bgArray, organismCode, sources, pThreshold);
const enrichmentDf = buildEnrichmentDf(results, queryArray, pThreshold);
```
D-10 change: split significant genes into `upGenes` (fc > 0) / `downGenes` (fc < 0), call `gGOSt`
twice with the **same `bgArray`**, then merge: build per-direction arrays, concat, add a `Direction`
categorical string column ('Up'/'Down') to the `buildEnrichmentDf` output shape (lines 134-185 —
mirror the `df.columns.addNewString(...).init(...)` bulk pattern).

**WP source — one-line additive change (D-12), analog `showEnrichmentDialog` lines 348-396):**
```typescript
// existing source checkboxes (lines 349-353)
const goBpInput = ui.input.bool('GO: Biological Process', {value: true});
...
const reactomeInput = ui.input.bool('Reactome Pathways', {value: true});
// existing OK-handler source build (lines 391-396)
if (goBpInput.value) selectedSources.push('GO:BP');
...
if (reactomeInput.value) selectedSources.push('REAC');
```
Add `const wpInput = ui.input.bool('WikiPathways', {value: true});` + `.add(wpInput)` and
`if (wpInput.value) selectedSources.push('WP');` (literal `'WP'`, default-on). `gGOSt` passes
`sources` through unchanged (lines 90-117).

---

### `src/viewers/enrichment-viewers.ts` (MODIFY — viewer: Direction-aware split charts)

**Analog: itself.** D-11 extends the Phase-9 dot/bar + cross-link, does not rebuild it.

**Cross-link to preserve unchanged (lines 92-134) — the Phase-9 `onCurrentRowChanged` pattern:**
```typescript
return enrichDf.onCurrentRowChanged.subscribe(() => {
  const rowIdx = enrichDf.currentRowIdx;
  if (rowIdx < 0) return;
  const intersectionCol = enrichDf.col('Intersection');
  if (!intersectionCol) return;
  const memberGenesStr = intersectionCol.get(rowIdx) as string;
  ...
  proteinDf.selection.setAll(false, false);
  for (const gene of memberGenes) {
    const rows = geneToRows.get(gene);
    if (rows) for (const row of rows) proteinDf.selection.set(row, true, false);
  }
  proteinDf.selection.fireChanged();
});
```
D-10 says this is **unchanged** — the merged DataFrame keeps the `Intersection` column; the
`Direction` column is an extra categorical. Do not break the `proteinDf.selection` wiring
(Anti-pattern: "Mutating the source DataFrame in the enrichment cross-link").

**Dock idiom for the side-by-side Up/Down split (lines 158-170) — the established package idiom:**
```typescript
const tv = existing ?? grok.shell.addTableView(enrichDf);
const dotPlot = createEnrichmentDotPlot(topDf);
const barChart = createEnrichmentBarChart(topDf);
const dotNode = tv.dockManager.dock(dotPlot, DG.DOCK_TYPE.RIGHT, null, 'Dot Plot', 0.5);
tv.dockManager.dock(barChart, DG.DOCK_TYPE.DOWN, dotNode, 'Bar Chart', 0.5);
```
D-11: produce split Up/Down dot+bar (filter `topDf` by `Direction`, or facet). The
`tv.dockManager.dock(viewer, DG.DOCK_TYPE.*, anchorNode, title, ratio)` call is the **exact docking
idiom R4 also reuses** for the Comparison Filter viewer in `package.ts`. `createEnrichmentDotPlot`
(lines 58-69) / `createEnrichmentBarChart` (lines 74-85) are the factories to clone per direction;
`createTopNEnrichmentDf` (lines 14-52) is the clone-for-isolation pattern (`enrichDf.clone(mask)`).
Subscription cleanup pattern: module-level `activeSubscriptions` array, `unsubscribe()` on re-open
(lines 7-8, 145-147, 170).

---

### `src/parsers/spectronaut-candidates-parser.ts` (MODIFY — parser: R3/D-08 sign-flip, stay PURE)

**Analog: itself.** Extend the existing bulk `getRawData()` block; **no `grok.shell`** (Pitfall 3 —
breaks the entire `SpectronautCandidates` test category).

**Existing bulk-read block to extend (lines 138-146):**
```typescript
const sigCol = df.columns.addNewBool('significant');
const fcRaw = df.col('log2FC')!.getRawData() as Float32Array | Float64Array;
const adjPRaw = df.col('adj.p-value')!.getRawData() as Float32Array | Float64Array;
sigCol.init((i) => {
  const fc = fcRaw[i];
  const adjP = adjPRaw[i];
  if (fc === DG.FLOAT_NULL || adjP === DG.FLOAT_NULL) return false;
  return Math.abs(fc) >= DEFAULT_FC_THRESHOLD && adjP <= DEFAULT_P_THRESHOLD;
});
```
D-08: before this block, parse each row's `Comparison (group1/group2)` string (`A / B`), build a
per-row sign-multiplier `Float32Array`, write a flipped `log2FC` array + (conditionally) swapped
`AVG Group Quantity Numerator`/`Denominator` arrays, then `col.init((i) => arr[i])` from them
(bulk pattern, memory `feedback_dg_column_bulk_init`; FLOAT_NULL handling exactly as the `sigCol`
block above — Pitfall 6). **Per-row decision only** (Pitfall 4 — unconditional inversion
re-introduces the mirror defect; A2 is HIGH-risk: verify against a real BP DMD/WT Candidates file
in Wave-0). Guard the quantity-column swap with `if (col present)` (A6, mirrors CK-omics).

**Column-detection / canonical-rename idiom to reuse for the AVG Group Quantity columns
(lines 34-42, 61-69):**
```typescript
function findCol(df: DG.DataFrame, names: readonly string[]): DG.Column | null {
  for (const name of names) { const col = df.col(name); if (col) return col; }
  return null;
}
function renameToCanonical(df: DG.DataFrame, src: DG.Column, canonical: string): void {
  if (src.name === canonical) return;
  const existing = df.col(canonical);
  if (existing && existing !== src) return;
  src.name = canonical;
}
```
Add `AVG_GROUP_QTY_*` name-variant arrays at the top alongside `LOG2FC_COLUMNS` (lines 12-22) and
detect via `findCol`. Verbatim flip mechanics in RESEARCH.md §"CK-omics Algorithm Port:
create_subset_data Flip". RESEARCH.md note: the candidates parser sets
`proteomics.de_complete='true'` (line 149) — D-09's report-DE fix is the *report* path, NOT this one
(Anti-pattern: "Re-running DE on a Candidates frame").

---

### `src/analysis/differential-expression.ts` (MODIFY — service: D-09 declared-contrast default)

**Analog: itself.** Direction-only change (memory: "does not move hit-list metrics").

**The defect — default Comparison = alphabetical `pairs[0]` (lines 288-307):**
```typescript
const pairs = [`${g2.name} vs ${g1.name}`, `${g1.name} vs ${g2.name}`];
const comparisonInput = ui.input.choice('Comparison', {
  value: pairs[0],                              // ← alphabetical default = the parked mirror defect
  items: pairs, nullable: false,
});
...
const isReversed = comparisonInput.value === pairs[1];
const numerator = isReversed ? g1.name : g2.name;
```
**OK-handler that consumes the choice (lines 355-361) — keep this index-based mapping intact:**
```typescript
const reversed = comparisonInput.value === pairs[1];
const numeratorCols = reversed ? g1.columns : g2.columns;
const denominatorCols = reversed ? g2.columns : g1.columns;
const numeratorName = reversed ? g1.name : g2.name;
const denominatorName = reversed ? g2.name : g1.name;
```
D-09: set the **default** `value:` to the declared/intended orientation instead of `pairs[0]`; the
dropdown stays as an override. **A3 (MEDIUM):** the origin of report group order
(`spectronaut-parser.ts` vs `analysis/experiment-setup.ts` `getGroups`/`setGroups` — interface at
`experiment-setup.ts:8-23`) is **not yet traced** and was not in the read set — plan a Wave-0 trace
task before the fix. `getGroups(df)` (line 269) is the read point; group order is JSON in
`proteomics.groups` via `setGroups` (always go through these per package convention, never parse
the tag).

---

### `src/panels/uniprot-panel.ts` (MODIFY — refactor to delegate to shared module)

**Analog: itself.** D-02 fold: `fetchUniProtData` (lines 77-91) becomes a thin caller of the new
`subcellular-location.ts` shared fetch+cache. `parseAccession` (lines 40-71) **stays here and stays
exported** (already imported by `enrichment.ts:4` and reused by the new module — do not move it).
The widget render path (lines 120-206, `ui.wait(async …)`) is unchanged.

---

### `src/utils/proteomics-types.ts` + `detectors.js` (MODIFY — config: new SEMTYPE, mirrored)

**Analog: itself.** Package contract: a new `SEMTYPE.*` requires a mirrored `detectors.js` entry.

`proteomics-types.ts` (entire file — add one line):
```typescript
export const SEMTYPE = {
  PROTEIN_ID: 'Proteomics-ProteinId',
  GENE_SYMBOL: 'Proteomics-GeneSymbol',
  LOG2FC: 'Proteomics-Log2FC',
  P_VALUE: 'Proteomics-PValue',
  INTENSITY: 'Proteomics-Intensity',
  // add: SUBCELLULAR_LOCATION: 'Proteomics-SubcellularLocation',
} as const;
```
`detectors.js` — mirror the `detectGeneSymbol` shape (lines 22-34; plain JS, `//meta.role:
semTypeDetector`, string literal must equal the new `SEMTYPE` value — these are NOT shared at
runtime, the literal is duplicated by contract):
```javascript
//meta.role: semTypeDetector
//input: column col
//output: string semType
detectGeneSymbol(col) {
  if (col.type !== DG.TYPE.STRING) return null;
  const name = col.name.toLowerCase();
  if (name === 'gene names' || name === 'gene name' || name === 'gene symbol' || name === 'gene') {
    col.semType = 'Proteomics-GeneSymbol';
    return col.semType;
  }
  return null;
}
```
The location column is assigned its semType in code (`volcano.ts` color path:
`col.semType = SEMTYPE.SUBCELLULAR_LOCATION`); the detector exists so vendor tables that already
carry a "Subcellular location" column get typed. Name-hint: `'subcellular location'` /
`'subcellular'`.

---

### `src/package.ts` (MODIFY — controller: Filter-viewer docking + volcano toggle menu)

**Analog: itself.** Shell orchestration only — parser stays pure (Pattern 1).

**Import handler to extend for R4/D-07 (lines 156-172):**
```typescript
static async importSpectronautCandidates(): Promise<void> {
  DG.Utils.openFile({
    accept: '.tsv,.txt,.csv',
    open: async (file: File) => {
      try {
        const text = await file.text();
        const df = parseSpectronautCandidatesText(text);
        df.name = file.name.replace(/\.[^.]+$/, '');
        grok.shell.addTableView(df);          // ← after this: count distinct Comparison,
        grok.shell.info(`Imported ${df.rowCount} candidates from Spectronaut`);
      } catch (e: any) {
        grok.shell.error(`Failed to import Spectronaut Candidates file: ${e?.message ?? e}`);
      }
    },
  });
}
```
D-07: after `addTableView(df)`, count distinct non-null Comparison values; if > 1 dock a native
Filter viewer using the **enrichment-viewers dock idiom**:
`tv.dockManager.dock(filterViewer, DG.DOCK_TYPE.RIGHT, null, 'Comparison', 0.3)` scoped to the
Comparison column via the Filters viewer `columnNames` look property. Single-comparison files: skip.

**Volcano menu handler to extend for D-05/D-06 (lines 236-247):**
```typescript
@grok.decorators.func({'top-menu': 'Proteomics | Visualize | Volcano Plot...'})
static async showVolcanoPlot(): Promise<void> {
  const tv = grok.shell.tv;
  const df = tv?.dataFrame;
  if (!tv || !df) { grok.shell.warning('No table open'); return; }
  if (!requireDifferentialExpression(df, 'Run Differential Expression first ...')) return;
  const groups = getGroups(df);
  const title = groups ? `Volcano: ${groups.group2.name} vs ${groups.group1.name}` : 'Volcano';
  const sp = createVolcanoPlot(df, {title});
  tv.addViewer(sp);
}
```
RESEARCH.md recommends approach (A): keep `createVolcanoPlot` a `ScatterPlotViewer` factory; add a
`Proteomics | Visualize | Volcano Options…` dialog (`ui.input.choice` metric + color-by) or a
`sp.onContextMenu` item that calls one `recomputeVolcano(df, sp, …)`. Defaults: metric
`adj.p-value`, color `significance`. Menu paths that open a dialog end with `...` (package
naming convention).

---

## Shared Patterns

### External HTTP — `grok.dapi.fetchProxy` (NEVER raw `fetch`)
**Source:** `src/panels/uniprot-panel.ts:81`, `src/analysis/enrichment.ts:61`
**Apply to:** `subcellular-location.ts` (UniProt stream + D-03 reviewed-by-gene fallback)
CLAUDE.md hard rule; CONCERNS.md documents the raw-fetch CORS defect this phase fixes. The
`enrichment.ts` `fetchWithTimeout` (lines 55-69) is the timeout-race precedent if needed.

### Bulk `Column.init` from a typed array — never per-row `col.set`
**Source:** `src/parsers/spectronaut-candidates-parser.ts:138-146`, `src/viewers/volcano.ts:22-27`,
`src/analysis/enrichment.ts:279-285` (`geneArr` → `geneCol.init`)
**Apply to:** R1 location string column (8329 rows), R3 log2FC sign-flip + AVG-quantity swap,
D-10 `Direction` column.
Memory `feedback_dg_column_bulk_init` (255× measured). Landmine: `col.init(() => null)` on a numeric
column leaves the FLOAT_NULL sentinel read back as finite positive (memory
`feedback_dg_column_init_null_sentinel`) — for numeric flips write a real `Float32Array` and
`isNone()`/explicit FLOAT_NULL exactly as the `sigCol` block does (Pitfall 6). Moot for the location
*string* column.

### In-place column update to preserve viewer/grid bindings
**Source:** `src/viewers/volcano.ts:48-50` (the comment IS the contract)
**Apply to:** the D-06 metric toggle (`ensureNegLog10Column`/`ensureDirectionColumn` must re-`init`
the existing column under a fixed name, never remove+re-add).

### Categorical color via `col.meta.colors.setCategorical(ARGB)`
**Source:** `src/viewers/volcano.ts:60-62`
**Apply to:** D-05 11-color subcellular-location map (ARGB ints in RESEARCH.md §"Categorical color
map on the location column"; convert `#RRGGBB` → `0xFF000000 | parseInt(hex.slice(1),16)`).

### `tv.dockManager.dock(viewer, DG.DOCK_TYPE.*, anchorNode, title, ratio)`
**Source:** `src/viewers/enrichment-viewers.ts:164-165`
**Apply to:** R4 Comparison Filter viewer (`package.ts`), D-11 side-by-side Up/Down charts.

### Pure parser, shell orchestration in `package.ts`
**Source:** `src/parsers/spectronaut-candidates-parser.ts` (zero `grok.shell`) +
`src/package.ts:156-172` (`addTableView` / dock here)
**Apply to:** R3 sign-flip = parser; R4 Filter dock = `package.ts`. Test-enforced (Phase 12):
`src/tests/spectronaut-candidates-parser.ts:26` calls `parseSpectronautCandidatesText(text)` with
no shell — any `grok.shell` in the parser breaks the category (Pitfall 3).

### Test file structure (`@datagrok-libraries/test`)
**Source:** `src/tests/spectronaut-candidates-parser.ts:1-2,20-21`
```typescript
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {parseSpectronautCandidatesText} from '../parsers/spectronaut-candidates-parser';
category('SpectronautCandidates', () => {
  test('parses standard rows', async () => { ... expect(df.rowCount, 2); });
});
```
**Apply to:** NEW `src/tests/subcellular-location.ts` (pure: classifier + palette + TSV positional
parse), NEW `src/tests/volcano.ts` (parameterized helpers + in-place recompute). Extend existing
`src/tests/{enrichment,spectronaut-candidates-parser,spectronaut-candidates-e2e,analysis}.ts`.
Register through `src/package-test.ts`. Memory `feedback_grok_test_skipbuild_stale`: rebuild after
adding test files (`--skip-build` reuses stale bundle → new tests return null/0).

---

## No Analog Found

None. Every file in scope has an exact in-repo analog (itself or a sibling). The only genuinely
new module, `src/analysis/subcellular-location.ts`, is a composition of three verified existing
patterns (fetchProxy client + `parseAccession` + `userDataStorage` KV) plus a verbatim Python port
whose contract is reproduced in full in RESEARCH.md — the planner copies, does not invent.

---

## Metadata

**Analog search scope:** `packages/Proteomics/src/{viewers,analysis,parsers,panels,utils,tests}/`,
`packages/Proteomics/src/package.ts`, `packages/Proteomics/detectors.js`
**Files scanned:** 13 (10 integration-surface files read in full + experiment-setup/shared-utils/
column-detection signatures + test directory listing)
**Pattern extraction date:** 2026-05-16
**Key risks flagged for the planner (from RESEARCH.md Assumptions Log):**
A2 (HIGH — verify real BP DMD/WT Candidates sign convention before locking R3),
A3 (MEDIUM — trace report-DE group-order origin in Wave-0),
A4 (LOW-MED — validate userDataStorage value size in Wave-0).
