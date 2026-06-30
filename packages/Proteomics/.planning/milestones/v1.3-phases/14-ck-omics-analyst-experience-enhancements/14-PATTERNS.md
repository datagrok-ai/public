# Phase 14: CK-omics Analyst-Experience Enhancements — Pattern Map

**Mapped:** 2026-05-28
**Files analyzed:** 15 (4 NEW, 11 MODIFIED)
**Analogs found:** 15 / 15 (every change has an in-tree precedent)

> Phase 14 is a layered enhancement on the Phase 13 foundation. Every new
> capability extends an existing file, or — for the four NEW files — copies
> an exact in-tree precedent. There is no greenfield path: the planner's job
> is wiring, not invention. The pattern excerpts below are LOAD-BEARING:
> deviating from them re-introduces the bugs that Phase 13's UAT closed.

---

## File Classification

| # | New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|---|-------------------|------|-----------|----------------|---------------|
| 1 | `src/utils/gene-label-resolver.ts` **(NEW)** | util / external-REST resolver | request-response + cross-session cache + bounded-concurrency worker pool | `src/analysis/subcellular-location.ts` (lines 165-393) | EXACT — copy cache / chunk / concurrency / fallback shape verbatim |
| 2 | `src/viewers/group-mean-correlation.ts` **(NEW)** | viewer factory | derived-column transform + scatter render | `src/viewers/volcano.ts:createVolcanoPlot` (lines 193-222) for factory + `src/viewers/qc-computations.ts:ensureFreshFloat/computeMA` (lines 28-55) for derived means + `src/viewers/pca-plot.ts:createPcaPlot` (lines 99-169) for `setOptions({title})` annotation | EXACT for factory shape, EXACT for ensureFreshFloat pattern |
| 3 | `src/tests/gene-label-resolver.ts` **(NEW)** | test (unit) | fixture + mocked `grok.dapi.fetchProxy` + assertions | `src/tests/spectronaut-parser.ts` (lines 1-80) for category/test/expect framing + cache-mock pattern | role-match (no existing Ensembl mock yet) |
| 4 | `src/tests/group-mean-correlation.ts` **(NEW)** | test (unit) | known-fixture correctness | `src/tests/parsers.ts` (lines 1-60) for compact `category/test/expect` shape | role-match |
| 5 | `src/viewers/volcano.ts` **(MODIFY)** | viewer factory | extend in place — direction strings, color migration, counter overlay, default top-N labels, axis-label overlay, dialog-state surface | self (lines 62-86 `ensureDirectionColumn`, lines 193-222 `createVolcanoPlot`, lines 231-250 `recomputeVolcano`) | self-extension |
| 6 | `src/package.ts` **(MODIFY)** | entry / shell-orchestration | extend `dockComparisonFilterIfMultiContrast` (lines 117-131) + extend `volcanoOptions` handler (lines 274-355) + register new `Group-Mean Correlation` menu item | self | self-extension |
| 7 | `src/panels/uniprot-panel.ts` **(MODIFY)** | context panel | append per-group SVG bar chart inside `renderUniProtWidget` (lines 109-176) after the GO Terms section | self (lines 152-173 GO Terms section structure) + `viewers/qc-computations.ts:groupMean` (lines 15-25) for the mean math | self-extension |
| 8 | `src/analysis/enrichment.ts` **(MODIFY)** | analysis dialog | extend `showEnrichmentDialog` (lines 371-467) with one `ui.input.bool` + add `applySmartPathwayFilter` transform between `gGOSt` (line 344) and `buildEnrichmentDf` (line 345) | self | self-extension |
| 9 | `src/parsers/maxquant-parser.ts` **(MODIFY)** | parser | add `await resolveGeneLabels(df)` post-process call before `return df` (line 158) | self + sibling parsers | self-extension (one-line addition × 5 parsers) |
| 10 | `src/parsers/spectronaut-parser.ts` **(MODIFY)** | parser | add `await resolveGeneLabels(df)` before `return df` (line 185) | self | self-extension |
| 11 | `src/parsers/spectronaut-candidates-parser.ts` **(MODIFY)** | parser | add `await resolveGeneLabels(df)` before `return df` (line 267) | self | self-extension |
| 12 | `src/parsers/fragpipe-parser.ts` **(MODIFY)** | parser | add `await resolveGeneLabels(df)` before `return df` (line 196) | self | self-extension |
| 13 | `src/parsers/generic-parser.ts` **(MODIFY)** | parser | call `await resolveGeneLabels(df)` inside the dialog onOK before `addTableView` | self (dialog onOK path) | self-extension |
| 14 | `src/utils/proteomics-types.ts` **(MODIFY)** | constants — single source of truth for SEMTYPE strings | add 4 new SEMTYPE constants: `DISPLAY_NAME`, `SOURCE_ID`, `NUMERATOR_MEAN`, `DENOMINATOR_MEAN` | self (lines 1-9 existing SEMTYPE block) | self-extension |
| 15 | `detectors.js` **(MODIFY — plain JS, not TS)** | semantic-type detector | add 4 mirror detectors matching the new SEMTYPE constants | self (lines 25-49 `detectGeneSymbol` / `detectSubcellularLocation`) | self-extension |
| — | `src/tests/volcano.ts` **(MODIFY)** | test (unit) | update assertions: new direction strings, color hex, top-N selection, counter math | self | self-extension |
| — | `src/tests/enrichment.ts` **(MODIFY)** | test (unit) | add smart-filter cap+drop assertions | self | self-extension |
| — | `src/tests/uniprot-panel.ts` **(MAYBE NEW)** | test (unit) | per-group bar chart math; existence depends on whether a panel test file already exists — planner verifies | `src/tests/parsers.ts` shape | role-match |
| — | `src/package-test.ts` **(MODIFY)** | test entry | one `import './tests/gene-label-resolver';` + one `import './tests/group-mean-correlation';` line per new test file | self | self-extension |

---

## Pattern Assignments

### 1. `src/utils/gene-label-resolver.ts` (NEW — util, external-REST resolver)

**Analog:** `src/analysis/subcellular-location.ts` (the ONLY in-tree reference for the
exact cache + chunk + worker-pool + fallback shape Phase 14 D-09 requires).

#### Imports pattern (subcellular-location.ts:1-2)

```typescript
import * as grok from 'datagrok-api/grok';
// NO datagrok-api/dg import here — this is a pure data layer.
// All df-mutation belongs in the CALLERS (the parsers), not in the resolver.
```

For Phase 14 the resolver returns a `Map<rawId, {displayName, sourceId, marker}>`;
the parsers do the `df.columns.addNewString(...).init(...)` writes. This mirrors
how `subcellular-location.ts` returns a `Map<acc, category>` and `volcano.ts`
(`ensureLocationColumn`, lines 120-168) does the column write.

#### Cache shape (subcellular-location.ts:165-253)

Copy verbatim — only the constant names and SCHEMA_V change:

```typescript
// VERBATIM PATTERN — subcellular-location.ts:165-174 with renamed constants
export const STORE_GENE_LABELS = 'proteomics-gene-labels';
const SCHEMA_V_GENE_LABELS = '14-r1-1';  // bump only if the resolver contract changes
const SCHEMA_KEY = '__schema_v';
const ENS_CHUNK = 1000;  // Ensembl /lookup/id hard cap (Pitfall 3 in 14-RESEARCH.md)

async function loadCache(): Promise<Record<string, GeneLabelEntry>> {
  try {
    const raw = (await grok.dapi.userDataStorage.get(STORE_GENE_LABELS)) ?? {};
    if (raw[SCHEMA_KEY] !== SCHEMA_V_GENE_LABELS)
      return {[SCHEMA_KEY]: SCHEMA_V_GENE_LABELS}; // stale → discard
    return raw as Record<string, GeneLabelEntry>;
  } catch {
    return {[SCHEMA_KEY]: SCHEMA_V_GENE_LABELS};
  }
}
```

Key invariants the analog enforces — copy them all:
- `userDataStorage` is fresh-instance-per-access (project memory `reference_dapi_fresh_instance_patch`). Never cache the storage reference.
- The `__schema_v` key is OUTSIDE the cached values — schema mismatch means discard everything except a fresh schema entry.
- Failed reads return an empty cache + the schema marker — they do not throw.

#### Worker pool — REUSE not duplicate (subcellular-location.ts:219-237)

```typescript
// Already exported for tests. Import it directly:
import {runWithConcurrency} from '../analysis/subcellular-location';

// Use the SAME concurrency cap. Ensembl's documented limit is 15 req/sec
// (55k/hour). UniProt's effective working cap was set to 6 — same cap here.
const FETCH_CONCURRENCY = 6;
```

#### Chunked POST pattern (subcellular-location.ts:301-331)

The analog does GET; Phase 14 needs POST with a JSON body. The wrapper shape is identical
— chunk → `runWithConcurrency` → `grok.dapi.fetchProxy` → warn-and-continue on failure.

```typescript
// PATTERN — mirror subcellular-location.ts:306-331, swap GET for POST.
const idChunks = chunk(misses, ENS_CHUNK);
let done = 0;
await runWithConcurrency(idChunks, FETCH_CONCURRENCY, async (group) => {
  try {
    const resp = await grok.dapi.fetchProxy('https://rest.ensembl.org/lookup/id', {
      method: 'POST',
      headers: {'Content-Type': 'application/json', 'Accept': 'application/json'},
      body: JSON.stringify({ids: group}),
    });
    if (resp.status === 429) {
      const retryAfter = parseFloat(resp.headers.get('Retry-After') ?? '1');
      await new Promise((r) => setTimeout(r, retryAfter * 1000));
      // single retry — exactly one extra request per chunk on rate-limit
    }
    if (!resp.ok) {
      console.warn(`Ensembl /lookup/id returned ${resp.status}; continuing`);
      return;
    }
    const json = await resp.json();
    // Unknown IDs come back as null values keyed by ID (Assumption A2 in RESEARCH).
    for (const [id, entry] of Object.entries(json)) {
      if (entry && typeof entry === 'object') fetched[id] = projectEntry(entry);
    }
  } catch (e: any) {
    console.warn(`Ensembl batch failed: ${e?.message ?? e}; continuing`);
  }
  done++;
  progress?.(done, idChunks.length, 'fetch-ensembl');
});
```

#### Timer-driven flush (subcellular-location.ts:289-298)

Mirror EXACTLY. The previous draft of `subcellular-location.ts` folded the flush into the
worker body and hit a double-flush race; the documented fix is a single `setInterval` that
the orchestrator owns, NOT the workers. Same applies here.

#### Fallback chain (CK-omics contract, R1 / D-09)

Three-level fallback per locked client contract:
1. `external_name` (Ensembl) → use as Display Name
2. `description` (Ensembl) → run `extractReadableDescription` (CK-omics `extract_readable_description` port)
3. raw ID → keep verbatim, mark with `†`

Acceptance gate (CKomics_tool2.py:977-988): reject candidate name if it equals the raw ID
or starts with any predicted-pattern prefix. Fall through to description if rejected.

#### What the resolver MUST write (R1 / D-08)

For every row, populate two columns the parsers will create:
- `Display Name` (string, semType=`SEMTYPE.DISPLAY_NAME`) — the canonical label everywhere
- `Source ID` (string, semType=`SEMTYPE.SOURCE_ID`) — raw ID for hover/disambiguation

**Pitfall 9 from RESEARCH:** the resolver MUST always create `Display Name`, even if no
predicted IDs are detected — populate with the unmodified gene name as a no-op. Otherwise
the volcano label binding silently breaks for files without predicted IDs.

---

### 2. `src/viewers/group-mean-correlation.ts` (NEW — viewer factory)

**Analog A — factory shape:** `src/viewers/volcano.ts:createVolcanoPlot` (lines 193-222).
**Analog B — derived columns:** `src/viewers/qc-computations.ts:ensureFreshFloat` + `groupMean` + `computeMA` (lines 15-55).
**Analog C — title annotation:** `src/viewers/pca-plot.ts:createPcaPlot` (lines 99-169 — uses `sp.setOptions({title})` at line 166).

#### Imports pattern (volcano.ts:1-9)

```typescript
import * as DG from 'datagrok-api/dg';
import {findColumn} from '../utils/column-detection';
import {SEMTYPE} from '../utils/proteomics-types';
import {getGroups, GroupAssignment} from '../analysis/experiment-setup';
import {DIRECTION_COL} from './volcano'; // direction column the volcano populates
```

#### Factory shape (volcano.ts:193-222)

```typescript
// PATTERN — mirror volcano.ts:193-222 EXACTLY (signature + return shape + setOptions title)
export function createGroupMeanCorrelation(
  df: DG.DataFrame,
  options?: {title?: string; topNLabels?: number},
): DG.ScatterPlotViewer {
  const groups = getGroups(df);
  if (!groups) throw new Error('Annotate experimental groups first');

  // Derived columns — re-run-safe via ensureFreshFloat (qc-computations.ts:28-32).
  computeGroupMeans(df, groups);  // adds 'Numerator Mean' / 'Denominator Mean' columns

  const sp = df.plot.scatter({
    x: 'Numerator Mean',
    y: 'Denominator Mean',
    color: DIRECTION_COL,  // reuses the magenta/cyan/gray map set by volcano's ensureDirectionColumn
  });

  const labelCol = findColumn(df, SEMTYPE.DISPLAY_NAME, ['display name', 'gene name']);
  if (labelCol) {
    sp.props.labelColumnNames = [labelCol.name];
    sp.props.showLabelsFor = 'MouseOverRow';
  }

  // y = x diagonal — same idiom as volcano.ts:183-188
  df.meta.formulaLines.addLine({
    formula: '${Numerator Mean} = ${Denominator Mean}',
    color: '#888888', width: 1, visible: true,
  });
  sp.props.showViewerFormulaLines = true;

  const {r, rho} = computePearsonSpearman(df, 'Numerator Mean', 'Denominator Mean');
  const baseTitle = options?.title ?? 'Group-Mean Correlation';
  sp.setOptions({title: `${baseTitle} — r=${r.toFixed(2)} (Pearson), ρ=${rho.toFixed(2)} (Spearman)`});

  return sp;
}
```

#### ensureFreshFloat pattern (qc-computations.ts:28-32)

```typescript
// VERBATIM — copy this helper into the new file (or import; it's currently file-private).
function ensureFreshFloat(df: DG.DataFrame, name: string): DG.Column {
  if (df.columns.contains(name))
    df.columns.remove(name);
  return df.columns.addNewFloat(name);
}
```

#### groupMean math (qc-computations.ts:15-25 + computeMA:35-55)

```typescript
// PATTERN — qc-computations.ts:35-55 shows the bulk-init writeback after per-row mean.
function computeGroupMeans(df: DG.DataFrame, groups: GroupAssignment): void {
  const numCols = groups.group1.columns.map((n) => df.col(n)!).filter(Boolean);
  const denCols = groups.group2.columns.map((n) => df.col(n)!).filter(Boolean);

  const n = df.rowCount;
  const numArr = new Float32Array(n);
  const denArr = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    const m1 = groupMean(numCols, i);
    const m2 = groupMean(denCols, i);
    numArr[i] = isNaN(m1) ? DG.FLOAT_NULL : m1;
    denArr[i] = isNaN(m2) ? DG.FLOAT_NULL : m2;
  }
  const numC = ensureFreshFloat(df, 'Numerator Mean');
  const denC = ensureFreshFloat(df, 'Denominator Mean');
  numC.init((i) => numArr[i]);
  denC.init((i) => denArr[i]);
  numC.semType = SEMTYPE.NUMERATOR_MEAN;
  denC.semType = SEMTYPE.DENOMINATOR_MEAN;
}
```

#### Inline Pearson / Spearman (RESEARCH §Pattern 9; no in-tree analog)

`@datagrok-libraries/statistics` only exports `kendallsTau` (RESEARCH source citation
`libraries/statistics/src/correlation-coefficient.ts:68`). Inline both formulas inside the
factory file — match the `getRawData()` + skip-null pattern used in volcano.ts:241 and
qc-computations.ts:68-69:

```typescript
function computePearsonSpearman(df: DG.DataFrame, xName: string, yName: string)
    : {r: number; rho: number} {
  const xRaw = df.col(xName)!.getRawData() as Float32Array | Float64Array;
  const yRaw = df.col(yName)!.getRawData() as Float32Array | Float64Array;
  const xs: number[] = []; const ys: number[] = [];
  for (let i = 0; i < df.rowCount; i++) {
    if (!df.filter.get(i)) continue;
    if (xRaw[i] === DG.FLOAT_NULL || yRaw[i] === DG.FLOAT_NULL) continue;
    if (isNaN(xRaw[i]) || isNaN(yRaw[i])) continue;
    xs.push(xRaw[i]); ys.push(yRaw[i]);
  }
  return {r: pearson(xs, ys), rho: pearson(rank(xs), rank(ys))};
}
```

---

### 3. `src/viewers/volcano.ts` (MODIFY — extend in place, NEVER fork)

**Analog:** self. Phase 13 D-05 LOCKED the "one volcano with property-toggle" model.

#### D-04 color + legend migration (extend `ensureDirectionColumn`, lines 62-86)

The full migration excerpt is in RESEARCH.md §"Code Examples > Verifying the volcano color-string migration":

```typescript
// REPLACES lines 62-86. Read group names from getGroups(df); emit "Enriched in <g1>" /
// "Enriched in <g2>" / "Not significant" strings; set ARGB ints magenta/cyan/gray.
export const DIRECTION_COLORS_BASE = {
  enrichedG1: 0xFFFF00FF,  // magenta
  enrichedG2: 0xFF00FFFF,  // cyan
  notSig:     0xFFAAAAAA,  // gray (matches existing volcano.ts:84 'not significant')
} as const;

// ... see RESEARCH.md §"Code Examples > Verifying the volcano color-string migration" for the full body.
```

**Pitfall 5 from RESEARCH:** any test in `src/tests/volcano.ts` asserting on literal
`'up'` / `'down'` / `'not significant'` will go red. Update tests as part of the SAME
task that flips the strings — do not let red tests linger across task boundaries.

#### D-03 default top-N labels (extend `createVolcanoPlot`, lines 209-213)

Current code binds `Gene name` only as a hover label:

```typescript
// CURRENT (volcano.ts:209-213)
const geneCol = findColumn(df, SEMTYPE.GENE_SYMBOL, ['gene name', 'gene symbol']);
if (geneCol) {
  sp.props.labelColumnNames = [geneCol.name];
  sp.props.showLabelsFor = 'MouseOverRow';
}
```

D-03 swaps the column binding to `Display Name` and adds a `showLabelsFor: 'Selected'` mode
driven by seeding `df.selection` with the top-N indices. Full pattern at RESEARCH.md §"Pattern 3":

```typescript
// EXTEND (preserve hover; add top-N selection-driven labels)
const labelCol = findColumn(df, SEMTYPE.DISPLAY_NAME, ['display name']) ??
                 findColumn(df, SEMTYPE.GENE_SYMBOL, ['gene name', 'gene symbol']);
if (labelCol) {
  sp.props.labelColumnNames = [labelCol.name];
  sp.props.showLabelsFor = 'Selected';  // platform RowSet enum value
  setTopNLabels(df, sp, options?.topNLabels ?? 15);  // seeds df.selection top-N
}
```

**Pitfall 7 from RESEARCH:** D-05 search-match ALSO writes `df.selection`. The two write
paths must union, not overwrite. `setTopNLabels` takes `mode: 'replace' | 'union'`;
default `'replace'` on metric/filter/color toggle; `'union'` when re-applying after a
search match.

#### D-06 live counter overlay (NEW — attach to viewer root)

Full code at RESEARCH.md §"Pattern 2: Live counter overlay anchored to a viewer".

Cleanup pattern: mirror `enrichment-viewers.ts:9` `activeSubscriptions` module-level array,
or store subscriptions in a `WeakMap<DG.ScatterPlotViewer, rxjs.Subscription[]>` (Pitfall 6
in RESEARCH).

```typescript
// PATTERN — enrichment-viewers.ts:9, 168-170
let activeSubscriptions: rxjs.Subscription[] = [];
// ... on re-entry: for (const s of activeSubscriptions) s.unsubscribe(); activeSubscriptions = [];
```

#### G2 dialog-state preload (extend `volcanoOptions` in package.ts:274-355)

Read `sp.getOptions()` at dialog-open time. Full pattern at RESEARCH.md §"Pattern 5". Critical
note from Pitfall 2: getOptions returns a snapshot; do not assume it is reactive.

#### G3 progress indicator (already in volcanoOptions, package.ts:337-352)

Already wraps `recomputeVolcano` with `DG.TaskBarProgressIndicator.create(...)` + `try/finally`.
No new code needed here — verify the `finally` clears on error (already does). The G3 fix is
purely a UX wording / immediate-feedback concern; the lifecycle is correct.

---

### 4. `src/package.ts` (MODIFY — extend in place)

#### D-05 + G4 — extend `dockComparisonFilterIfMultiContrast` (lines 117-131)

Current code uses `columnNames` only:

```typescript
// CURRENT (package.ts:117-131)
export function dockComparisonFilterIfMultiContrast(
  tv: DG.TableView, df: DG.DataFrame,
): boolean {
  const cmpCol = df.col('Comparison (group1/group2)') ?? df.col('Comparison');
  if (!cmpCol) return false;
  // ... distinct-count guard ...
  const filters = DG.Viewer.filters(df, {columnNames: [cmpCol.name]});
  tv.dockManager.dock(filters, DG.DOCK_TYPE.RIGHT, null, 'Comparison', 0.3);
  return true;
}
```

The G4 root-cause fix is the typed `filters` array, NOT the `columnNames` allowlist
(RESEARCH §"Pitfall 1"). Migrate to:

```typescript
// EXTEND — RESEARCH §"Pattern 1: Filters viewer with typed per-column filters"
const displayNameCol = findColumn(df, SEMTYPE.DISPLAY_NAME, ['display name']);
const sourceIdCol = findColumn(df, SEMTYPE.SOURCE_ID, ['source id']);
const filterSpecs: Array<{type: string; column: string}> = [
  {type: 'categorical', column: cmpCol.name},
];
if (displayNameCol) filterSpecs.push({type: 'free-text', column: displayNameCol.name});
if (sourceIdCol)    filterSpecs.push({type: 'free-text', column: sourceIdCol.name});

const filters = DG.Viewer.filters(df, {
  filters: filterSpecs,
  showBoolCombinedFilter: false,  // Flags exclusion — see Pitfall 1
  showHeader: true,
  showSearchBox: true,
});
tv.dockManager.dock(filters, DG.DOCK_TYPE.RIGHT, null, 'Filters', 0.3);
```

**Verification gate (RESEARCH §"Pattern 1"):** after dock, read
`viewer.getOptions().look.filters` and assert `Flags` is absent. If the platform still
auto-includes it, fall back to a post-create `viewer.setOptions({filters, columnNames})`.

#### D-12 new menu item — `Group-Mean Correlation…`

Mirror the existing `Visualize | Volcano Plot` handler at package.ts:262-272:

```typescript
// PATTERN — package.ts:262-272 (volcano launcher). Use `createX` factory convention.
@grok.decorators.func({'top-menu': 'Proteomics | Visualize | Group-Mean Correlation...'})
static async showGroupMeanCorrelation(): Promise<void> {
  const tv = grok.shell.tv;
  const df = tv?.dataFrame;
  if (!tv || !df) { grok.shell.warning('No table open'); return; }
  if (!requireDifferentialExpression(df, 'Run Differential Expression first ...')) return;
  const sp = createGroupMeanCorrelation(df);
  tv.addViewer(sp);
}
```

The correlation viewer stays on the protein DataFrame (no separate table view); D-12
explicitly notes this is different from PCA. See CLAUDE.md's "PCA opens in a separate
table view" note for the contrast.

---

### 5. `src/panels/uniprot-panel.ts` (MODIFY — extend `renderUniProtWidget`)

**Analog:** self (the existing GO Terms section, lines 152-173, is the template) + `qc-computations.ts:groupMean` (lines 15-25) for the math.

#### Section-header pattern (uniprot-panel.ts:153-157)

```typescript
// CURRENT — copy this idiom for the new "Per-Group Quantities" section
const goHeader = ui.divText('GO Terms');
goHeader.style.fontWeight = 'bold';
goHeader.style.marginTop = '8px';
goHeader.style.marginBottom = '4px';
container.appendChild(goHeader);
```

D-11 adds a `"Per-Group Quantities"` header in exactly the same shape, placed AFTER the
GO Terms block (lines 152-173), BEFORE the final `return container` (line 175).

#### Bar chart (full SVG body at RESEARCH §"Pattern 6: Inline SVG bar chart in a panel")

Key constraints from UI-SPEC §"UniProt panel: per-group quantity chart":
- SVG (no canvas), no new viewer.
- Two bars (one per group from `getGroups(df)`).
- Fills: `#FF00FF` (group1) / `#00FFFF` (group2) — LOCKED palette.
- Whiskers: mean ± 1 SD.
- Caption: `"{group} (n={N})"` below each bar.
- Numeric `"mean={x.xx}  SD={x.xx}"` text below caption at `fontSize: 0.85em` (matches
  `uniprot-panel.ts:170` font-size precedent).
- Empty state: `ui.divText('No per-group quantities available for this protein')`.

**Critical wiring detail (RESEARCH §"Pattern 6" footnote):** the panel receives only a
`proteinId` string, not the host DataFrame. Add a `findHostDataFrameForProtein(accession)`
helper that walks `grok.shell.tables` to find the protein DF (one whose `getGroups()` is
set AND whose primary-protein-id column resolves the accession). Preference: the DF on the
current `grok.shell.tv` if multiple match.

---

### 6. `src/analysis/enrichment.ts` (MODIFY — extend `showEnrichmentDialog` + insert R5 transform)

**Analog:** self.

#### Dialog checkbox (extend showEnrichmentDialog, lines 371-467)

Mirror the existing `ui.input.bool` lines (393-398) — six source checkboxes are already
the established input pattern in this dialog:

```typescript
// PATTERN — enrichment.ts:393-398
const goBpInput = ui.input.bool('GO: Biological Process', {value: true});
// ...

// NEW (after line 398, before the countDiv block at line 401)
const smartFilterInput = ui.input.bool('Apply smart pathway filter', {value: true});
smartFilterInput.setTooltip(
  'Drops generic GO parent terms when more specific child terms are present; ' +
  'caps each source at top-N by FDR. Recommended for cleaner results; ' +
  'uncheck for raw g:Profiler output.');
// ... then .add(smartFilterInput) inside the ui.dialog(...) chain
```

#### Smart filter transform (NEW — between gGOSt and buildEnrichmentDf at lines 344-345)

Full verbatim port at RESEARCH §"Pattern 7: Smart pathway filter — verbatim CK-omics port".
The locked constant lists are:

```typescript
// LOCKED CLIENT CONTRACT — port verbatim from ~/Downloads/ck/CKomics_tool2.py:4685-4736
const GENERIC_PARENT_TERMS = [
  'localization', 'cellular component organization', 'transport',
  'cellular process', 'biological process', 'metabolic process',
];
const SPECIFIC_CHILD_TERMS = ['actin', 'vesicle', 'endocytosis', 'cytoskeleton'];
```

**Invocation order (RESEARCH §"Pitfall 8"):** run `applySmartPathwayFilter` BEFORE
`buildEnrichmentDf` so the per-row Intersection strings don't reindex.

```typescript
// EXTEND enrichment.ts:343-350 (runEnrichmentPipeline)
if (upGenes.length > 0) {
  let upResults = await gGOSt(upGenes, bgArray, organismCode, sources, pThreshold);
  if (smartFilterEnabled) {
    const {kept, stats} = applySmartPathwayFilter(upResults, 15);
    upResults = kept;
    // capture stats so the banner can read them later
  }
  directionDfs.push(buildEnrichmentDf(upResults, upGenes, pThreshold, 'Up'));
}
```

#### Banner rendering (UI-SPEC §"Enrichment table: smart pathway filter banner")

Open Question 5 in RESEARCH recommends:

```typescript
// PATTERN — enrichment-viewers.ts:189-196 (existing dock pattern in the same view)
tv.dockManager.dock(
  ui.divText(bannerCopy),
  DG.DOCK_TYPE.TOP, null, 'Smart Filter Info', 0.1);
```

Set `df.setTag('proteomics.enrichment_smart_filtered', 'true')` when the filter ran
AND dropped rows (UI-SPEC §"Set tag").

---

### 7. Five parsers (MODIFY — one-line addition each)

**Analog:** self. Each parser has a single `return df` statement at the very end.

| File | Line | Existing return | Add before |
|------|------|-----------------|------------|
| `src/parsers/maxquant-parser.ts` | 158 | `return df;` | `await resolveGeneLabels(df);` |
| `src/parsers/spectronaut-parser.ts` | 185 | `return df;` | `await resolveGeneLabels(df);` |
| `src/parsers/spectronaut-candidates-parser.ts` | 267 | `return df;` | `await resolveGeneLabels(df);` |
| `src/parsers/fragpipe-parser.ts` | 196 | `return df;` | `await resolveGeneLabels(df);` |
| `src/parsers/generic-parser.ts` | dialog onOK (no synchronous `return df` — handler-driven) | n/a | call inside the onOK before `addTableView` |

**Critical (Pitfall 9):** every parser MUST call the resolver — even those whose typical
inputs (e.g. MaxQuant) rarely contain predicted IDs. The resolver itself is a no-op when
no predicted IDs are detected, but it ALWAYS adds the `Display Name` column so downstream
volcano label bindings stay consistent.

**Signature implication:** the parsers' return type changes from sync `DG.DataFrame` to
`Promise<DG.DataFrame>` for the four sync parsers. The package.ts import handlers (e.g.
`importMaxQuant` at lines 138-154) already use `async (file) => ...` for I/O; the
`await parseMaxQuantText(text)` change is local. Generic-parser is already inside an
async onOK, so it's only the synchronous form-of-call that changes.

Alternative: extract the resolver to be optionally invoked from the package.ts handlers
instead of inside the parsers. The planner picks. The RESEARCH recommendation is
"shared helper called from all 5 parsers" because it keeps the contract obvious — every
parser produces `Display Name` (Pitfall 9).

---

### 8. `src/utils/proteomics-types.ts` (MODIFY)

**Analog:** self.

#### Current shape (lines 1-9)

```typescript
export const SEMTYPE = {
  PROTEIN_ID: 'Proteomics-ProteinId',
  GENE_SYMBOL: 'Proteomics-GeneSymbol',
  LOG2FC: 'Proteomics-Log2FC',
  P_VALUE: 'Proteomics-PValue',
  INTENSITY: 'Proteomics-Intensity',
  SUBCELLULAR_LOCATION: 'Proteomics-SubcellularLocation',
} as const;
```

#### Add 4 new constants (D-08, D-12)

```typescript
DISPLAY_NAME: 'Proteomics-DisplayName',
SOURCE_ID: 'Proteomics-SourceId',
NUMERATOR_MEAN: 'Proteomics-NumeratorMean',
DENOMINATOR_MEAN: 'Proteomics-DenominatorMean',
```

**Convention from CLAUDE.md:** this file is the SINGLE source of truth for every
`Proteomics-*` string. The mirroring in `detectors.js` MUST track these constants byte-for-byte.

---

### 9. `detectors.js` (MODIFY — plain JS, mirror SEMTYPE additions)

**Analog:** self — every existing detector (lines 5-99) is a template.

#### detectGeneSymbol pattern (detectors.js:25-34)

```javascript
//meta.role: semTypeDetector
//input: column col
//output: string semType
detectGeneSymbol(col) {
  if (col.type !== DG.TYPE.STRING)
    return null;
  const name = col.name.toLowerCase();
  if (name === 'gene names' || name === 'gene name' || name === 'gene symbol' || name === 'gene') {
    col.semType = 'Proteomics-GeneSymbol';
    return col.semType;
  }
  return null;
}
```

#### Add 4 new detectors

| Function | Column name match | Sets semType |
|----------|-------------------|--------------|
| `detectDisplayName(col)` | `'display name'` exact (string col) | `'Proteomics-DisplayName'` |
| `detectSourceId(col)` | `'source id'` exact (string col) | `'Proteomics-SourceId'` |
| `detectNumeratorMean(col)` | `'numerator mean'` exact (float col) | `'Proteomics-NumeratorMean'` |
| `detectDenominatorMean(col)` | `'denominator mean'` exact (float col) | `'Proteomics-DenominatorMean'` |

**Constraint from CLAUDE.md:** `detectors.js` is plain JS (no ES module `import` syntax),
loaded by the platform BEFORE the webpack bundle. Mirror the existing function shape
EXACTLY — no TypeScript, no imports, no helper functions outside the class.

---

### 10. `src/tests/gene-label-resolver.ts` (NEW)

**Analog:** `src/tests/spectronaut-parser.ts` (lines 1-80) for the `category/test/expect`
shell + the fixture-builder helper pattern. No existing test mocks `grok.dapi.fetchProxy`;
the planner picks the mock approach (a stub on `grok.dapi` for the duration of the test,
restored in a `try/finally`).

#### Test-framework imports (spectronaut-parser.ts:1-7)

```typescript
// PATTERN
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {resolveGeneLabels} from '../utils/gene-label-resolver';
import {SEMTYPE} from '../utils/proteomics-types';
```

#### Tests required (RESEARCH §"Phase Requirements → Test Map")

- `geneLabelResolver` — Display Name populated; markers correct; predicted prefixes detected; LOC/RGD/AABR keep raw ID + `†`
- `geneLabelResolverCache` — `__schema_v` bump invalidates; cache miss → POST; cache hit → short-circuit
- `geneLabelEnsemblPost` — one POST per ≤1000 IDs; 429 with Retry-After → single retry
- `geneLabelDuplicates` — identical Display Names with different Source IDs → disambiguated; one warning toast

#### Mock pattern (no in-tree precedent; planner picks)

Recommendation: patch `grok.dapi.fetchProxy` to a stub inside each test, restore in
`finally`. Note that `grok.dapi.userDataStorage` returns a fresh instance per access
(`reference_dapi_fresh_instance_patch` in project memory) — patch the prototype, not the
instance.

---

### 11. `src/tests/group-mean-correlation.ts` (NEW)

**Analog:** `src/tests/parsers.ts` (lines 1-60) for the compact `category/test/expect`
shape. Use the `pcaDf`-style fixture pattern in `src/tests/analysis.ts` if it exists.

#### Tests required (RESEARCH §"Phase Requirements → Test Map")

- `groupMeanCorrelation` — known fixture (linear → r=1, ρ=1; uncorrelated → r≈0); rank-tie handling
- `groupMeanColumns` — Numerator Mean / Denominator Mean populated; semTypes set; re-run replaces via ensureFreshFloat

---

### 12. `src/tests/volcano.ts` (MODIFY — update assertions to new strings)

**Analog:** self.

#### Tests to update or add (RESEARCH §"Wave 0 Gaps")

- Update any literal `'up'`/`'down'`/`'not significant'` assertions to `'Enriched in ${g1}'`/`'Enriched in ${g2}'`/`'Not significant'`
- `volcanoDirectionColors` — assert color map contains `0xFFFF00FF`/`0xFF00FFFF`/`0xFFAAAAAA`
- `volcanoDirectionStrings` — assert categories derived from `getGroups(df)`
- `volcanoTopN` — top-15 rows by ascending metric + tiebreaker `|log2FC|` desc
- `volcanoCounter` — counter math: `df.filter` mask → category counts; `df.selection` change doesn't affect counts; metric toggle drives recompute

---

## Shared Patterns

### Pattern S-1: External REST + cross-session cache

**Source:** `src/analysis/subcellular-location.ts:165-393`
**Apply to:** `src/utils/gene-label-resolver.ts` (only new file using this pattern in Phase 14).

All seven pillars copy from the analog:

| Pillar | Excerpt source |
|--------|----------------|
| Store key + schema version | `subcellular-location.ts:169-173` |
| `loadCache` with `__schema_v` invalidation | `subcellular-location.ts:244-253` |
| `chunk(arr, size)` helper | `subcellular-location.ts:205-209` (private; copy or extract to shared util) |
| `runWithConcurrency` worker pool — REUSE | `subcellular-location.ts:219-237` (already exported for tests) |
| `FETCH_CONCURRENCY = 6` cap | `subcellular-location.ts:192` |
| Timer-driven flush every 5s + final flush in `finally` | `subcellular-location.ts:289-298, 384-388` |
| Warn-and-continue on per-batch failure | `subcellular-location.ts:312-325` |

### Pattern S-2: Bulk-init column writes (Phase 11 carry-forward)

**Source:** project memory `feedback_dg_column_bulk_init` (137-255× speedups measured on
the volcano / DE writeback paths). In-tree examples at
`volcano.ts:50-54` (ensureNegLog10Column), `volcano.ts:74-81` (ensureDirectionColumn),
`qc-computations.ts:53-54` (computeMA), `enrichment-viewers.ts:48-50` (negLogCol).

**Apply to:** every new column write in Phase 14 (`Display Name`, `Source ID`,
`Numerator Mean`, `Denominator Mean`, top-N label selection bitset).

```typescript
// PATTERN — never use per-row col.set(); always bulk-init.
const arr = new Float32Array(n);  // or string[] for string columns
for (let i = 0; i < n; i++) arr[i] = compute(i);
const col = df.columns.addNewFloat(name);  // or ensureFreshFloat for re-run safety
col.init((i) => arr[i]);
```

**Pitfall (`feedback_dg_column_init_null_sentinel`):** for FLOAT_NULL writes via init,
use `DG.FLOAT_NULL` explicitly. `col.init(() => null)` leaves the sentinel unsynced and
`get()` reads it as a finite positive — use `isNone()` / explicit FLOAT_NULL in the array.

### Pattern S-3: Find-column via semantic type, never raw name

**Source:** `src/utils/column-detection.ts:1-18` + CLAUDE.md "Conventions specific to this package".

```typescript
// PATTERN — always use findColumn / findProteomicsColumns
const labelCol = findColumn(df, SEMTYPE.DISPLAY_NAME, ['display name', 'gene name']);
// NEVER df.col('Display Name') without the SEMTYPE fallback path
```

**Apply to:** all new code reading the Display Name / Source ID / Numerator Mean /
Denominator Mean / direction / location columns. The mirrored detectors in `detectors.js`
make the SEMTYPE path the primary lookup; the name-hint array is the fallback for
DataFrames where the detector hasn't run.

### Pattern S-4: rxjs subscription cleanup on viewer re-entry

**Source:** `src/viewers/enrichment-viewers.ts:9, 167-170`.

```typescript
// PATTERN — module-level subscriptions array, cleared on every re-entry.
let activeSubscriptions: rxjs.Subscription[] = [];

// at re-entry:
for (const sub of activeSubscriptions) sub.unsubscribe();
activeSubscriptions = [];

// add new subscriptions:
activeSubscriptions.push(sub1, sub2, sub3);
```

**Apply to:** the volcano counter overlay (D-06) and any new cross-DataFrame wiring.
Pitfall 6 in RESEARCH calls out the consequence of skipping this: stacked overlays + frozen
counters when the user re-opens the viewer.

### Pattern S-5: TaskBarProgressIndicator for long async paths

**Source:** `src/package.ts:337-352` (volcanoOptions) and `src/package.ts:365-371` (heatmap).

```typescript
// PATTERN
const pi = DG.TaskBarProgressIndicator.create('<label>');
try {
  await someLongAsyncWork(progress);
} catch (e: any) {
  grok.shell.error(`<op> failed: ${e?.message ?? e}`);
} finally {
  pi.close();
}
```

**Apply to:** R1 resolver invocation in parsers (long Ensembl fetch on cold cache),
G3 Color → Location (already wrapped — just verify the cache-aware pre-OK toast wording).

### Pattern S-6: ensureFreshFloat for re-run-safe derived columns

**Source:** `src/viewers/qc-computations.ts:28-32`.

```typescript
function ensureFreshFloat(df: DG.DataFrame, name: string): DG.Column {
  if (df.columns.contains(name))
    df.columns.remove(name);
  return df.columns.addNewFloat(name);
}
```

**Apply to:** `Numerator Mean` / `Denominator Mean` in the new correlation viewer (D-12).
The helper is currently file-private to `qc-computations.ts`; the planner picks between
copying it into `group-mean-correlation.ts` and lifting it to `src/utils/`.

---

## No Analog Found

None. Every Phase 14 capability has at least one in-tree precedent.

The one "soft" gap is the **mocked-fetchProxy test pattern** for the new
`gene-label-resolver.ts` test file: no existing test in `src/tests/` patches
`grok.dapi.fetchProxy`. RESEARCH §"Wave 0 Gaps" flags this; planner picks the mock approach
(stub on `grok.dapi` for test duration, restore in `finally`).

---

## Metadata

**Analog search scope:**
- `src/analysis/` (5 files) — found cache + worker-pool analog in `subcellular-location.ts`
- `src/viewers/` (7 files) — found viewer-factory analog in `volcano.ts`, derived-column analog in `qc-computations.ts`, title-annotation analog in `pca-plot.ts`, subscription-cleanup analog in `enrichment-viewers.ts`
- `src/panels/` (1 file) — found section-header analog in `uniprot-panel.ts` itself
- `src/parsers/` (5 files) — confirmed every parser ends with `return df` at a clearly identifiable line; 4 lines need a single `await resolveGeneLabels(df)` insertion
- `src/utils/` (2 files) — `proteomics-types.ts` and `column-detection.ts` are the two single-sources-of-truth this phase extends
- `src/tests/` (7 files) — found `category/test/expect` shape in every test file
- `detectors.js` (root) — found 6 existing detector patterns; new SEMTYPEs need 4 mirror detectors

**Files scanned:** 21 source files + 1 root JS detector file + 4 phase-14 planning docs.

**Pattern extraction date:** 2026-05-28

**Confidence:** HIGH — every excerpt is verified against the in-tree file at the cited
line range; cross-references to RESEARCH.md and 14-UI-SPEC.md were re-checked against
14-CONTEXT.md D-01..D-14.
