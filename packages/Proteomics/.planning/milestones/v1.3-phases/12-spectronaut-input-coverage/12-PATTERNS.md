# Phase 12: Spectronaut Input Coverage - Pattern Map

**Mapped:** 2026-05-15
**Files analyzed:** 7 (2 modified, 1 new module, 4 new/extended assets)
**Analogs found:** 7 / 7 (all in-repo; no RESEARCH.md-only fallbacks needed)

> Every analog cited below is a real file in `packages/Proteomics/`. The phase is
> **additive** (new input path); the highest-leverage pattern is *reuse the
> existing pivot/tag/log2 tail verbatim* and only replace the text→rows front.

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|----------------|---------------|
| `src/parsers/spectronaut-parser.ts` (MODIFY: extract `finalizeSpectronaut`, add `parseSpectronautStream`, add `aggToPivotResult`) | parser | streaming + transform | itself (`parseSpectronautText` / `pivotSpectronaut` / `buildWideDataFrame`) | exact (in-file refactor) |
| `src/package.ts` (MODIFY: `importSpectronaut` header-sniff branch ~L99-115) | controller (menu handler) | file-I/O / request-response | `importSpectronaut` (current) + `importMaxQuant` (L81-97) | exact |
| `src/tests/spectronaut-parser.ts` (EXTEND `makeLongFormatTsv`; add streaming/golden/parity tests) | test | transform | itself (existing `category('Spectronaut')` suite) | exact |
| `tools/spectronaut-aggregate.sql` (NEW, from `/tmp`) | tooling (dev/CI) | batch / transform | `tools/generate-spectronaut-candidates-fixture.mjs` | role-match |
| `tools/spectronaut-aggregate.sh` (NEW, from `/tmp`) | tooling (dev/CI) | batch | `/tmp/spectronaut-aggregate.sh` (move as-is) + fixture-gen one-liner convention | exact (relocate) |
| `files/demo/spectronaut-hye-precursor.tsv` + `…-precursor-golden.tsv` (NEW fixtures) | fixture (data asset) | file-I/O | `files/demo/spectronaut-hye-candidates.tsv` (derived fixture) | role-match |
| `files/demo/README.md` (EXTEND: 2 new fixture rows + regen one-liner + flip caveat) | docs | — | existing `spectronaut-hye-candidates.tsv` section (L76-93) | exact |

---

## Pattern Assignments

### `src/parsers/spectronaut-parser.ts` (parser, streaming + transform)

**Analog:** itself — the streaming path reuses the entire post-pivot tail; only the
front (text→`DG.DataFrame`→`pivotSpectronaut`) is replaced by a streaming aggregator
that emits the *same* `PivotResult`.

**Imports pattern** (current L1-7 — keep, the new code needs no new imports):
```typescript
import * as DG from 'datagrok-api/dg';
import {SEMTYPE} from '../utils/proteomics-types';
import {
  log2TransformColumns, copyAsLog2Columns, detectLog2Status,
  addPrimaryColumnIfNeeded,
} from './shared-utils';
import {setGroups} from '../analysis/experiment-setup';
```

**`PivotResult` contract — the seam** (L17-23, unchanged; the streaming aggregator's
output target):
```typescript
interface PivotResult {
  proteinMap: Map<string, Map<string, number>>;
  sampleKeys: string[];
  sampleFileNames: Map<string, string>;
  organisms: Map<string, string>;
}
```

**Tail to EXTRACT into `finalizeSpectronaut(result: PivotResult)` — lift L174-203
verbatim** (this is the single most important reuse; both paths must call it so
output is byte-identical — addresses Pitfall 2 / R3):
```typescript
// currently inline at the end of parseSpectronautText (L174-203):
const df = buildWideDataFrame(result);
addPrimaryColumnIfNeeded(df, 'PG.ProteinGroups', 'Primary Protein ID', SEMTYPE.PROTEIN_ID);
const intensityColNames = result.sampleKeys;
const log2Status = detectLog2Status(df, intensityColNames);
if (log2Status.isLog2) copyAsLog2Columns(df, intensityColNames);
else log2TransformColumns(df, intensityColNames);
df.setTag('proteomics.preNormalized', 'true');   // tag contract — DO NOT change
df.setTag('proteomics.source', 'spectronaut');    // tag contract — DO NOT change
autoPopulateGroups(df, result.sampleKeys);
return df;
```
After extraction, `parseSpectronautText` ends with `return finalizeSpectronaut(result);`
and `parseSpectronautStream` ends with `return finalizeSpectronaut(aggToPivotResult(agg));`.

**Bulk wide-column construction pattern — KEEP as-is in `buildWideDataFrame`**
(L106-117; the established `Float32Array` + `Column.fromFloat32Array` idiom — never
regress to per-row `col.set`, per `feedback_dg_column_bulk_init` / commit `a7ae0e8f42`):
```typescript
for (const sampleKey of result.sampleKeys) {
  const values = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    const sampleMap = result.proteinMap.get(proteins[i])!;
    values[i] = sampleMap.has(sampleKey) ? sampleMap.get(sampleKey)! : DG.FLOAT_NULL;
  }
  const col = DG.Column.fromFloat32Array(sampleKey, values);
  col.semType = SEMTYPE.INTENSITY;
  if (result.sampleFileNames.has(sampleKey))
    col.setTag('spectronaut.fileName', result.sampleFileNames.get(sampleKey)!);
  cols.push(col);
}
```

**`sampleKey` format — MUST match `pivotSpectronaut` L65 exactly** (`aggToPivotResult`
has to produce identical keys or `autoPopulateGroups`'s `key.lastIndexOf('_')` split
breaks):
```typescript
const sampleKey = `${condition}_${replicate}`;   // pivotSpectronaut L65 — single source of truth
```

**Filter semantics — the streaming aggregator must reproduce `pivotSpectronaut`
L47-78** (q-value: numeric `> threshold` excluded, non-numeric/null passes; drop
`!protein`, `CON__`/`REV__`; quantity-column resolution order `PG.IBAQ` → `PG.Quantity`
from `QUANTITY_COLUMNS` L15):
```typescript
const QUANTITY_COLUMNS = ['PG.IBAQ', 'PG.Quantity'];   // L15 — resolve ONCE from header, not per row
// L60-61:
if (!protein) continue;
if (protein.startsWith('CON__') || protein.startsWith('REV__')) continue;
```
> Parity note (RESEARCH Pitfall 4): `pivotSpectronaut` uses loose `Number(raw)`/`isNaN`;
> the streaming path must match **duckdb `TRY_CAST`/`nullstr`** (the D-04 oracle) — add a
> `tryCastDouble(s): number|null` helper (`''`/`'NaN'`/`'NA'`/null → null). The two
> differ only on degenerate inputs the real/fixture data does not contain.

**Required-column validation pattern** (L159-169 — replicate header-side in the stream
once the header line is parsed):
```typescript
for (const colName of REQUIRED_COLUMNS) {           // ['R.Condition','R.Replicate','PG.ProteinGroups']
  if (!longDf.col(colName)) throw new Error(`Missing required Spectronaut column: ${colName}`);
}
const quantityColName = QUANTITY_COLUMNS.find((n) => longDf.col(n) !== null);
if (!quantityColName) throw new Error(`Missing protein-group quantity column ...`);
```

**Streaming front (NEW) — model on RESEARCH Pattern 1 + Code Examples**, using the
package's progress + explicit-yield idiom (see `differential-expression.ts` analog
below for the `TaskBarProgressIndicator` create/try/finally shape; the *explicit
`setTimeout(0)` macrotask yield* is the new piece — RESEARCH Pitfall 3).

---

### `src/package.ts` (controller / menu handler, file-I/O)

**Analog:** the current `importSpectronaut` (L99-115) — the branch is **additive**;
the `no` branch must remain byte-identical. `importMaxQuant` (L81-97) is the sibling
showing the identical `DG.Utils.openFile` → `try/catch` → `df.name` → `addTableView`
→ `grok.shell.info`/`error` shape every importer follows.

**Handler skeleton to preserve (current L99-115)** — keep the decorator, `openFile`,
`accept`, naming, success/error toasts; only swap the body of `open`:
```typescript
@grok.decorators.func({'top-menu': 'Proteomics | Import | Spectronaut Report...'})
static async importSpectronaut(): Promise<void> {
  DG.Utils.openFile({
    accept: '.tsv,.txt,.csv',
    open: async (file: File) => {
      try {
        // D-01 header-sniff branch goes HERE (replaces the two `file.text()` lines):
        //   const df = (await sniffIsPrecursor(file))
        //     ? await parseSpectronautStream(file)
        //     : parseSpectronautText(await file.text());   // ← UNCHANGED legacy path
        df.name = file.name.replace(/\.[^.]+$/, '');
        grok.shell.addTableView(df);
        grok.shell.info(`Imported ${df.rowCount} protein groups from Spectronaut`);
      } catch (e: any) {
        // D-05 failure-path hint goes HERE — point at tools/spectronaut-aggregate.sh
        grok.shell.error(`Failed to import Spectronaut file: ${e?.message ?? e}`);
      }
    },
  });
}
```

**Failure-path hint pattern** — the package already has a precedent for a multi-line
actionable `grok.shell.warning` that names the exact menu/command to use
(`requireSampleLevelData`, package.ts L61-73). Mirror its tone for the D-05 hint
pointing at `tools/spectronaut-aggregate.sh` / the README:
```typescript
grok.shell.warning(
  `${action} needs per-sample intensities, but this table is a Spectronaut Candidates ` +
  `report ... Import the matching Spectronaut Report (Proteomics | Import | Spectronaut Report) ...`,
);
```

**Header-sniff helper (NEW)** — model on RESEARCH "Code Examples / Header-sniff": one
`file.stream().pipeThrough(new TextDecoderStream('utf-8')).getReader()`, read until
first `\n`, `await reader.cancel()` in `finally`, split header on `\t`, test
`['EG.ModifiedPeptide','FG.Charge','PEP.StrippedSequence'].some(...)`. No in-repo
analog for the stream-sniff itself — RESEARCH is the source (citation HIGH).

---

### `src/tests/spectronaut-parser.ts` (test, transform)

**Analog:** itself — the existing `category('Spectronaut')` suite (19 tests). New
tests slot into the same file/category; `makeLongFormatTsv` is *extended in place*
per D-03 (it already emits `PEP.StrippedSequence` and 2 peptide rows per
protein×sample — L33-37).

**`makeLongFormatTsv` shape to extend** (L9-42 — add `EG.ModifiedPeptide`/`FG.Charge`/
`FG.Id` columns + multi-precursor/`CON__`/`REV__`/sub-threshold/non-numeric/empty
`EG.Qvalue` rows; keep `CondA`/`CondB` defaults so the duckdb DMD↔WT flip is a
structural no-op — RESEARCH Pitfall 1):
```typescript
function makeLongFormatTsv(opts: {
  proteins: {id: string; organism?: string}[];
  conditions: string[];
  replicates: number[];
  qValues?: Map<string, number | string>;
  ibaqValues?: Map<string, number>;
  fileNameTemplate?: string;
  quantityColumn?: string;            // 'PG.IBAQ' default; 'PG.Quantity' alias
}): string {
  const headers = ['R.FileName','R.Condition','R.Replicate','PG.ProteinGroups',
    'PG.Organisms', quantityCol, 'EG.Qvalue', 'PEP.StrippedSequence'];
  // ... two peptide rows per protein+sample to test deduplication ...
}
```
> Keep carry-along columns (`PG.Organisms`, `R.FileName`) **constant within each
> (protein×condition×replicate) group** so duckdb `any_value()` == streaming
> first-non-null (RESEARCH Pitfall 5). The hye-mix demo has **no** `PG.Genes`/
> `PG.ProteinAccessions` columns — the fixture/SQL/test must not assume them.

**Test structure pattern** (the `category`/`test`/`expect` idiom — L44-58, mirror for
all new tests: streaming-vs-text equivalence, streaming-vs-golden, filter-branch
parity, tag set, tools-file presence):
```typescript
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {parseSpectronautText} from '../parsers/spectronaut-parser';
// add: import {parseSpectronautStream} from '../parsers/spectronaut-parser';

category('Spectronaut', () => {
  const baseOpts = {proteins: baseProteins, conditions: ['CondA','CondB'], replicates: [1,2]};

  test('pivot produces correct protein count', async () => {
    const tsv = makeLongFormatTsv(baseOpts);
    const df = parseSpectronautText(tsv);
    expect(df.rowCount, 3);
  });
});
```

**Tag/semType/group assertion patterns to reuse for the new equivalence test**
(L109-143, 183-187 — assert the streaming DataFrame matches these exactly):
```typescript
expect(df.col('PG.ProteinGroups')?.semType, SEMTYPE.PROTEIN_ID);   // L112
expect(df.col('CondA_1')?.semType, SEMTYPE.INTENSITY);             // L118
expect(df.getTag('proteomics.source'), 'spectronaut');             // L186
expect(df.getTag('proteomics.preNormalized'), 'true');             // L159
const groups = getGroups(df);                                      // L137
expect(groups!.group1.columns.length, 2);
```

> RESEARCH Open Q1 recommendation: **keep the 19 existing tests on
> `parseSpectronautText`** (locks the text path) and **add** a focused test asserting
> `parseSpectronautStream(fixture)` produces an equal DataFrame (rows, sample cols,
> tags, groups) — do not rewrite the existing suite.

---

### `tools/spectronaut-aggregate.sql` + `tools/spectronaut-aggregate.sh` (tooling, batch/transform)

**Analog:** `/tmp/spectronaut-aggregate.{sql,sh}` (move verbatim) + the established
`tools/` precedent `tools/generate-spectronaut-candidates-fixture.mjs` (committed dev
generator with a documented one-liner + README row).

**`.sh` already references a sibling SQL by `dirname` (L16)** — moving both into
`tools/` keeps the wrapper working unchanged:
```bash
SQL="$(dirname "$0")/spectronaut-aggregate.sql"
duckdb < "$TMP_SQL"            # L35
```

**The `R.Condition` DMD↔WT flip (SQL L18-25) is reference-file-specific — KEEP it in
the committed `tools/` SQL but document it LOUDLY as reference-only** (RESEARCH
Pitfall 1, the single highest-risk parity trap — the streaming aggregator must NOT
port it; the synthetic fixture's `CondA`/`CondB` make it a no-op so the same SQL
serves as both fallback and D-04 oracle):
```sql
-- Flip mis-labeled R.Condition values. Cross-tab against R.FileName showed
-- every WT* filename tagged DMD and every DMD* filename tagged WT (24/24);
CASE "R.Condition" WHEN 'DMD' THEN 'WT' WHEN 'WT' THEN 'DMD' ELSE "R.Condition" END
```

**Aggregate/filter logic the streaming TS must match** (SQL L28-49 — `max(TRY_CAST
PG.Quantity)`, `min(TRY_CAST EG.Qvalue)`, q≤0.01-or-null, `CON__`/`REV__`/null-PG
drop, `GROUP BY PG, Condition, Replicate`, `nullstr=['','NaN','NA']`,
`ignore_errors=true`).

**`tools/` generator conventions to mirror** (`generate-spectronaut-candidates-fixture.mjs`
L1-25, 207-209): top-of-file doc block with a `Run: node tools/...` line, `PKG_ROOT`
resolution, writes into `files/demo/`, `console.log` summary. The D-04 regen one-liner
is the duckdb equivalent: `tools/spectronaut-aggregate.sh files/demo/<fixture>.tsv files/demo/<golden>.tsv`.

**`package.json` scripts convention** — there is no per-fixture npm target today
(fixture gen is a documented bare command, not an npm script). Follow that: document
the regen one-liner in `files/demo/README.md`, do **not** add an npm script unless the
planner decides otherwise (existing scripts L are build/lint/test/publish only).

---

### `files/demo/spectronaut-hye-precursor.tsv` + `…-precursor-golden.tsv` (fixture, file-I/O)

**Analog:** `files/demo/spectronaut-hye-candidates.tsv` — a **derived** committed
fixture generated by a `tools/` script, documented with a "Regenerate with …" line
and a `**Derived**`/license note. The new pair follows the same lifecycle: synthetic
precursor TSV (from extended `makeLongFormatTsv` logic or a small `tools/` emitter) +
its duckdb-derived golden `.tsv`, regenerated together via the one-liner so they
never drift (RESEARCH "Runtime State Inventory").

**Fixture column shape** must carry the D-01 precursor signature so header-sniff
routes it to the streaming path (from hye-mix header positions, verified): include
`EG.ModifiedPeptide`, `FG.Charge`, `PEP.StrippedSequence` plus `R.FileName`,
`R.Condition`, `R.Replicate`, `PG.ProteinGroups`, `PG.Organisms`, and a quantity
column. **Use `PG.IBAQ` or `PG.Quantity` consistently** (hye-mix uses `PG.IBAQ` at
col 11; the streaming resolver prefers `PG.IBAQ`). Keep conditions `CondA`/`CondB`
(not `DMD`/`WT`) so the committed-SQL flip is a no-op.

---

### `files/demo/README.md` (docs)

**Analog:** the existing `## spectronaut-hye-candidates.tsv` section (L76-93) — the
exact table + prose template for a derived fixture, including the regen one-liner and
the "not a verbatim export / values are ours" caveat:
```markdown
| **Source** | **Derived** from `spectronaut-hye-mix.tsv` via this package's own
Welch's t-test + BH FDR. Regenerate with `node tools/generate-spectronaut-candidates-fixture.mjs`. |
```
Add two analogous rows (precursor fixture + golden) and a **flip caveat** noting the
committed `tools/spectronaut-aggregate.sql` carries a reference-file-only DMD↔WT
correction that is a structural no-op on this `CondA`/`CondB` fixture. Append matching
**License** lines (L95-101 pattern: synthetic → "no restrictions").

---

## Shared Patterns

### Progress indicator + close-in-finally
**Source:** `src/analysis/differential-expression.ts:339-403`
**Apply to:** `parseSpectronautStream` (and the header-sniff if it can be slow)
```typescript
const pi = DG.TaskBarProgressIndicator.create(`Running ${method} analysis...`);
try {
  // ... work; periodic pi.update(pct, message) ...
} finally {
  pi.close();
}
```
> The DE precedent shows the create/try/finally lifecycle and `pi.update(pct, msg)`
> cadence. It does **not** show an event-loop yield (commit `a7ae0e8f42` was a
> bulk-init fix, not a yield — RESEARCH Pitfall 3). The streaming path must add an
> **explicit macrotask yield** `await new Promise((r) => setTimeout(r, 0))` every N
> lines/~16 ms, alongside `pi.update`. There is no in-repo yield analog; RESEARCH
> (CONCERNS.md `setTimeout(0)` idiom) is the source.

### Importer handler shape (openFile → parse → name → addTableView → toast)
**Source:** `src/package.ts` `importMaxQuant` L81-97, `importSpectronaut` L99-115,
`importSpectronautCandidates` L117-133 (three identical instances)
**Apply to:** the modified `importSpectronaut` — preserve this exact skeleton; the
header-sniff is the only insert.

### `proteomics.*` tag contract (set on completion, gated downstream)
**Source:** `spectronaut-parser.ts:195-198`
**Apply to:** `finalizeSpectronaut` (the single place both paths set tags)
```typescript
df.setTag('proteomics.preNormalized', 'true');
df.setTag('proteomics.source', 'spectronaut');
```
Per CLAUDE.md tag table: streaming path sets the **same** `proteomics.source=spectronaut`,
`proteomics.preNormalized=true`, and auto-populates `proteomics.groups` via
`autoPopulateGroups` → `setGroups` (experiment-setup.ts) — no new tag introduced.

### Bulk `Column.init` / `fromFloat32Array` (no per-row `col.set`)
**Source:** `buildWideDataFrame` L106-117; `shared-utils.ts` `log2TransformColumns`
L15-22 / `copyAsLog2Columns` L38-39 / `addPrimaryColumnIfNeeded` L62-67 (all use
`col.init((i) => …)` / `Float32Array`)
**Apply to:** all DataFrame construction on the streaming path — it already routes
through `buildWideDataFrame`, so this is satisfied by reuse; do not add row loops.
(`feedback_dg_column_bulk_init`, commit `a7ae0e8f42`.)

### `findColumn`/`SEMTYPE` convention
**Source:** `src/utils/proteomics-types.ts` `SEMTYPE.*`; parser uses
`SEMTYPE.PROTEIN_ID` / `SEMTYPE.INTENSITY` (spectronaut-parser.ts L95, 113)
**Apply to:** any semType assignment in new code — never inline `'Proteomics-*'`
strings. (Note: the parser reads Spectronaut columns by fixed vendor name
`R.Condition` etc. — that is correct *here* because it is a vendor-specific parser,
not generic pipeline code; `findColumn` applies to pipeline/viewer code, not the
vendor parser's own required-column lookup.)

---

## No Analog Found

| File / concern | Role | Data Flow | Reason |
|----------------|------|-----------|--------|
| Browser stream loop (`Blob.stream()` + `TextDecoderStream` + carry-over buffer) | parser | streaming | No streaming reader exists anywhere in the package (every importer uses `await file.text()`). Use RESEARCH Pattern 1 + Code Examples (HIGH confidence, WHATWG Baseline APIs). |
| Explicit `setTimeout(0)` macrotask yield in a hot loop | parser | streaming | DE "unblock" precedent was a bulk-init fix, not a yield (RESEARCH Pitfall 3). No copyable in-repo yield; RESEARCH (CONCERNS.md kNN/QC idiom) is the source. |
| Header-sniff via partial stream read + `reader.cancel()` | controller | file-I/O | No partial-read precedent; RESEARCH "Code Examples / Header-sniff" is the source. |

> These three are confined to the **new streaming front**. The pivot/tag/log2 tail,
> the importer handler, the test harness, and the tooling/README all have exact
> in-repo analogs and must reuse them — that is where the regression risk lives.

## Metadata

**Analog search scope:** `packages/Proteomics/src/{parsers,analysis,tests}/`,
`packages/Proteomics/{package.ts}`, `packages/Proteomics/tools/`,
`packages/Proteomics/files/demo/`, `/tmp/spectronaut-aggregate.{sql,sh}`,
`packages/Proteomics/package.json`
**Files scanned:** 11 (spectronaut-parser.ts, shared-utils.ts, package.ts §imports,
spectronaut-parser test, differential-expression.ts §progress,
generate-spectronaut-candidates-fixture.mjs, files/demo/README.md, demo dir listing,
hye-mix header, /tmp .sql, /tmp .sh)
**Pattern extraction date:** 2026-05-15
