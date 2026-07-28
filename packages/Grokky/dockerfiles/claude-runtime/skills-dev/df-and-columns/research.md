# DataFrame & Column API research — for `datagrok-df-and-columns` skill

Citations use `<absolute-path>:<line>` form. All `.ts` paths are under
`/Users/aleksashka/Desktop/datagrok/reddata/public/js-api/src/` (abbreviated `js-api/`)
and `/Users/aleksashka/Desktop/datagrok/reddata/public/packages/ApiSamples/scripts/` (abbreviated `samples/`).

## 1. Summary table

| # | Area | Key entrypoints | Canonical idiom | Wrapper-worthiness |
|---|------|-----------------|-----------------|--------------------|
| A | Finding columns | `df.col`, `df.getCol`, `df.columns.byName/byIndex/byTags/bySemType/bySemTypeAll/firstWhere/numerical/categorical/dateTime/numericalNoDateTime/boolean/names/contains/toMap` | `df.col(name)` for nullable, `df.getCol(name)` to throw; iterate `df.columns.numerical` etc.; `bySemType(DG.SEMTYPE.MOLECULE)` | **High** — picking the "right" col when several plausible exist |
| B | Column metadata | `col.semType`, `col.tags`, `col.getTag/setTag`, `col.meta.{friendlyName,description,format,units,formula,choices,cellRenderer}`, `col.meta.colors`, `col.meta.markers`, `DG.SEMTYPE`, `DG.TAGS`, `DG.COLUMN_TYPE` | `col.meta.units = 'kg'`, `col.meta.format = '0.00'`, `col.semType = DG.SEMTYPE.MOLECULE`; **never** `setTag` for color coding | **High** — Claude doesn't know about `col.semType`, `meta.colors.setLinear` |
| C | Column stats | `col.stats.{min,max,avg,stdev,med,q1,q2,q3,sum,valueCount,missingValueCount,uniqueCount,totalCount,variance,skew,kurt,corr}`, `col.min`, `col.max`, `col.categories`, `col.toList()`, `col.values()`, `col.getRawData()` | `col.stats.avg`, iterate `col.categories` (string only) | **Medium** — already concise, mostly teaching |
| D | Cloning DataFrames | `df.clone(rowMask?, columnIds?, saveSelection?, saveTags?)`, `col.clone(mask?)` | `df.clone(null, ['a','b'])` for column subset; `df.clone(df.filter)` for filtered copy. `df.copy()` does **not** exist | **Medium** — wrapper to clarify args |
| E | Adding columns | `df.columns.addNew(name,type)`, `addNewInt/Float/String/Bool/DateTime/Qnum/Bytes`, `addNewCalculated(name, formula, type?)`, `addNewVirtual`, `add(col)`, `insert(col, index)`, `getOrCreate(name, type)`, `replace`, `getUnusedName` | typed shorthand + `.init(scalar)` or `.init(i => …)`; `addNewCalculated` is async | **High** — picking right variant |
| F | Removing / renaming | `df.columns.remove(name|index|Column)`, `col.name = 'X'`, `getUnusedName` | rename: just assign `col.name`; remove: `df.columns.remove('x')` (chainable; case-insensitive) | **Medium** — wrapper around bulk remove + name collisions |
| G | Setting values | `col.set(i, v, notify?)`, `col.init(scalar|fn)`, `col.fireValuesChanged()`, `col.getRawData()/setRawData()`, `col.setString` | bulk: `col.init(i => …)`; performance: write to `getRawData()` directly + `fireValuesChanged` | **Low** — well-shaped already |
| H | Rows | `df.rows.removeAt`, `removeWhere`, `removeWhereIdx`, `insertAt`, `addNew(values?)`, `setValues`, `match`, `select`, `filter`, `df.row(i).get/set` | row deletion: `df.rows.removeWhere(r => r.get('make')==='Honda')`. **Anti-pattern**: looping `df.rows` to read/write column data | **Low** — but worth warning about |
| I | Current state on DF | `df.currentRow`/`currentRowIdx`, `df.currentCol`, `df.currentCell`, `df.selection` (BitSet), `df.filter` (BitSet), `df.mouseOverRowIdx` | `df.currentRowIdx` for index, `df.selection.trueCount` for selected-row count | **Medium** — covered in separate "state" skill |
| J | Existing grokky helpers | `addCalculatedColumn`, `filterRows`, `sortRows`, `selectRows`, `aggregateBy`, `pivot`, `unpivot`, `joinTables`, `addViewer`, `configureViewer` | only `addCalculatedColumn` overlaps with this skill | **N/A** — see §10 for gaps |

## 2. Finding columns (Area A)

Entrypoints live in `ColumnList`:

- `byName(name)` — case-insensitive lookup, returns `Column` or `null`-like dart wrapper (`js-api/dataframe/column-list.ts:40`)
- `byIndex(i)` (`column-list.ts:46`)
- `bySemType(semType)` returns first match or `null` (`column-list.ts:49`)
- `bySemTypeAll(semType)` returns array (`column-list.ts:55`)
- `bySemTypesExact(semTypes[])` returns array or `null` if any missing (`column-list.ts:60`)
- `byTags({tag: value})` — `null`/`undefined` value matches any (`column-list.ts:74`, sample: `samples/data-frame/find-columns.js:17`)
- `firstWhere(predicate)` (`column-list.ts:79`)
- `numerical` / `categorical` / `dateTime` / `numericalNoDateTime` / `boolean` / `selected` getters (`column-list.ts:92-115`)
- `names()`, `toList()`, `toMap()` (`column-list.ts:119,129,134`)
- `contains(name)` case-insensitive (`column-list.ts:257`)
- Top-level shortcuts: `df.col(nameOrIndex)` returns nullable (`data-frame.ts:210`); `df.getCol(name)` throws on missing (`data-frame.ts:230`).

Idiom from `samples/data-frame/find-columns.js:5-21`:

```js
for (let column of demog.columns.categorical) ...
for (let column of demog.columns.byTags({'tag1': 'value1', 'tag2': 'value2'})) ...
```

Iteration via `for (const col of df.columns)` works because `ColumnList` implements `Symbol.iterator` (`column-list.ts:276`).

**Gotchas**:
- `byName` and `contains` are **case-insensitive**; this can mask bugs when user passes camelCase but DF has PascalCase.
- `bySemType` returns first match — Claude must call `bySemTypeAll` when multiple molecule columns exist.
- `byTags({k: undefined})` matches any value (this is a *feature*, easily confused with "tag not present").
- No fuzzy/Levenshtein matching out of the box — Claude needs help when user says "mw" but column is "Mol Weight".

**Wrapper rationale**: `findColumns(df, query)` that combines name (fuzzy), semType, type, and tag filters, returning ranked candidates. Reduces Claude picking the wrong "mass-y" column when there are three.

## 3. Column metadata (Area B)

Per-column metadata splits across three places:

**Direct `Column` properties** (`js-api/dataframe/column.ts`):
- `col.type` (read-only `ColumnType`, `column.ts:188`)
- `col.semType` get/set string (`column.ts:212`/`216`)
- `col.name` get/set (`column.ts:231`/`235`)
- `col.layoutColumnId` (`column.ts:222`)
- `col.isNumerical`, `col.isCategorical`, `col.isVirtual` (`column.ts:176-195`)
- `col.tags` is a `MapProxy` — supports `tags.foo`, `tags['x']`, `tags.has`, `tags.delete`, `tags.set`, `tags.clear`, iteration (`column.ts:33`; sample: `samples/data-frame/metadata/metadata.js:38-55`).
- `col.getTag/setTag` direct API (`column.ts:369`/`377`)
- `col.temp` — non-serialized scratch space (`column.ts:32`).

**`col.meta` helper** (`js-api/dataframe/column-helpers.ts:178`):
- `friendlyName`, `description`, `format`, `formula`, `units`, `cellRenderer` — all string get/set backed by the corresponding `DG.TAGS.*` (`column-helpers.ts:213-238`).
- `format` falls back to `grok_Column_GetAutoFormat` if unset (`column-helpers.ts:222`).
- `choices` / `autoChoices` (string-typed combo editor) — `choices` is JSON-encoded under the hood (`column-helpers.ts:243-252`).
- `multiValueSeparator`, `includeInCsvExport`, `includeInBinaryExport`, `allowColorPicking`, `linkClickBehavior` (`column-helpers.ts:254-274`).
- Setting any of these to `null` removes the underlying tag via `setNonNullTag` (`column-helpers.ts:206`).

**`col.meta.colors`** = `ColumnColorHelper` (`column-helpers.ts:48`):
- `setLinear(range?, {min, belowMinColor, max, aboveMaxColor})` (`column-helpers.ts:75`)
- `setLinearAbsolute({value: hexColor}, options)` (`column-helpers.ts:92`)
- `setCategorical(colorMap?, {fallbackColor, matchType})` (`column-helpers.ts:103`)
- `setConditional({'20-170': '#00FF00'})` (`column-helpers.ts:113`)
- `setDisabled()` (`column-helpers.ts:124`)
- `getType()`, `getColor(i)`, `getColors()` (`column-helpers.ts:55,128,132`)

Idiom from `samples/grid/color-coding/color-coding.js:4-13`:

```js
t.col('height').meta.colors.setConditional({'20-170': '#00FF00', '170-190': '#220505'});
t.col('age').meta.colors.setLinear(['#ff0000', '#ffff00', '#00ff00'], {min: 19, max: 70});
t.col('race').meta.colors.setCategorical({'Asian': 4278190335, 'Black': 4286578816});
```

**`col.meta.markers`** = `ColumnMarkerHelper`:
- `assign(category, marker)`, `default(marker)`, `setAll(map)`, `reset()` (`column-helpers.ts:137-175`)

**Constants Claude should always cite, not invent**:
- `DG.SEMTYPE` — string consts: `MOLECULE`, `MACROMOLECULE`, `MOLECULE3D`, `IC50`, `EC50`, `Ki`, `CONCENTRATION`, `VOLUME`, `EMAIL`, `URL`, `LATITUDE`, `LONGITUDE`, `IMAGE`, `FILE`, `CHEMICAL_REACTION`, ... (`js-api/const.ts:206-246`).
- `DG.COLUMN_TYPE` — `STRING='string'`, `INT='int'`, `FLOAT='double'` ⚠ (not `'float'`), `BOOL='bool'`, `BYTE_ARRAY='byte_array'`, `DATE_TIME='datetime'`, `BIG_INT='bigint'`, `QNUM='qnum'`, `DATA_FRAME='dataframe'`, `OBJECT='object'` (`const.ts:67-78`).
- `DG.TAGS` — leading-dot keys are hidden/system tags (e.g. `.color-coding-type`, `.choices`, `.tooltip`). Plain keys appear in details panel (`description`, `format`, `units`, `friendlyName`, `formula`, `quality`=SEMTYPE) (`const.ts:282-342`).
- `DG.UNITS.Molecule` — `SMILES`, `MOLBLOCK`, `V3K_MOLBLOCK`, `INCHI` (`const.ts:248-255`).
- `DG.COLOR_CODING_TYPE` — `CATEGORICAL`, `CONDITIONAL`, `LINEAR`, `OFF` (`const.ts:797-802`).

**Gotchas**:
- `DG.COLUMN_TYPE.FLOAT === 'double'` — Claude mistypes this constantly. Always prefer the constant over a literal.
- The semType tag key is actually `'quality'` not `'semType'` (`const.ts:318`). Use `col.semType` accessor, **never** `setTag('semType', ...)`.
- Setting `col.meta.choices = null` deletes the tag (`column-helpers.ts:244`).
- `col.tags` is a `MapProxy`; values must be **strings** — `setTag` throws otherwise on the DF (`data-frame.ts:179`); on `Column.setTag` no explicit guard, but Dart side will reject non-strings.

**Wrapper rationale**: `setColumnMeta(col, {semType, format, units, friendlyName, description, colorCoding})` collapses ~5 lines of mixed accessor calls into one call.

## 4. Stats (Area C)

`col.stats` returns `Stats` (`column.ts:434`). All metrics are cached on the column (`column.ts:421-424` for min/max also expose `col.min`/`col.max` as shorthands).

Surface (`js-api/dataframe/stats.ts:35-114`):
`totalCount`, `valueCount`, `missingValueCount`, `uniqueCount`, `min`, `max`, `sum`, `avg`, `stdev`, `variance`, `skew`, `kurt`, `med`, `q1`, `q2`, `q3`. Plus `corr(other)`, `spearmanCorr(other)` (`stats.ts:116-119`) and the static `histogramsByCategories` (`stats.ts:122`).

Sample (`samples/data-frame/stats.js:5-23`) builds a per-column stats table.

`col.categories` returns unique sorted strings for string columns (`column.ts:395`); `col.toList()` is documented as expensive (`column.ts:390`); `col.values()` is a generator (`column.ts:440`); `col.getRawData()` returns the underlying typed array — `Int32Array` for int and **string indexes**, `Float64Array` for floats/qnum/datetime, `Uint32Array` for bools (`column.ts:294-304`).

**Gotcha**: `col.stats.uniqueCount` is the cheap way to count distinct values; `col.toList().filter(…).length` is the wrong way.

## 5. Cloning DataFrames (Area D)

`df.clone(rowMask?, columnIds?, saveSelection=false, saveTags=true)` returns a new `DataFrame` (`data-frame.ts:290`).

- `df.clone()` — full copy.
- `df.clone(df.filter)` — only filtered rows.
- `df.clone(null, ['a', 'b'])` — only these columns.
- `df.clone(df.filter, ['a', 'b'], true)` — also preserves selection bits.

`col.clone(mask?)` returns a **detached** column (`column.ts:280`). The note explicitly says: "the cloned column is not added to this column's dataframe."

`df.copy()` — **does not exist**. Claude sometimes invents it. CSV/parquet/byte-array round-trips exist: `df.toCsv` (`data-frame.ts:241`), `df.toByteArray` / `DataFrame.fromByteArray` (`data-frame.ts:269` + `:75`), `df.toParquet`, `df.toArrow`.

**Wrapper rationale**: `cloneDf(df, {rows, cols, withSelection})` removes the positional arg confusion.

## 6. Adding columns (Area E)

All in `ColumnList`:
- `addNew(name, type)` (`column-list.ts:169`) — generic, requires `DG.COLUMN_TYPE.*`
- `addNewString` / `addNewInt` / `addNewFloat` / `addNewBool` / `addNewDateTime` / `addNewQnum` / `addNewBytes` (`column-list.ts:189-226`) — preferred when type is known.
- `addNewCalculated(name, formula, type='auto', treatAsString=false, subscribeOnChanges=true)` (`column-list.ts:180`) — **async, returns `Promise<Column>`**. Formula DSL uses `${col}` placeholders: `'${x}+${y}-${z}'` (sample `samples/data-frame/modification/calculated-columns/add-calculated-column.js:8`).
- `addNewVirtual(name, getValue, type=TYPE.OBJECT, setValue?)` (`column-list.ts:238`) — lazy; sample at `samples/data-frame/advanced/virtual-columns.js:23`.
- `add(col, notify?)` (`column-list.ts:143`) — insert prebuilt column at end.
- `insert(col, index?, notify?)` (`column-list.ts:163`) — at a specific position.
- `getOrCreate(name, type)` (`column-list.ts:152`) — idempotent add.
- `replace(old, new, notify?)` (`column-list.ts:263`) — preserves position; preserves nothing else.
- `getUnusedName(name, choices?)` (`column-list.ts:271`) — collision-free naming.

Init idioms (`samples/data-frame/modification/add-columns.js:15-20`):

```js
t.columns.addNewInt('int').init(5);              // scalar
t.columns.addNewFloat('float').init(i => i/10);  // function
```

`Column.fromList`, `fromStrings`, `fromInt32Array`/`fromFloat32Array`/`fromFloat64Array`/`fromBigInt64Array`, `fromIndexes`, `fromBitSet` plus typed factories `Column.int/float/string/bool/dateTime/qnum/dataFrame` exist (`column.ts:54-172`). Note `fromInt32Array` etc. **share** memory — the array is not copied (`column.ts:71-83`).

**Decision tree for Claude**:
1. Know type + want empty? → `addNewInt('x')` etc.
2. Need formula based on other columns? → `addNewCalculated` (`await`!).
3. Have raw data already? → `Column.fromList(type, name, list)` then `df.columns.add(col)`.
4. Need computed-on-demand cells (e.g. wraps JS objects)? → `addNewVirtual`.
5. Need to ensure unique name? → `getUnusedName` then `addNew`.

**Gotcha**: `addNewCalculated` is the only addNew that's `async`. Forgetting `await` causes downstream code to receive a `Promise`.

**Wrapper rationale**: `addColumn(df, {name, type?, formula?, values?, init?, semType?, format?, units?})` — chooses the right backing call and applies metadata in one go.

## 7. Removing / renaming (Area F)

`df.columns.remove(column, notify?)` accepts string, index, or `Column` (`column-list.ts:246`). Case-insensitive. Returns `this` (chainable).

Renaming: `col.name = 'newName'` (`column.ts:235`). Fires `onColumnNameChanged` (`data-frame.ts:451`). No batch helper exists.

**Gotchas**:
- Grid columns/viewer property bindings that reference by name **break silently** on rename. Layout matching uses `layoutColumnId` (`column.ts:222`) which is more stable.
- Removing a column that's actively shown in a viewer doesn't crash — viewer reacts to `onColumnsRemoved` (`data-frame.ts:463`).

**Wrapper rationale**: `removeColumns(df, names[])` for bulk remove; pass missing names to skip silently (or warn).

## 8. Setting values (Area G)

- `col.set(i, value, notify=true)` (`column.ts:354`). With `notify=false`, the DF event is suppressed until `col.fireValuesChanged()` (`column.ts:343`).
- `col.init(scalarOrFn)` (`column.ts:264`) — clean bulk init. Function arg gets row index.
- `col.setString(i, str, notify?)` (`column.ts:338`) — parses string into the column's type, returns success.
- `col.getRawData()` / `col.setRawData(arr, notify?)` (`column.ts:294`/`311`) — fastest path, raw typed array.
- `col.values()` generator (`column.ts:440`).

For bulk writes, use either `init` or write directly into `getRawData()` and then call `setRawData` + `fireValuesChanged`. Looping `col.set(i, v, true)` for every row fires a notification per write, which is slow.

Sample for init: `samples/data-frame/modification/init-column.js:11-15`.

## 9. Rows (Area H)

- `df.rows.removeAt(idx, count=1, notify?)` (`row.ts:226`)
- `df.rows.removeWhere(row => …)` (`row.ts:232`)
- `df.rows.removeWhereIdx(i => …)` (`row.ts:238`)
- `df.rows.insertAt(idx, count?, notify?)` (`row.ts:245`)
- `df.rows.addNew(values?, notify?)` (`row.ts:253`)
- `df.rows.setValues(idx, values, notify?)` (`row.ts:267`)
- `df.rows.match(query|object)` (`row.ts:274`), `df.rows.select(pred)` (`row.ts:289`), `df.rows.filter(pred)` (`row.ts:297`)
- `df.row(i).get(name)` / `df.row(i).set(name, value)` via proxy (`row.ts:60-117`)

Sample (`samples/data-frame/modification/remove-where.js:7`):
```js
t.rows.removeWhere(row => row.get('make') === 'Honda');
```

**Anti-pattern, explicit in source** (`row.ts:191-192`):
> Refrain from accessing data via `RowList` and `Row` in performance-critical scenarios. To maximize performance, get values via `DataFrame.columns`, instead.

So `for (const row of df.rows) { row.foo = bar(row.bar) }` is the wrong shape; the right shape is

```js
const src = df.getCol('bar');
const dst = df.getCol('foo');
for (let i = 0; i < df.rowCount; i++)
  dst.set(i, bar(src.get(i)), false);
dst.fireValuesChanged();
```

## 10. Current state on DataFrame (Area I)

- `df.currentRowIdx` get/set (`data-frame.ts:300`)
- `df.currentRow` get/set (Row wrapper) (`data-frame.ts:296`)
- `df.currentCol` get/set (`data-frame.ts:309-317`)
- `df.currentCell` (`data-frame.ts:321`)
- `df.mouseOverRowIdx` (`data-frame.ts:304`)
- `df.selection` — `BitSet` of selected rows (`data-frame.ts:161`)
- `df.filter` — `BitSet` of rows that pass filters (`data-frame.ts:51`, constructed in ctor `:62`)

Events: `onCurrentRowChanged`, `onCurrentColChanged`, `onCurrentCellChanged`, `onSelectionChanged`, `onFilterChanged` (`data-frame.ts:427-493`).

Full deep-dive belongs in a separate "state" skill. For this skill, surface only:
- "current row" = `df.currentRowIdx` (not `0`).
- "selected rows" = iterate `df.selection.getSelectedIndexes()` (BitSet) or use `df.rows.indexes({onlySelected: true})` (`row.ts:218`).
- "filtered rows" = `df.filter` (BitSet) — iterate via `df.rows.indexes({onlyFiltered: true})`.

## 11. Existing grokky helpers (Area J)

From `Grokky/src/claude/grokky-api.ts`:

| Helper | Lines | Relevance to this skill |
|--------|-------|-------------------------|
| `addCalculatedColumn(df, name, formula, type?)` | `:4-11` | Wraps `addNewCalculated` + perf log. Already covered. |
| `filterRows(view, column, criteria)` | `:20-31` | Belongs in a state/filter skill, not here. |
| `sortRows(df, columns, orders)` | `:33-36` | Adjacent, not central. |
| `selectRows(df, predicate)` | `:38-48` | State skill. |
| `aggregateBy(df, groupCols, aggMap)` | `:50-56` | Stats/aggregation skill (out of scope). |
| `pivot(df, opts)` / `unpivot(df, opts)` | `:58-69` | Adjacent. |
| `joinTables(left, right, opts)` | `:78-88` | Adjacent. |
| `addViewer` / `configureViewer` | `:162-173` | Viewer skill. |

**Nothing in `grokky.*` covers**: column finding by intent, semantic-type discovery, column-metadata bulk set, df cloning, column removal, descriptive snapshot, color coding. All of §3, §4, §5, §6, §7 are unwrapped.

## 12. Proposed grokky helpers

```ts
/** Find columns by free-form intent. Tries (in order): exact name (case-insensitive),
 *  semType match, type filter, tag predicate, fuzzy name match. Returns ranked candidates
 *  with a short "why" so the caller can disambiguate. */
function findColumns(df: DG.DataFrame, query: {
  name?: string;
  semType?: DG.SemType | DG.SemType[];
  type?: DG.ColumnType | 'numerical' | 'categorical' | 'dateTime';
  tag?: Record<string, string | null | undefined>;
  fuzzy?: boolean;       // default true
  limit?: number;        // default 5
}): Array<{ column: DG.Column; score: number; reason: string }>;

/** Concise JSON-friendly snapshot of a column: type, semType, units, format, length,
 *  stats (min/max/avg/missing/unique), top 5 categories for string cols. Cheap; uses cached stats. */
function describeColumn(col: DG.Column): {
  name: string; type: DG.ColumnType; semType: string | null;
  length: number; missing: number; unique: number;
  units?: string; format?: string; friendlyName?: string;
  numerical?: { min: number; max: number; avg: number; stdev: number };
  categorical?: { topCategories: string[] };
};

/** Bulk-set column metadata. `null` clears the tag (matches `meta.setNonNullTag` semantics).
 *  Color coding is handled here instead of forcing callers to know the helper path. */
function setColumnMeta(col: DG.Column, meta: {
  semType?: string | null;
  format?: string | null;
  units?: string | null;
  friendlyName?: string | null;
  description?: string | null;
  choices?: string[] | null;
  colorCoding?:
    | { kind: 'linear'; range?: number[]; min?: number; max?: number; belowMinColor?: string; aboveMaxColor?: string }
    | { kind: 'categorical'; map?: Record<string, string | number>; fallbackColor?: string | number }
    | { kind: 'conditional'; rules: Record<string, string | number> }
    | { kind: 'off' };
  tags?: Record<string, string | null>;  // raw tag escape hatch
}): DG.Column;

/** Pure rewrite of df.clone with named args. Returns a new DataFrame. */
function cloneDf(df: DG.DataFrame, opts?: {
  rows?: DG.BitSet | 'filtered' | 'selected';
  cols?: string[];
  withSelection?: boolean;
  withTags?: boolean;            // default true
}): DG.DataFrame;

/** Add a column choosing the right backend call. Supports raw values, init function,
 *  formula (async path), or empty + typed. Applies metadata if provided. */
async function addColumn(df: DG.DataFrame, spec: {
  name: string;
  type?: DG.ColumnType;
  formula?: string;              // mutually exclusive with values/init
  values?: any[];                // mutually exclusive with formula/init
  init?: (i: number) => any;     // mutually exclusive with formula/values
  insertAt?: number;
  ensureUniqueName?: boolean;    // default true — uses getUnusedName
  meta?: Parameters<typeof setColumnMeta>[1];
}): Promise<DG.Column>;

/** Bulk remove. Missing names: 'skip' (default) | 'throw' | 'warn'. */
function removeColumns(df: DG.DataFrame, names: string[], onMissing?: 'skip' | 'throw' | 'warn'): string[];

/** Helper sugar: rename with collision-check using getUnusedName. */
function renameColumn(df: DG.DataFrame, from: string, to: string, opts?: { ensureUnique?: boolean }): string;

/** One-line: top-N categories with counts. */
function topCategories(col: DG.Column, n?: number): Array<{ value: string; count: number }>;
```

### Rationale per helper

- `findColumns` — directly addresses the audit gap "picks wrong column when several plausible exist". The `reason` field is the key piece: it lets Claude show the user *why* a column was chosen and offer alternates.
- `describeColumn` — single call returns everything Claude needs to talk about a column. Replaces the ad-hoc 5-10 line pattern Claude generates today.
- `setColumnMeta` — fixes the "uses setTag for color coding" anti-pattern by making the correct path the *easy* path. Also clarifies tag-vs-meta distinction.
- `cloneDf` — kills positional-args confusion (`df.clone(null, null, false, true)`).
- `addColumn` — kills the "pick the right `addNew*`" decision tree; one function dispatches based on which of `formula/values/init` is provided.
- `removeColumns` — the missing-name policy is the real value; users say "remove the temp columns" and Claude doesn't have to guard.
- `renameColumn` — minor sugar, mostly the unique-name guard.
- `topCategories` — cheap, very common ask, currently requires `col.toList()` + `reduce`.

## 13. Anti-patterns Claude should be warned against

1. **Looping over `df.rows` for column work.** Source comment: `row.ts:191-192`. Use `df.getCol(name)` + index loop, with `notify=false` and a final `fireValuesChanged`.
2. **`setTag('color-coding-…', …)` directly.** Use `col.meta.colors.setLinear/setCategorical/setConditional/setDisabled`. See `column-helpers.ts:48-135`.
3. **`setTag('semType', 'Molecule')`.** Use `col.semType = DG.SEMTYPE.MOLECULE`. The tag key is actually `'quality'` (`const.ts:318`).
4. **Hard-coding `'double'` / `'float'` strings.** `DG.COLUMN_TYPE.FLOAT === 'double'` — always use the constant.
5. **`df.copy()`.** Doesn't exist. Use `df.clone()`.
6. **`col.toList().length` for value count.** Use `col.length` (total) or `col.stats.valueCount` (non-null).
7. **`new Set(col.toList()).size` for unique count.** Use `col.stats.uniqueCount`. For string cols, also `col.categories.length`.
8. **Mutating raw arrays without notify.** `setRawData(arr, notify=true)` works; if you bypass it via `getRawData()` direct write, call `col.fireValuesChanged()`.
9. **Forgetting `await` on `addNewCalculated`.** It's the only `addNew*` that returns a Promise (`column-list.ts:180`).
10. **`df.columns.byName(...)` followed by `if (col == null)` without `getCol`.** When you *want* to throw on missing, use `df.getCol(name)` (`data-frame.ts:230`).
11. **Renaming a column expected by a saved layout.** Use `col.layoutColumnId` for stable identity (`column.ts:222`).
12. **`fromInt32Array`/`fromFloat32Array` assuming the array is copied.** It is **not** — `column.ts:71-83` says explicitly "will be used as column's storage". Mutations to the source array show up in the column.

## 14. Open questions

1. **Skill scope on state**: Area I (current row / selection / filter) is mentioned here but earmarked for a separate skill. Confirm whether `findColumns` should also support `{onlyFiltered: true}` or whether that lives in the state skill.
2. **Fuzzy matching budget**: `findColumns` needs a similarity function — should it reuse the Levenshtein already in `grokky-api.ts:100-114`, or pull in a token-based matcher (Jaccard on word tokens) for multi-word names like "Mol Weight"?
3. **`describeColumn` cost**: `col.stats.uniqueCount` for very wide string columns can be expensive. Cap at, say, 100k rows? Or always run — Datagrok caches the result on the column anyway.
4. **`setColumnMeta` and color coding off**: Should passing `colorCoding: undefined` leave existing coding alone, while `colorCoding: { kind: 'off' }` clears it? Proposed: yes.
5. **`addColumn` and async**: Should we always return `Promise<DG.Column>` for API uniformity (even synchronous paths), or only when `formula` is used? Lean toward always async for predictability.
6. **`renameColumn` and viewers**: Should we attempt to update viewer property bindings that reference the old name (very fragile), or just warn?
7. **Bulk tag clearing**: The skill should likely teach "clear all custom tags" as an idiom. The current API has `col.tags.clear()` (`MapProxy` supports it per `samples/data-frame/metadata/metadata.js:76`), but that also nukes system-set tags like `.choices`. Worth a `clearUserTags(col)` helper? Marked low priority pending feedback.
