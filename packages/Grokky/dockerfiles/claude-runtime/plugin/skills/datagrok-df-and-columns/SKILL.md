---
name: datagrok-df-and-columns
description: Find, describe, add, remove, rename, clone, or set metadata on columns of a Datagrok DataFrame inside a datagrok-exec block. Use whenever the user asks to locate "the X column", summarize a column, add a typed/empty/values-filled/virtual column, set semantic type / units / format / friendly name, apply linear or categorical or conditional color coding, drop or rename columns, or copy a DataFrame. Covers everything in DataFrame.columns and Column.meta — but not row filtering/selection (datagrok-table-ops) and not formula-only columns (datagrok-calc-column).
---

# datagrok-df-and-columns

Use the `grokky.*` column helpers inside a `datagrok-exec` block. They wrap the
js-api so Claude doesn't have to remember whether a column is added via
`addNewFloat` vs `addNewCalculated` vs `Column.fromList`, that the float type
constant is the string `'double'`, or that the semType tag key is `quality`
under the hood.

If the user wants a **formula-driven** column (recomputes on source change),
use the `datagrok-calc-column` skill — `grokky.addCalculatedColumn` is the
right entry point there. This skill handles every *other* column operation.

## What this skill covers

DataFrame-level: locating columns by intent, cloning a DataFrame, removing
columns in bulk. Column-level: describing a column (type, stats, top
categories), adding columns (typed-empty, values, init function, virtual),
setting metadata (semType, format, units, friendly name, color coding) and
renaming. Cell-value mutation is touched briefly — `col.init(fn)` for bulk,
plus a perf note. Current row / selection / filter live in a separate state
skill.

## Quick reference

| Helper                                  | One-liner                                                |
|-----------------------------------------|----------------------------------------------------------|
| `grokky.findColumns(df, query)`         | Ranked column candidates by name/semType/type/tag/fuzzy. |
| `grokky.describeColumn(col)`            | JSON snapshot: type, stats, top categories.              |
| `grokky.cloneDf(df, {rows, cols, ...})` | Named-args wrapper over `df.clone(...)`.                 |
| `grokky.addColumn(df, spec)`            | Dispatches to typed / values / init / formula / virtual. |
| `grokky.removeColumns(df, names, ...)`  | Bulk remove with missing-name policy.                    |
| `grokky.renameColumn(df, from, to)`     | Rename + optional unique-name guard. Warns about viewers.|
| `grokky.setColumnMeta(col, meta)`       | Bulk metadata set, including color coding.               |
| `grokky.topCategories(col, n)`          | Top-N categories with counts.                            |

Globals available inside every `datagrok-exec` block: `grok`, `ui`, `DG`,
`view`, `t` (the current `DG.DataFrame`, when the view is a TableView),
`grokky`.

## Finding the right column

When the user says "the molecule column" or "the activity values", pick the
narrowest tool that works:

| User said...                              | Right call                                              |
|-------------------------------------------|---------------------------------------------------------|
| "column named X" (must exist)             | `t.getCol('X')` — throws if missing                     |
| "column named X" (may be absent)          | `t.col('X')` — returns `null` if missing                |
| "the molecule column"                     | `t.columns.bySemType(DG.SEMTYPE.MOLECULE)`              |
| "every molecule column"                   | `t.columns.bySemTypeAll(DG.SEMTYPE.MOLECULE)`           |
| "all numeric columns"                     | iterate `t.columns.numerical`                           |
| "the mass column" (vague name, no semType)| `grokky.findColumns(t, {name: 'mass', fuzzy: true})`    |
| "the activity column" (multiple hits)     | `grokky.findColumns(t, {...})` — ranked + reason field  |
| "an IC50 column"                          | `grokky.findColumns(t, {semType: DG.SEMTYPE.IC50})`     |
| "columns tagged `quality=Molecule`"       | `t.columns.byTags({quality: 'Molecule'})`               |

`t.col` and `t.columns.byName` are **case-insensitive**. Iteration order
follows column order in the DataFrame.

`findColumns` is column-discovery only — it does not look at row filters or
selection. It returns at most `limit` (default 5) candidates, each with a
`score` (0–1) and a short `reason` string ("exact name", "semType match",
"fuzzy name distance 2", "tag match `units=kg`") so Claude can disambiguate
in the response.

If the query is empty (all fields `undefined`), `findColumns` returns `[]` —
it deliberately does not dump the whole column list.

```datagrok-exec
// Find the molecule column. If there are several, prefer the one the user
// most likely meant — molecules with structure are SMILES/MOLBLOCK semType.
const hits = grokky.findColumns(t, {semType: DG.SEMTYPE.MOLECULE});
return hits;
```

```datagrok-exec
// User said "the mw column" but the DF has "Molecular Weight". Fall back to
// fuzzy match. Score and reason let Claude phrase the response well.
const hits = grokky.findColumns(t, {name: 'mw', fuzzy: true, limit: 3});
return hits;
```

## Describing a column

`grokky.describeColumn(col)` returns a JSON-friendly snapshot — type, semType,
units, format, length, missing count, unique count, plus `numerical` stats
(min/max/avg/stdev) for numeric columns and `topCategories` for string
columns. Cheap: it relies on `col.stats` which DG caches per column.

```datagrok-exec
const col = grokky.findColumns(t, {name: 'activity', fuzzy: true})[0]?.column ?? t.columns.numerical[0];
return grokky.describeColumn(col);
```

To describe many columns at once, map `grokky.describeColumn` over a column
selector. Always use the helper rather than inlining `col.stats` access — the
helper handles missing values, distinguishes numeric vs categorical output,
and keeps response shapes consistent.

```datagrok-exec
// Summarize every numeric column
return t.columns.numerical.map(grokky.describeColumn);
```

```datagrok-exec
// Summarize every molecule column
return t.columns.bySemTypeAll(DG.SEMTYPE.MOLECULE).map(grokky.describeColumn);
```

For raw stats access, `col.stats` exposes `totalCount`, `valueCount`,
`missingValueCount`, `uniqueCount`, `min`, `max`, `sum`, `avg`, `stdev`,
`variance`, `med`, `q1`, `q2`, `q3`, `skew`, `kurt`, `corr(other)`,
`spearmanCorr(other)`. `col.min` / `col.max` are shorthands.

| Need                              | Use                                      |
|-----------------------------------|------------------------------------------|
| total row count                   | `col.length`                             |
| non-null value count              | `col.stats.valueCount`                   |
| missing-value count               | `col.stats.missingValueCount`            |
| distinct value count              | `col.stats.uniqueCount`                  |
| string-column sorted distinct     | `col.categories`                         |
| typed-array raw view              | `col.getRawData()`                       |
| value iterator                    | `col.values()` (generator)               |
| top-N category counts             | `grokky.topCategories(col, n)`           |

Do **not** compute these by hand:

- `col.length` for total rows (not `col.toList().length`).
- `col.stats.valueCount` for non-null count.
- `col.stats.uniqueCount` for distinct values (not `new Set(col.toList()).size`).
- For string columns, `col.categories` is the sorted unique-string array.

## Cloning DataFrames

`df.copy()` does **not** exist. Use `df.clone(...)` or `grokky.cloneDf(...)`.

`df.clone(rowMask?, columnIds?, saveSelection=false, saveTags=true)` is the
raw API — positional args, easy to get wrong. `grokky.cloneDf` accepts named
args and supports two ergonomic strings for `rows`: `'filtered'` (uses
`df.filter`) and `'selected'` (uses `df.selection`).

```datagrok-exec
// Full copy
const copy = grokky.cloneDf(t);
return copy;
```

```datagrok-exec
// Just two columns, filtered rows only
const small = grokky.cloneDf(t, {cols: ['smiles', 'activity'], rows: 'filtered'});
return small;
```

```datagrok-exec
// Raw df.clone is fine for one specific case — no named args needed
const justSelected = t.clone(t.selection);
return justSelected;
```

`col.clone(mask?)` exists too — note that **the cloned column is not attached
to a DataFrame**. You need `df.columns.add(clonedCol)` if you want it back in.

## Adding columns

Pick the right shape first; the helper just dispatches.

| You have...                            | Right path                                              |
|----------------------------------------|---------------------------------------------------------|
| nothing — want empty typed column      | `grokky.addColumn(df, {name, type})`                    |
| an array of values                     | `grokky.addColumn(df, {name, values: [...]})`           |
| a function of row index                | `grokky.addColumn(df, {name, type, init: i => ...})`    |
| a formula using other columns          | use `datagrok-calc-column` skill (`addCalculatedColumn`)|
| compute-on-demand (no storage)         | `grokky.addColumn(df, {name, virtual: i => ..., type})` |

`grokky.addColumn` is **always async** for return-type uniformity, even when
the underlying path is synchronous. Always `await` it.

Rules enforced by the helper:

- Exactly zero or one of `formula` / `values` / `init` / `virtual` may be set.
- `name` is run through `df.columns.getUnusedName(name)` when
  `ensureUniqueName !== false` (default true).
- If `meta` is supplied, it's applied via `setColumnMeta` before returning.

```datagrok-exec
// Typed empty column, then apply metadata in the same call
const col = await grokky.addColumn(t, {
  name: 'Ki',
  type: DG.COLUMN_TYPE.FLOAT,
  meta: {semType: DG.SEMTYPE.Ki, units: 'nM', format: '0.00'},
});
return col;
```

```datagrok-exec
// Formula-driven column with recompute. Prefer datagrok-calc-column, but
// addColumn delegates correctly when {formula} is set.
const lipE = await grokky.addColumn(t, {
  name: 'LipE',
  formula: '${pIC50} - ${cLogP}',
  type: DG.COLUMN_TYPE.FLOAT,
});
return lipE;
```

```datagrok-exec
// Virtual column — computed on each access, not stored. Useful for huge tables.
const labels = await grokky.addColumn(t, {
  name: 'Label',
  virtual: (i) => `${t.get('compound', i)}@${t.get('target', i)}`,
  type: DG.COLUMN_TYPE.STRING,
});
return labels;
```

Underneath, the dispatch is:

- `formula` → `df.columns.addNewCalculated(name, formula, type ?? 'auto')` (async).
- `values` → `DG.Column.fromList(type ?? 'string', name, values)` then `df.columns.add(...)`.
- `init` → typed `addNew*` then `col.init(fn)`.
- `virtual` → `df.columns.addNewVirtual(name, fn, type)`.
- none of the above → typed `addNew*` (empty).

## Removing & renaming

`grokky.removeColumns(df, names, onMissing='skip'|'throw'|'warn')` accepts a
list of names and returns the list of names that were actually removed.
Default is `'skip'`: missing names are silently ignored, which matches the
"remove the temp columns we added" intent.

```datagrok-exec
// Drop bulk by name, skipping any that don't exist
const removed = grokky.removeColumns(t, ['_tmp1', '_tmp2', 'doesNotExist']);
return removed;  // ['_tmp1', '_tmp2']
```

```datagrok-exec
// Strict mode — throw if any are missing
grokky.removeColumns(t, ['mustExist'], 'throw');
```

To target **calculated columns** specifically (e.g. "remove all the calculated
columns we added"), use `col.isCalculated` — a first-class accessor on
`DG.Column` that's true when the column has a formula. Combine with
`removeColumns`:

```datagrok-exec
// Drop every calculated column
const names = [...t.columns].filter(c => c.isCalculated).map(c => c.name);
return grokky.removeColumns(t, names);
```

The skill does **not** ship a session log — there's no built-in way to know
which columns were added "this session." `col.isCalculated` (or, equivalently,
checking `col.meta.formula`) is the right heuristic. The accessor is faster
and reads cleaner.

Rename: just assign `col.name = 'newName'`, or use `grokky.renameColumn(df,
from, to, {ensureUnique: true})` for a collision-safe version. The helper
returns the final name actually used (may differ if `ensureUnique` had to
disambiguate).

```datagrok-exec
const finalName = grokky.renameColumn(t, 'pIC50', 'pIC50_old', {ensureUnique: true});
return finalName;
```

Rename caveat: **viewer property bindings referencing a column by name break
silently**. Layouts use `col.layoutColumnId` for stable identity instead.
`grokky.renameColumn` issues a `console.warn` if the column appears bound to
any viewer in the current view — it does **not** attempt to rewrite viewer
bindings (too fragile).

## Setting column metadata

Three things live under "column metadata" and they all have different
accessors. Cheat sheet:

| Where metadata lives                     | Read                                       | Write                                                |
|------------------------------------------|--------------------------------------------|------------------------------------------------------|
| Semantic type                            | `col.semType`                              | `col.semType = DG.SEMTYPE.MOLECULE`                  |
| Friendly name / description / format / units / formula / choices / cellRenderer | `col.meta.friendlyName` etc.        | `col.meta.units = 'nM'`; assign `null` to clear      |
| Color coding                             | `col.meta.colors.getType()`                | `col.meta.colors.setLinear/setCategorical/setConditional/setDisabled` |
| Arbitrary tags                           | `col.getTag(k)`, `col.tags.has(k)`         | `col.setTag(k, v)`, `col.tags.delete(k)`             |
| Layout-stable identity                   | `col.layoutColumnId`                       | `col.layoutColumnId = '...'`                         |

`grokky.setColumnMeta(col, meta)` collapses all of the above into one call.
Fields:

- `semType` — `col.semType = ...` (uses the right accessor, never `setTag('semType', ...)`).
- `format`, `units`, `friendlyName`, `description`, `choices` — go via `col.meta.*`.
- `colorCoding` — discriminated union; see below.
- `tags` — raw escape hatch, sets `col.setTag(k, v)` or `col.tags.delete(k)` if `v === null`.

`null` clears a metadata field (removes the underlying tag). `undefined`
**leaves it alone** — this is the key distinction. Pass `colorCoding:
undefined` to leave existing coding untouched; pass `colorCoding: {kind:
'off'}` to actively clear it.

Color coding shapes:

| `colorCoding` value                                                            | Effect                                  |
|--------------------------------------------------------------------------------|-----------------------------------------|
| `{kind: 'linear', range, min, max, belowMinColor, aboveMaxColor}`              | Numerical only — throws on string col.  |
| `{kind: 'categorical', map, fallbackColor}`                                    | Per-category colors.                    |
| `{kind: 'conditional', rules: {'20-170': '#00FF00'}}`                          | Range-rule based.                       |
| `{kind: 'off'}`                                                                | Disable color coding.                   |

```datagrok-exec
// Linear color coding on a numeric column
grokky.setColumnMeta(t.getCol('age'), {
  colorCoding: {kind: 'linear', range: ['#ff0000', '#ffff00', '#00ff00'], min: 19, max: 70},
});
```

```datagrok-exec
// Categorical color coding
grokky.setColumnMeta(t.getCol('race'), {
  colorCoding: {kind: 'categorical', map: {'Asian': '#0000FF', 'Black': '#FF0000'}},
});
```

```datagrok-exec
// Conditional color coding by range
grokky.setColumnMeta(t.getCol('height'), {
  colorCoding: {kind: 'conditional', rules: {'20-170': '#00FF00', '170-190': '#220505'}},
});
```

```datagrok-exec
// Set semType + units + format in one shot, also clear any old description
grokky.setColumnMeta(t.getCol('Ki'), {
  semType: DG.SEMTYPE.Ki,
  units: 'nM',
  format: '0.00',
  description: null,
});
```

```datagrok-exec
// Turn off color coding
grokky.setColumnMeta(t.getCol('activity'), {colorCoding: {kind: 'off'}});
```

## Working with cell values

For bulk init, `col.init(scalar | (i) => value)` is the right shape:

```datagrok-exec
// Bulk init via function — fastest correct path
const col = await grokky.addColumn(t, {name: 'rowSquared', type: DG.COLUMN_TYPE.INT});
col.init((i) => i * i);
return col;
```

Avoid looping `col.set(i, v, true)` for every row — each call notifies. For
high-volume mutation, either use `col.init(fn)`, or write values with
`notify=false` and call `col.fireValuesChanged()` once at the end:

```datagrok-exec
const src = t.getCol('input');
const dst = t.getCol('output');
for (let i = 0; i < t.rowCount; i++)
  dst.set(i, src.get(i) * 2, false);
dst.fireValuesChanged();
```

For truly hot loops, `col.getRawData()` returns the underlying typed array
(`Int32Array`/`Float64Array`/`Uint32Array`). Mutate it directly, then call
`col.fireValuesChanged()`. Note: `Column.fromInt32Array` and siblings **share
memory** with the source array — mutating the array mutates the column.

## Constants Claude must always cite, never invent

| Constant family             | Members                                                              |
|-----------------------------|----------------------------------------------------------------------|
| `DG.COLUMN_TYPE`            | `STRING='string'`, `INT='int'`, **`FLOAT='double'`** (!), `BOOL='bool'`, `BYTE_ARRAY='byte_array'`, `DATE_TIME='datetime'`, `BIG_INT='bigint'`, `QNUM='qnum'`, `DATA_FRAME='dataframe'`, `OBJECT='object'` |
| `DG.SEMTYPE`                | `MOLECULE`, `MACROMOLECULE`, `MOLECULE3D`, `IC50`, `EC50`, `Ki`, `CONCENTRATION`, `VOLUME`, `EMAIL`, `URL`, `LATITUDE`, `LONGITUDE`, `IMAGE`, `FILE`, `CHEMICAL_REACTION`, ... |
| `DG.UNITS.Molecule`         | `SMILES`, `MOLBLOCK`, `V3K_MOLBLOCK`, `INCHI`                        |
| `DG.COLOR_CODING_TYPE`      | `CATEGORICAL`, `CONDITIONAL`, `LINEAR`, `OFF`                        |
| `DG.TAGS`                   | tag keys — leading-dot for system tags (`.color-coding-type`, `.choices`). The semType tag is actually `quality` — never `setTag('semType', ...)` |

**The single biggest footgun is `DG.COLUMN_TYPE.FLOAT === 'double'`.** Never
write the literal `'float'` or `'double'` in code — always use the constant.

```datagrok-exec
// CORRECT
await grokky.addColumn(t, {name: 'x', type: DG.COLUMN_TYPE.FLOAT});

// WRONG — there is no 'float' type
// t.columns.addNew('x', 'float');
```

## Anti-patterns — and the fix

1. **Looping `df.rows` to read/write column data.** Wrong shape. Use
   `df.getCol(name)` and an index loop with `notify=false` + a final
   `col.fireValuesChanged()`. The source code explicitly calls this out
   (`row.ts:191-192`).
2. **`col.setTag('color-coding-type', ...)`.** Use
   `col.meta.colors.setLinear/setCategorical/setConditional/setDisabled` — or
   `grokky.setColumnMeta({colorCoding: ...})`.
3. **`col.setTag('semType', 'Molecule')`.** The tag key is `quality`, not
   `semType`. Use `col.semType = DG.SEMTYPE.MOLECULE` or
   `grokky.setColumnMeta({semType: DG.SEMTYPE.MOLECULE})`.
4. **`'float'` or `'double'` string literals.** Always `DG.COLUMN_TYPE.FLOAT`.
5. **`df.copy()`.** Doesn't exist. Use `df.clone()` or `grokky.cloneDf(...)`.
6. **`col.toList().length` for row count.** Use `col.length`, or
   `col.stats.valueCount` for non-null count.
7. **`new Set(col.toList()).size` for unique count.** Use
   `col.stats.uniqueCount`. For string cols, `col.categories.length` also works.
8. **Mutating `col.getRawData()` and forgetting to fire.** Always end with
   `col.fireValuesChanged()` — otherwise viewers/grid don't refresh.
9. **Forgetting `await` on `addNewCalculated` / `grokky.addColumn`.** Both
   return Promises. Missing the `await` gives you a `Promise` object as the
   "column".
10. **`df.columns.byName(name)` plus `if (col == null)` when you want an
    error.** Use `df.getCol(name)` — it throws on missing.
11. **Renaming a column that a saved layout references by name.** Layouts
    track `col.layoutColumnId` for stable identity; renames break only the
    in-view bindings, which our helper warns about.
12. **`Column.fromInt32Array(arr)` assuming `arr` is copied.** It is not. The
    column shares memory with the array. Either `slice()` the input or
    treat it as the column's storage.

## State

This skill is column-shape, not row-shape. Current row (`df.currentRowIdx`),
selection (`df.selection`), and filter (`df.filter`) belong to the separate
state skill. `findColumns` deliberately ignores them. For row work, see
`datagrok-table-ops` (filter / sort / select) and `datagrok-calc-column`
(formula columns).
