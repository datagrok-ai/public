---
name: datagrok-df-and-columns
description: Find, describe, add, remove, rename, clone, or set metadata on columns of a Datagrok DataFrame inside a datagrok-exec block. Use whenever the user asks to locate "the X column", summarize a column, add a typed/empty/values-filled/virtual column, set semantic type / units / format / friendly name, apply linear or categorical or conditional color coding, drop or rename columns, or copy a DataFrame. Covers everything in DataFrame.columns and Column.meta — but not row filtering/selection (datagrok-filtering, datagrok-selection) and not formula-only columns (datagrok-calc-column).
---

# datagrok-df-and-columns

Use the `DG.DataFrame` / `DG.Column` js-api inside a `datagrok-exec` block.
Globals available in the block: `grok`, `ui`, `DG`, `view`, `t` (the current
`DG.DataFrame`, when `view.type === 'TableView'`).

For **formula-driven** columns (recompute on source change), use
`datagrok-calc-column`. Row filtering / selection lives in `datagrok-filtering`
and `datagrok-selection`.

## Quick reference

| Need                                         | Call                                                            |
|----------------------------------------------|-----------------------------------------------------------------|
| Get column by name (must exist)              | `t.getCol(name)` — throws on missing                            |
| Get column by name (may be absent)           | `t.col(name)` — returns `null`                                  |
| First column with a given semType            | `t.columns.bySemType(DG.SEMTYPE.MOLECULE)`                      |
| Every column with a given semType            | `t.columns.bySemTypeAll(DG.SEMTYPE.MOLECULE)`                   |
| All numeric / categorical / dateTime columns | `t.columns.numerical` / `.categorical` / `.dateTime`            |
| All names                                    | `t.columns.names()`                                             |
| Add typed empty column                       | `t.columns.addNewFloat(name)` (and `addNewInt/String/Bool/...`) |
| Add column from values                       | `t.columns.add(DG.Column.fromList(type, name, values))`         |
| Virtual column (compute-on-demand)           | `t.columns.addNewVirtual(name, i => ..., type)`                 |
| Collision-free name                          | `t.columns.getUnusedName(name)`                                 |
| Remove a column                              | `t.columns.remove(name)` (string, index, or `Column`)           |
| Rename                                       | `col.name = 'new'`                                              |
| Clone DataFrame (full)                       | `t.clone()`                                                     |
| Clone filtered rows                          | `t.clone(t.filter)`                                             |
| Clone column subset                          | `t.clone(null, ['a', 'b'])`                                     |
| Stats (cached on column)                     | `col.stats.{min,max,avg,stdev,med,q1,q2,q3,sum,valueCount,missingValueCount,uniqueCount}` |
| Set semType                                  | `col.semType = DG.SEMTYPE.MOLECULE`                             |
| Set friendly name / units / format / desc    | `col.meta.friendlyName / .units / .format / .description = ...` |
| Linear color coding (numeric)                | `col.meta.colors.setLinear(range, opts?)`                       |
| Categorical color coding                     | `col.meta.colors.setCategorical(map, opts?)`                    |
| Conditional color coding                     | `col.meta.colors.setConditional(rules)`                         |
| Disable color coding                         | `col.meta.colors.setDisabled()`                                 |

## Finding the right column

| User said...                                | Right call                                              |
|---------------------------------------------|---------------------------------------------------------|
| "column named X" (must exist)               | `t.getCol('X')` — throws if missing                     |
| "column named X" (may be absent)            | `t.col('X')` — returns `null` if missing                |
| "the molecule column"                       | `t.columns.bySemType(DG.SEMTYPE.MOLECULE)`              |
| "every molecule column"                     | `t.columns.bySemTypeAll(DG.SEMTYPE.MOLECULE)`           |
| "all numeric columns"                       | iterate `t.columns.numerical`                           |
| "an IC50 column"                            | `t.columns.bySemType(DG.SEMTYPE.IC50)`                  |
| "columns tagged `quality=Molecule`"         | `t.columns.byTags({quality: 'Molecule'})`               |
| "any of these columns: X, Y, Z (first hit)" | `t.columns.firstWhere(c => names.includes(c.name))`     |
| "does the DF have a column named X?"        | `t.columns.contains('X')` (case-insensitive)            |

`t.col` and `t.columns.byName` are case-insensitive. Always pass the semType
*constant* (`DG.SEMTYPE.MOLECULE`), never the string literal `'Molecule'`.

```datagrok-exec
// Find the molecule column. bySemType returns the first match or null.
const molCol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);
return ui.divText(molCol ? `Molecule column: ${molCol.name}` : 'No molecule column');
```

```datagrok-exec
// User said "mw" but the DF has "Molecular Weight" — try exact first, then
// fall back to a name substring match.
const wanted = 'mw';
let col = t.col(wanted);
if (!col) {
  const lc = wanted.toLowerCase();
  col = t.columns.firstWhere((c) => c.name.toLowerCase().includes(lc));
}
return ui.divText(col ? col.name : 'no match');
```

## Describing a column

`col.stats` is cached per column — repeated access is cheap.

| Need                              | Use                                      |
|-----------------------------------|------------------------------------------|
| total row count                   | `col.length`                             |
| non-null value count              | `col.stats.valueCount`                   |
| missing-value count               | `col.stats.missingValueCount`            |
| distinct value count              | `col.stats.uniqueCount`                  |
| min / max / avg / stdev           | `col.stats.{min,max,avg,stdev}`          |
| median / quartiles                | `col.stats.{med,q1,q2,q3}`               |
| sum                               | `col.stats.sum`                          |
| pairwise correlation              | `col.stats.corr(other)`                  |
| sorted distinct strings           | `col.categories`                         |
| typed-array raw view              | `col.getRawData()`                       |
| value iterator (generator)        | `col.values()`                           |

Do **not** compute these by hand:

- Total rows: `col.length`, not `col.toList().length` (allocates).
- Non-null count: `col.stats.valueCount`.
- Distinct count: `col.stats.uniqueCount`, not `new Set(col.toList()).size`.
- For string columns, `col.categories.length === col.stats.uniqueCount`.

```datagrok-exec
// Build a key-value summary of one numeric column.
const col = t.getCol('age');
const s = col.stats;
return ui.tableFromMap({
  name: col.name,
  type: col.type,
  semType: col.semType ?? '(none)',
  length: col.length,
  missing: s.missingValueCount,
  unique: s.uniqueCount,
  min: s.min, max: s.max, avg: s.avg, stdev: s.stdev,
});
```

```datagrok-exec
// Summarize every numeric column. t.columns.numerical is Iterable<Column>
// (not an Array) — convert via Array.from to get .map etc.
const rows = Array.from(t.columns.numerical, (c) => ({
  name: c.name, min: c.stats.min, max: c.stats.max, avg: c.stats.avg,
}));
return DG.Viewer.grid(DG.DataFrame.fromObjects(rows)).root;
```

### Top-N categories

No built-in top-N for a string column. Tally into a `Map<string, number>`:

```datagrok-exec
const col = t.getCol('class');
const counts = new Map();
for (let i = 0; i < col.length; i++) {
  const v = col.get(i);
  if (v === null || v === undefined) continue;
  counts.set(v, (counts.get(v) ?? 0) + 1);
}
const top = Array.from(counts.entries())
  .sort((a, b) => b[1] - a[1])
  .slice(0, 5);
return ui.tableFromMap(Object.fromEntries(top));
```

If you only need distinct values (no counts), `col.categories` is already
sorted alphabetically.

## Cloning DataFrames

`t.copy()` does **not** exist. Use
`t.clone(rowMask?, columnIds?, saveSelection?, saveTags?)`:

| Need                                | Call                                          |
|-------------------------------------|-----------------------------------------------|
| Full copy                           | `t.clone()`                                   |
| Filtered rows only                  | `t.clone(t.filter)`                           |
| Selected rows only                  | `t.clone(t.selection)`                        |
| Column subset (all rows)            | `t.clone(null, ['a', 'b'])`                   |
| Filtered rows + column subset       | `t.clone(t.filter, ['a', 'b'])`               |
| Preserve selection bits             | `t.clone(t.filter, null, true)`               |

```datagrok-exec
// Just SMILES and activity, filtered rows only.
const small = t.clone(t.filter, ['smiles', 'activity']);
return DG.Viewer.grid(small).root;
```

`col.clone(mask?)` clones a single column, but the cloned column is **not
attached to any DataFrame**. Call `df.columns.add(clonedCol)` to attach it.

## Adding columns

| You have...                              | Right call                                                      |
|------------------------------------------|-----------------------------------------------------------------|
| Nothing — want empty typed column        | `t.columns.addNewFloat(name)` (and `addNewInt/String/Bool/...`) |
| Generic typed empty column               | `t.columns.addNew(name, DG.COLUMN_TYPE.FLOAT)`                  |
| An array of values                       | `t.columns.add(DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, values))` |
| A function of row index                  | `t.columns.addNewFloat(name).init(i => ...)`                    |
| A formula using other columns            | See `datagrok-calc-column` skill                                |
| Compute-on-demand (no storage)           | `t.columns.addNewVirtual(name, i => ..., DG.TYPE.STRING)`        |
| Want a collision-free name first         | `const n = t.columns.getUnusedName(name);` then `addNew*`       |
| Insert at a specific position            | `t.columns.insert(col, index)`                                  |

> **`DG.COLUMN_TYPE.FLOAT === 'double'`.** Never write the literal `'float'`
> (no such type) or `'double'` (works but brittle). Use the constant, or the
> typed shorthand `addNewFloat`.

```datagrok-exec
// Typed empty column with metadata applied in the same block.
const col = t.columns.addNewFloat('Ki');
col.semType = DG.SEMTYPE.Ki;
col.meta.units = 'nM';
col.meta.format = '0.00';
return ui.divText(`Added ${col.name} (${col.type})`);
```

```datagrok-exec
// Add a column from an array of values. fromList infers length from the array.
const values = [1.2, 3.4, null, 5.6];
const col = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'score', values);
t.columns.add(col);
return ui.divText(`Added ${col.name}`);
```

```datagrok-exec
// Virtual column — recomputed on every access, not stored. Type argument is
// DG.TYPE.* (not DG.COLUMN_TYPE.*) — addNewVirtual takes the broader enum.
const labels = t.columns.addNewVirtual('Label',
  (i) => `${t.get('compound', i)}@${t.get('target', i)}`,
  DG.TYPE.STRING);
return ui.divText(`Added virtual column ${labels.name}`);
```

```datagrok-exec
// Avoid name collisions explicitly.
const name = t.columns.getUnusedName('Ki');  // returns 'Ki' or 'Ki (1)' etc.
const col = t.columns.addNewFloat(name);
return ui.divText(`Final name: ${col.name}`);
```

## Removing & renaming

```datagrok-exec
// Remove one column. Accepts string, index, or Column instance.
t.columns.remove('temp');
```

```datagrok-exec
// Remove several. Loop — no bulk remove call exists.
for (const name of ['_tmp1', '_tmp2', 'temp3'])
  if (t.columns.contains(name))
    t.columns.remove(name);
```

To target **calculated columns**, check `col.meta.formula`:

```datagrok-exec
// Drop every calculated column.
const names = t.columns.toList().filter((c) => c.meta.formula != null).map((c) => c.name);
for (const n of names) t.columns.remove(n);
return ui.divText(`Removed: ${names.join(', ') || '(none)'}`);
```

Rename: assign `col.name` directly.

```datagrok-exec
t.getCol('pIC50').name = 'pIC50_old';
```

## Setting column metadata

Use the right accessor — never `setTag` for these.

| Where metadata lives | Read                          | Write                                                       |
|----------------------|-------------------------------|-------------------------------------------------------------|
| Semantic type        | `col.semType`                 | `col.semType = DG.SEMTYPE.MOLECULE`                         |
| Friendly name        | `col.meta.friendlyName`       | `col.meta.friendlyName = 'Potency'` (or `null` to clear)    |
| Description          | `col.meta.description`        | `col.meta.description = '...'`                              |
| Format               | `col.meta.format`             | `col.meta.format = '0.00'`                                  |
| Units                | `col.meta.units`              | `col.meta.units = 'nM'`                                     |
| Choices (combo)      | `col.meta.choices`            | `col.meta.choices = ['A', 'B']`                             |
| Color coding         | `col.meta.colors.getType()`   | `col.meta.colors.setLinear/setCategorical/setConditional/setDisabled` |
| Arbitrary tags       | `col.getTag(k)`               | `col.setTag(k, v)`, `col.tags.delete(k)`                    |
| Layout-stable id     | `col.layoutColumnId`          | `col.layoutColumnId = '...'`                                |

Assign `null` to a `col.meta.*` accessor to **clear** the underlying tag.

> **Trap.** The semType tag key is `'quality'` in storage. Don't use
> `col.setTag('semType', ...)` or `col.setTag('quality', ...)`. Use the
> `col.semType` accessor — it sets the right tag and fires the right event.
> Friendly name / format / units / description go through `col.meta.*`, not
> through `setTag`.

```datagrok-exec
// Set several metadata fields at once.
const col = t.getCol('Ki');
col.semType = DG.SEMTYPE.Ki;
col.meta.units = 'nM';
col.meta.format = '0.00';
col.meta.description = null;  // null clears
```

```datagrok-exec
// Clear all custom (user-facing) metadata. Do NOT call col.tags.clear() —
// that also nukes system tags like .color-coding-type.
const col = t.getCol('activity');
col.meta.friendlyName = null;
col.meta.units = null;
col.meta.format = null;
col.meta.description = null;
```

## Color coding

Lives under `col.meta.colors`. Four shapes:

### Linear (numeric only — throws on string)

```datagrok-exec
// Red → yellow → green across the value range, with min/max pinned.
t.getCol('age').meta.colors.setLinear(
  ['#ff0000', '#ffff00', '#00ff00'],
  {min: 19, max: 70},
);
```

Hex strings and ARGB integers both work in `range`. Optional `belowMinColor`
/ `aboveMaxColor` paint out-of-range values.

### Categorical

```datagrok-exec
t.getCol('race').meta.colors.setCategorical(
  {'Asian': '#0000FF', 'Black': '#FF0000'},
  {fallbackColor: '#CCCCCC'},
);
```

### Conditional (range-based rules)

```datagrok-exec
// '20-170' means "value in the range 20..170". '<100' / '>50' also work.
t.getCol('height').meta.colors.setConditional({
  '20-170': '#00FF00',
  '170-190': '#220505',
});
```

### Off

```datagrok-exec
t.getCol('activity').meta.colors.setDisabled();
```

> **Trap.** Don't use `col.setTag('.color-coding-type', ...)` to enable /
> disable coloring. It bypasses rule validation and the companion tag writes
> that `setLinear` / `setCategorical` / `setConditional` perform, producing
> broken color coding.

## Working with cell values

Bulk init: `col.init(scalar | (i) => value)` is the right shape for an
already-allocated column.

```datagrok-exec
const col = t.columns.addNewInt('rowSquared');
col.init((i) => i * i);
```

For hot paths, mutate values with `notify=false` and call
`col.fireValuesChanged()` once at the end:

```datagrok-exec
const src = t.getCol('input');
const dst = t.getCol('output');
for (let i = 0; i < t.rowCount; i++)
  dst.set(i, src.get(i) * 2, false);
dst.fireValuesChanged();
```

## Constants — cite, never invent

| Constant family             | Members                                                              |
|-----------------------------|----------------------------------------------------------------------|
| `DG.COLUMN_TYPE`            | `STRING='string'`, `INT='int'`, **`FLOAT='double'`** (!), `BOOL='bool'`, `BYTE_ARRAY='byte_array'`, `DATE_TIME='datetime'`, `BIG_INT='bigint'`, `QNUM='qnum'`, `DATA_FRAME='dataframe'`, `OBJECT='object'` |
| `DG.SEMTYPE`                | `MOLECULE`, `MACROMOLECULE`, `MOLECULE3D`, `IC50`, `EC50`, `Ki`, `CONCENTRATION`, `VOLUME`, `EMAIL`, `URL`, `LATITUDE`, `LONGITUDE`, `IMAGE`, `FILE`, `CHEMICAL_REACTION`, ... |
| `DG.UNITS.Molecule`         | `SMILES`, `MOLBLOCK`, `V3K_MOLBLOCK`, `INCHI`                        |
| `DG.COLOR_CODING_TYPE`      | `CATEGORICAL`, `CONDITIONAL`, `LINEAR`, `OFF`                        |

## Anti-patterns — and the fix

1. **Looping `t.rows` to read / write column data.** Slow and
   notification-heavy. Use `t.getCol(name)` and an index loop with
   `set(..., false)` plus a final `col.fireValuesChanged()`.
2. **`col.setTag('color-coding-type', ...)` / `col.setTag('.color-coding-type', ...)`.**
   Use `col.meta.colors.setLinear/setCategorical/setConditional/setDisabled`.
3. **`col.setTag('semType', 'Molecule')`.** Use
   `col.semType = DG.SEMTYPE.MOLECULE`.
4. **`col.setTag('friendlyName', '...')`, `col.setTag('format', '...')`.**
   Go through `col.meta.friendlyName` / `col.meta.format`.
5. **`'float'` or `'double'` string literals for column type.** Always
   `DG.COLUMN_TYPE.FLOAT` (or just `addNewFloat`).
6. **`t.copy()`.** Doesn't exist. Use `t.clone()`.
7. **`col.toList().length` for row count.** Use `col.length`, or
   `col.stats.valueCount` for non-null count.
8. **`new Set(col.toList()).size` for unique count.** Use
   `col.stats.uniqueCount`. For string cols, `col.categories.length` also works.
9. **Mutating `col.getRawData()` and forgetting to fire.** End with
   `col.fireValuesChanged()` — otherwise viewers / grid don't refresh.
10. **Forgetting `await` on `t.columns.addNewCalculated`.** It's the only
    `addNew*` that returns a `Promise<Column>`. See `datagrok-calc-column`.
11. **`t.columns.byName(name)` + `if (col == null)` when you want an error.**
    Use `t.getCol(name)` — it throws on missing.
12. **`col.tags.clear()` to "reset metadata".** Nukes system tags too (color
    coding, choices, etc). Clear specific `col.meta.*` fields with `null`.
