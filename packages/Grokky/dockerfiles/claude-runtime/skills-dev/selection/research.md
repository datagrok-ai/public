# Row selection API research — for `datagrok-selection` skill

Citations use `<absolute-path>:<line>`. Paths abbreviated where unambiguous:
- `js-api/` = `/Users/aleksashka/Desktop/datagrok/reddata/public/js-api/src/`
- `samples/` = `/Users/aleksashka/Desktop/datagrok/reddata/public/packages/ApiSamples/scripts/`
- `chem/` = `/Users/aleksashka/Desktop/datagrok/reddata/public/packages/Chem/src/`
- `grokky/` = `/Users/aleksashka/Desktop/datagrok/reddata/public/packages/Grokky/src/`

Selection and filtering share the **exact same `BitSet` API** (`js-api/dataframe/bit-set.ts`). The full BitSet surface — `set/get`, `init`, `setAll`, `invert`, `and/or/xor/andNot`, `copyFrom`, `clone`, `findNext/Prev`, `getSelectedIndexes`, `trueCount/falseCount/length`, `onChanged`, `fireChanged` — is documented in topic 2's research (`grokky/dockerfiles/claude-runtime/skills-dev/filtering/research.md:24-65`). This note **cross-references** that surface rather than duplicating it, and focuses on what's *different* about selection: the opposite polarity for clearing, the lack of a "selection group" UI analog, and the cluster of "current row / cell / column" pointers that travel with it.

## 1. Summary table

| # | Area | Key entrypoints | Canonical idiom | Wrapper-worthiness |
|---|------|-----------------|-----------------|--------------------|
| A | The selection `BitSet` | `df.selection` (`js-api/dataframe/data-frame.ts:161`); shares full BitSet API (`bit-set.ts:18-222`) | `true` = "selected", `false` = "not selected". `trueCount` = selected count. | **High** — teach the polarity flip vs filter |
| B | Programmatic selection | `df.selection.init(pred)` (`bit-set.ts:170`); `df.rows.select(rowPred)` (`js-api/dataframe/row.ts:289`); `df.selection.set(i, true)` (`bit-set.ts:161`) | `df.selection.init(i => col.get(i) > 10)` for speed; `df.rows.select(row => row.x > 10)` for row-shape | **High** — wrap to absorb index-list, predicate, and BitSet inputs in one helper |
| C | Clearing & inverting | `df.selection.setAll(false)` (clears!); `.invert()`; `.and/or/andNot()` (`bit-set.ts:97-134`) | clear = `setAll(false)`; invert = `selection.invert()` | **High** — `clearSelection` exists to kill the polarity footgun |
| D | Reading selection state | `df.selection.trueCount/falseCount/anyTrue/anyFalse` (`bit-set.ts:71-84`); `.getSelectedIndexes()` (`bit-set.ts:188`); `.get(i)`; `df.rows.indexes({onlySelected: true})` (`row.ts:218`) | `df.selection.getSelectedIndexes()` for an array; `df.selection.trueCount` for a count | **Low** — wrap only as part of `selectedDf` etc. |
| E | Current row / cell / column | `df.currentRowIdx` get/set (`data-frame.ts:300`); `df.currentRow` (Row) (`:296`); `df.currentCol`/`currentCell` (`:309`,`:321`); `df.mouseOverRowIdx` (`:304`) | `df.currentRowIdx = 5` to focus; selection and "current" are independent | **Medium** — distinct concept; a small `setCurrentRow` helper + docs |
| F | Selection event lifecycle | `df.onSelectionChanged = df.selection.onChanged` (`data-frame.ts:490`); `.fireChanged()` (`bit-set.ts:198`); per-bit `notify` flags | One coalesced event per logical update; pass `notify: false` then `fireChanged()` for batches | **Low** — teach, don't wrap |
| G | Common idioms (df.clone, copy across masks) | `df.clone(df.selection, cols?, saveSelection?, saveTags?)` (`data-frame.ts:290`); `df.selection.copyFrom(df.filter)`; `df.filter.copyFrom(df.selection)` | "selected as DF" → `df.clone(df.selection)`; "select what's filtered" → `df.selection.copyFrom(df.filter)` | **High** — `selectedDf` + cross-skill `selectionFromFilter` / `filterFromSelection` |
| H | Existing `grokky.selectRows` | `grokky/src/claude/grokky-api.ts:38-48` — accepts `BitSet \| (i)=>bool \| number[]` | Already absorbs three shapes; missing: clear/invert/extend modes, "select all", view sync | **High** — extend with `mode: 'replace' \| 'add' \| 'remove'` |
| I | Multi-DF selection sync | `grok.data.linkTables(t1, t2, k1, k2, [SYNC_TYPE.SELECTION_TO_SELECTION, …])` (`js-api/data.ts:181`); `SYNC_TYPE` enum (`const.ts:42-53`) | One call wires master/detail sync of current/selection/filter | **Low** — out of scope for `datagrok-exec` lifecycle; document |
| J | Selected rows as DataFrame | `df.clone(df.selection)` (canonical) (`data-frame.ts:290`) | One-liner; identical shape to `df.clone(df.filter)` from topic 2 | **High** — `selectedDf(df, opts?)` mirroring `filteredDf` |
| K | Mouse / click selection (out of scope) | `BitSet.handleClick(rowPredicate, mouseEvent, modifiedSelectOnly)` (`bit-set.ts:220`); `df.mouseOverRowIdx` (`data-frame.ts:304`); `onMouseOverRowChanged` (`:430`) | Driven by grid; not for `datagrok-exec` | **None** — document only |

## 2. The selection BitSet (Area A)

`df.selection` is a `BitSet` exposed via getter (`data-frame.ts:161`):

```ts
get selection(): BitSet {
  return new BitSet(api.grok_DataFrame_Get_Selection(this.dart));
}
```

Same underlying class as `df.filter` — every method enumerated in topic 2's research table (`filtering/research.md:33-55`) applies verbatim. The convention difference is purely semantic:

| Mask | `true` means | `false` means | Clear (no-op state) |
|---|---|---|---|
| `df.filter` | "row passes filter" / visible | "row hidden" | `setAll(true)` — every row visible |
| `df.selection` | "row is selected" | "row not selected" | `setAll(false)` — nothing selected |

**This is the #1 footgun.** Claude will reflexively reach for `setAll(false)` to clear *both* masks, and will be right exactly half the time. The skill must lead with a side-by-side table.

Confirmed by samples and source:

- `samples/data-frame/bitset/select-rows-fast.js:8` selects evens via `selection.init((i) => i % 2 === 0)`.
- `samples/data-frame/aggregation/aggregate.js:3` does the same.
- `samples/data-frame/modification/manipulate.js:21-23` does `demog.selection.invert(); demog.selection.set(5, false); demog.selection.findNext(0, false);` — direct evidence the same BitSet surface is in active use for selection.
- `samples/data-frame/bitset/and-or-xor.js:1-13` demonstrates `clone().and/or/xor/andNot` chaining on a standalone `BitSet`, with `toBinaryString()` for inspection — the exact ops used to combine selections.

Reading "is row N selected?" → `df.selection.get(i)`. Counts via `df.selection.trueCount`/`falseCount` (`bit-set.ts:71-77`, both cached). Quick booleans `anyTrue`/`anyFalse` (`bit-set.ts:81-84`). Selected indices as `Int32Array` → `df.selection.getSelectedIndexes()` (`bit-set.ts:188`; result is cached, but the cache invalidates only on the public mutators — direct `getBuffer()` writes go stale, same warning as topic 2).

## 3. Programmatic selection (Area B)

Three idioms — same shape menu as filtering, with one extra "by index list" pattern:

1. **`df.selection.init(i => predicate(i))`** (`bit-set.ts:170`) — buffer-direct, single notification, fastest. Zeros the buffer first then writes — **replaces** the selection.
2. **`df.rows.select(row => predicate(row))`** (`row.ts:289`) — row-shape predicate. Implementation (`row.ts:280-291`):
   ```ts
   _applyPredicate(bitset, rowPredicate) {
     for (let row of this) bitset.set(row.idx, rowPredicate(row), false);
   }
   select(rowPredicate) {
     this._applyPredicate(this.table.selection, rowPredicate);
     this.table.selection.fireChanged();
   }
   ```
   Like filter's `rows.filter`, it iterates the `RowList` (slow) and `set(i, val, false)` per bit, then a single `fireChanged`. Functionally a "replace", but unlike `init` it does **not** zero the buffer — it overwrites each bit by the predicate's return value. Same end state, slower path. Source flags `RowList` as not-for-perf (`row.ts:191-192`).
3. **Per-bit `set` with deferred notify** (`samples/data-frame/bitset/select-rows-fast.js:12-14`):
   ```js
   for (let i = 0; i < 10; i++) selection.set(i, true, false);
   selection.fireChanged();
   ```
   The canonical "select these specific indices" pattern. The third arg `notify=false` (`bit-set.ts:161`) suppresses the per-bit event; `fireChanged()` (`bit-set.ts:198`) emits the coalesced one at the end.

Sample (`samples/data-frame/row-matching/select-rows.js:4-5`):

```js
demog.rows.select((row) => row.sex === 'M');
demog.rows.filter((row) => row.age > 42);
```

`RowMatcher` (`row.ts:128`) exposes `.select()`/`.filter()`/`.highlight()` on a Dart-parsed query string from `df.rows.match('age > 50')` (`row.ts:274`). Useful when the user phrasing is already a textual predicate, but less general than a JS function.

## 4. Clearing, inverting, combining (Area C)

| Intent | Code |
|---|---|
| clear selection (deselect all) | `df.selection.setAll(false)` |
| select all rows | `df.selection.setAll(true)` |
| invert | `df.selection.invert()` |
| AND another mask (intersect) | `df.selection.and(otherBitSet)` |
| OR | `df.selection.or(otherBitSet)` |
| XOR (symmetric difference) | `df.selection.xor(otherBitSet)` |
| subtract | `df.selection.andNot(otherBitSet)` |
| copy from another mask | `df.selection.copyFrom(otherBitSet)` |
| fresh standalone mask | `DG.BitSet.create(df.rowCount, i => …)` |

All mutating ops return `this`, chain freely. Each takes an optional `notify` flag (`bit-set.ts` all signatures).

`samples/demo/demo.js:16` and `samples/data-frame/modification/manipulate.js:22` both show `df.selection.invert()` in action — invert is a real, used idiom.

## 5. Reading selection state (Area D)

Three discriminated by performance and shape:

| Need | Code |
|---|---|
| count selected rows | `df.selection.trueCount` (cached) |
| count unselected | `df.selection.falseCount` |
| boolean "anything selected" | `df.selection.anyTrue` |
| boolean "any unselected" | `df.selection.anyFalse` |
| sorted indices | `const idx = df.selection.getSelectedIndexes()` → `Int32Array` |
| iterate as wu-stream | `for (const i of df.rows.indexes({onlySelected: true})) …` (`row.ts:218`) |
| "is row N selected?" | `df.selection.get(N)` |
| find next / prev set bit | `df.selection.findNext(i, true)` / `findPrev(i, true)` (`bit-set.ts:139-149`) |

Real usage in `chem/src/analysis/molecular-matched-pairs/mmp-viewer/mmp-grids.ts:388,400`:

```ts
const selectedRows = this.fpGrid.table.selection.getSelectedIndexes();
…
const indexesPairs = pairsGrid ? this.mmpGridTrans.dataFrame.selection.getSelectedIndexes() : …
```

`chem/src/utils/reaction-enumeration/enumerator-app.ts:639-648` uses the count guard before cloning, which is the standard pattern:

```ts
const sel = df.selection;
if (sel.trueCount === 0) { grok.shell.info('Select rows…'); return; }
if (sel.trueCount === df.rowCount) { grok.shell.info('All rows are selected — nothing to subset.'); return; }
const subset = df.clone(sel);
```

UI binding from `samples/ui/viewers/types/markup-advanced.js:50` and `samples/functions/custom-viewers/viewers.js:36`:

```
Selected: <span>#{t.selection.trueCount}</span>
Selected: ${this.dataFrame.selection.trueCount}<br>
```

— both confirm `trueCount` is the canonical "how many are selected" read.

## 6. Current row / cell / column (Area E)

Distinct concept from selection — Datagrok carries pointers to the *one* "focused" row, column, and cell:

| Property | Type | Setter | Source |
|---|---|---|---|
| `df.currentRowIdx` | `number` | yes | `data-frame.ts:300-301` |
| `df.currentRow` | `Row` (proxy) | yes (writes the idx) | `:296-297` |
| `df.currentCol` | `Column` | yes | `:309-317` |
| `df.currentCell` | `Cell` | yes | `:321-329` |
| `df.mouseOverRowIdx` | `number` | yes | `:304-305` |

Sample (`samples/data-frame/events/current-elements.js:1-12`):

```js
t.currentRow = 4;                  // ← Row idx assignment shortcut
t.currentCol = t.col('age');
t.currentCell = t.cell(5, 'sex');

t.onCurrentCellChanged.subscribe(_ => …);
t.onCurrentColChanged.subscribe(_ => …);
t.onCurrentRowChanged.subscribe(_ => …);
```

Note: `t.currentRow = 4` works because the setter takes a `Row` but accepts a number too via `row.idx`. Source signature is `set currentRow(row: Row)` (`data-frame.ts:297`) — Dart-side normalizes the `idx`. Safer to write `t.currentRowIdx = 4` from TS.

Cell has `.rowIndex` and `.column.name` (`row.ts:358-365`), which are the two properties Claude will reach for when answering "what cell is selected?".

**Selection ≠ current row.** They're independent state:
- Setting `df.currentRowIdx = 5` does **not** select row 5.
- Selecting rows does **not** move the current row.
- `SYNC_TYPE.CURRENT_ROW_TO_SELECTION` exists precisely because the link is opt-in (`const.ts:45`).

Three events:

| Event | Fires when |
|---|---|
| `df.onCurrentRowChanged` (`data-frame.ts:427`) | `currentRowIdx` mutates |
| `df.onCurrentColChanged` (`:433`) | `currentCol` mutates |
| `df.onCurrentCellChanged` (`:439`) | either of the above changes the cell |
| `df.onMouseOverRowChanged` (`:430`) | hover |

Subscribe lifecycle is the same as filter: in `datagrok-exec` the closure dies with the block (no leak); in widget/viewer code push onto `viewer.subs[]` (`js-api/viewer.ts:398`).

## 7. Selection event lifecycle (Area F)

Only one event — `df.onSelectionChanged`, which is just `df.selection.onChanged` (`data-frame.ts:490`):

```ts
get onSelectionChanged(): Observable<any> { return this.selection.onChanged; }
```

No `onSelectionChanging` / no two-phase lifecycle, **unlike filter**. The filter system has `onRowsFiltering` (running) and `onRowsFiltered` (done) because filters can collaborate via `onRowsFiltering` (custom filters AND their contribution onto `df.filter`). Selection doesn't have collaborators — it's just a piece of state, so there's no need.

`samples/data-frame/events/events.js:30`:
```js
demog.onSelectionChanged.subscribe(_ => info('ddt-selection-changed'));
```

**Coalescing**: every `BitSet` mutator (`set`, `setAll`, `init`, `invert`, `and/or/xor/andNot`, `copyFrom`) takes `notify: boolean = true`. The default fires immediately on each call (`bit-set.ts:97-196`). The deferral pattern is from `samples/data-frame/bitset/select-rows-fast.js:11-14`:

```js
for (let i = 0; i < 10; i++) selection.set(i, true, false);
selection.fireChanged();
```

So:
- `df.selection.init(...)` → 1 event (default).
- A loop of `set(i, true)` without the `false` flag → **N events**, slow.
- A loop of `set(i, true, false)` then a single `fireChanged()` → 1 event, fast.
- `df.selection.setAll(false).and(mask)` → 2 events (one per mutator). To coalesce, pass `notify: false` to the first.

The `BitSet` exposes `version` (`bit-set.ts:87`) that bumps per mutation — useful in widgets for "did anything change since I last looked?" without rebinding.

## 8. Common idioms from samples & packages (Area G)

| Intent | Code | Citation |
|---|---|---|
| Selected rows as a DF | `const sub = df.clone(df.selection);` | `chem/src/utils/reaction-enumeration/enumerator-app.ts:648` |
| Selected rows, subset columns | `df.clone(df.selection, ['smiles', 'activity'])` | `data-frame.ts:284-292` |
| Selected rows, keep selection on clone | `df.clone(df.selection, null, true)` | `data-frame.ts:290` (third arg `saveSelection`) |
| Select all currently visible | `df.selection.copyFrom(df.filter)` | derives from polarity match: filter `true` = visible = should be selected |
| Push selection into filter | `df.filter.copyFrom(df.selection)` | cross-skill (topic 2) |
| Combine two selections | `dfA.selection.or(dfB.selection)` | only safe if `rowCount` matches |
| "Select where" predicate | `df.selection.init(i => df.getCol('x').get(i) > 10)` | `bit-set.ts:170` |
| Toggle row | `df.selection.set(i, !df.selection.get(i))` | direct |
| Select indices | `for (i of arr) df.selection.set(i, true, false); df.selection.fireChanged();` | `samples/data-frame/bitset/select-rows-fast.js:11-14` |

The "selected rows as DF" pattern is overwhelmingly the common one in production code (`chem/utils/reaction-enumeration/enumerator-app.ts:639-704` repeats it three times with identical structure). Worth a one-liner helper.

## 9. Existing `grokky.selectRows` — current state and gaps (Area H)

Source (`grokky/src/claude/grokky-api.ts:38-48`):

```ts
function selectRows(df: DG.DataFrame, predicate: DG.BitSet | ((i: number) => boolean) | number[]): void {
  if (predicate instanceof DG.BitSet)
    df.selection.copyFrom(predicate);
  else if (typeof predicate === 'function')
    df.selection.init(predicate);
  else {
    df.selection.setAll(false, false);
    for (const i of predicate) df.selection.set(i, true, false);
    df.selection.fireChanged();
  }
}
```

What it does well:
- Three input shapes in one call (BitSet, predicate, index array).
- For the index-array path: zeros first, batches writes with `notify=false`, single `fireChanged` — the textbook fast path.
- Uses `copyFrom` for BitSet — preserves Dart-side identity.

What it lacks:
- **Replace-only semantics.** No `'add'` (OR), `'remove'` (AND NOT), or `'intersect'` (AND) modes. Common requests: "also select where X", "deselect rows where Y", "select the intersection of A and B".
- **No clear / invert / select-all helpers** — those have to be reached for through `df.selection.setAll(false)` (the polarity trap) or `df.selection.invert()`.
- **No "selected as DataFrame" wrapper.** Users will write `df.clone(df.selection)` ad hoc; standardizing kills two bugs: `saveSelection` flag forgotten, and column subsetting confusion.
- **No selection ↔ filter bridges.** Topic 2's `filterFromSelection` is unwritten; the inverse `selectionFromFilter` is equally common (chem flows, reaction enumeration).
- **Doesn't handle `Int32Array` from `getSelectedIndexes`.** The `number[]` branch loops via `for..of`, which works on `Int32Array` (it's iterable), but the type annotation excludes it — TS users get yellow squiggles when they pipe `getSelectedIndexes()` straight in.
- **Doesn't validate index bounds.** Out-of-range writes silently no-op at the BitSet level (no error), but a debug warning would help.

## 10. Multi-DF selection sync (Area I)

`grok.data.linkTables(t1, t2, keyCols1, keyCols2, linkTypes, initialSync=false, filterAllOnNoRowsSelected=false)` (`js-api/data.ts:181`). `SYNC_TYPE` enum (`const.ts:42-53`):

```
CURRENT_ROW_TO_ROW           = 'row to row'
CURRENT_ROW_TO_SELECTION     = 'row to selection'
CURRENT_ROW_TO_FILTER        = 'row to filter'
MOUSE_OVER_ROW_TO_SELECTION  = 'mouse-over to selection'
MOUSE_OVER_ROW_TO_FILTER     = 'mouse-over to filter'
FILTER_TO_FILTER             = 'filter to filter'
FILTER_TO_SELECTION          = 'filter to selection'
SELECTION_TO_FILTER          = 'selection to filter'
SELECTION_TO_SELECTION       = 'selection to selection'
```

Sample (`samples/data-frame/join-link/link-tables.js:20-26`) shows master/detail wiring via `[SYNC_TYPE.CURRENT_ROW_TO_FILTER, SYNC_TYPE.MOUSE_OVER_ROW_TO_SELECTION]`.

**Out of scope** for `datagrok-exec`-style one-shot ops (linking establishes long-lived behavior that persists past the block), but worth a footnote in the skill: "if the user says 'when I click in table A, select in table B', that's `linkTables(...)`, not a `selectRows` call". The SKILL.md can include a 4-line example for discoverability.

## 11. Selected rows as DataFrame (Area J)

The canonical pattern is `df.clone(df.selection)` (`data-frame.ts:290`):

```ts
clone(rowMask: BitSet | null = null, columnIds: string[] | null = null, saveSelection: boolean = false, saveTags: boolean = true): DataFrame
```

Identical shape to `df.clone(df.filter)` from topic 2 — same args, just a different mask. Topic 2 wraps this as `filteredDf(df, opts)`; for symmetry, topic 3 should wrap as `selectedDf(df, opts)`.

`saveSelection` (third arg): if `true`, the new DF inherits the selection mask scoped to its rows. For "give me the selected rows" the typical answer is `false` (the new DF starts clean), but expose the knob.

`saveTags` (fourth arg, default `true`): tag forwarding. The wrapper should pass through transparently.

Empty-selection guard: returns an empty DF (zero rows, columns intact). The chem code (`chem/src/utils/reaction-enumeration/enumerator-app.ts:640-647`) emits a user-facing toast before calling clone — the wrapper should match that pattern when called from a `view` context, or at least return early without allocating.

## 12. Mouse / click / interactive selection — out of scope (Area K)

`BitSet.handleClick(rowPredicate, mouseEvent, modifiedSelectOnly)` (`bit-set.ts:220-222`):

```ts
handleClick(rowIndexPredicate, mouseEvent, modifiedSelectOnly = false) {
  api.grok_Utils_SelectRowsWhere(this.dart, rowIndexPredicate, mouseEvent.ctrlKey, mouseEvent.shiftKey, mouseEvent.metaKey, modifiedSelectOnly);
}
```

This is the implementation hook for custom viewers that want grid-style ctrl/shift/cmd-click semantics. Out of scope for `datagrok-exec` blocks (no mouse event in hand), but worth one bullet in the skill so Claude knows it exists if the user describes "shift-click to extend".

`df.mouseOverRowIdx` (`data-frame.ts:304-305`) is also runtime, not script-time.

## 13. Proposed grokky helpers

Signatures only — rationale below. The shape mirrors topic 2's `helper.ts` where possible.

```ts
// 13.1 — extended; absorbs replace / add / remove / intersect modes and Int32Array input
export type SelectMode = 'replace' | 'add' | 'remove' | 'intersect';

export type SelectInput =
  | DG.BitSet
  | ((i: number) => boolean)
  | ArrayLike<number>;   // number[] | Int32Array | Uint32Array

export function selectRows(
  df: DG.DataFrame,
  input: SelectInput,
  opts?: {mode?: SelectMode},
): void;

// 13.2 — deselect everything (polarity trap killer)
export function clearSelection(df: DG.DataFrame): void;

// 13.3 — select all rows
export function selectAll(df: DG.DataFrame): void;

// 13.4 — flip mask in place
export function invertSelection(df: DG.DataFrame): void;

// 13.5 — explicit "select these indices" with batched writes
export function selectByIndices(
  df: DG.DataFrame,
  indices: ArrayLike<number>,
  opts?: {mode?: SelectMode},   // default 'replace'
): void;

// 13.6 — predicate variant; thin sugar that telegraphs polarity in the name
export function selectByPredicate(
  df: DG.DataFrame,
  pred: (i: number) => boolean,
  opts?: {mode?: SelectMode},   // default 'replace'
): void;

// 13.7 — "give me a real DF of just the selected rows"
export type SelectedDfOpts = {
  cols?: string[];
  saveSelection?: boolean;
};
export function selectedDf(df: DG.DataFrame, opts?: SelectedDfOpts): DG.DataFrame;

// 13.8 — cross-skill bridges
export function selectionFromFilter(df: DG.DataFrame): void;
export function filterFromSelection(df: DG.DataFrame): void;

// 13.9 — current row / cell pointer (distinct from selection!)
export function setCurrentRow(df: DG.DataFrame, rowIdx: number): void;

// 13.10 — snapshot for debugging / reporting
export function describeSelection(df: DG.DataFrame): {
  rowCount: number;
  selectedCount: number;
  unselectedCount: number;
  selectedIndexes: number[];      // truncated to first ~100 with elided marker
  currentRowIdx: number;
};
```

### Rationale

- **13.1 `selectRows` (extended).** Keep the existing three-shape input; add `opts.mode`:
  - `'replace'` (default) — current behavior. Predicate → `init`; BitSet → `copyFrom`; index list → `setAll(false)` + batched sets.
  - `'add'` — OR onto current selection. Predicate → build a fresh `DG.BitSet.create(df.rowCount, pred)` then `df.selection.or(bs)`. BitSet → `df.selection.or(bs)`. Indices → loop `set(i, true, false)` then `fireChanged()`.
  - `'remove'` — AND NOT. Build fresh BitSet, `df.selection.andNot(bs)`. For indices, loop `set(i, false, false)`.
  - `'intersect'` — AND. Build fresh BitSet, `df.selection.and(bs)`.
  Accepts `Int32Array` (typed annotation `ArrayLike<number>`) so `getSelectedIndexes()` results pipe back in cleanly. **No `view` overload** — selection is a DF concept; if the caller has a view, `view.dataFrame` is one property away.
- **13.2 `clearSelection`.** One-liner (`df.selection.setAll(false)`). Justified entirely by the polarity trap: every Claude turn that wants to "clear" something will reach the same generic shape and pick the wrong sign for half of them. Named helper kills that.
- **13.3 `selectAll`.** Symmetric to clear. `df.selection.setAll(true)`. Useful for "select all then deselect X" patterns.
- **13.4 `invertSelection`.** `df.selection.invert()`. Tiny but discoverable.
- **13.5 `selectByIndices`.** Explicit shape for the most common "I have these row indices, select them" call. Internally same as `selectRows(df, idx, opts)`, but documenting separately catches users who don't realize `selectRows` takes arrays.
- **13.6 `selectByPredicate`.** Same logic, predicate-shaped. The polarity note here is **the same direction as `filterByPredicate`** — predicate returns `true` for "include in selection". That symmetry between the two skills is helpful; the only difference is the BitSet polarity for clearing, not for predicates.
- **13.7 `selectedDf`.** Mirrors `filteredDf` from topic 2. Wraps `df.clone(df.selection, opts?.cols ?? null, opts?.saveSelection ?? false)`. If `df.selection.trueCount === 0`, return an empty clone without warning (caller may legitimately want that); the SKILL.md should advise an explicit `if (df.selection.anyTrue)` guard for UX.
- **13.8 `selectionFromFilter` / `filterFromSelection`.** Cross-skill bridges. Both are one-line `copyFrom` calls but answer a frequent verbal request: "select what's currently visible" / "filter down to selection". Worth naming.
- **13.9 `setCurrentRow`.** `df.currentRowIdx = rowIdx`. Justified because the user request "highlight row 5" is ambiguous between "select row 5" and "navigate to row 5"; named helpers force the disambiguation. Out-of-bounds index is a no-op in the Dart layer; the wrapper can `console.warn` if `rowIdx < 0 || rowIdx >= df.rowCount`.
- **13.10 `describeSelection`.** Read-only snapshot. Truncated index array (first 100, plus `…N more`) to keep output sane in agent context. Includes `currentRowIdx` so the agent doesn't conflate the two.

**Helpers I considered and rejected:**

- `selectIntersect(df, other)` / `selectUnion(df, other)` — folded into `selectRows(..., {mode})`.
- `getSelectedIndexes(df)` — `df.selection.getSelectedIndexes()` is already terse; `describeSelection` covers the discovery angle.
- `onSelectionChanged(df, cb)` subscription wrapper — same reasoning as topic 2: sub leaks are a viewer-scope concern, not a `datagrok-exec` concern.
- `selectMatching(df, query)` — wraps `df.rows.match(query).select()`. Tempting, but the DSL syntax is Dart-parsed and undocumented from TS; better to teach predicate form.
- `linkTablesSelection(t1, t2, …)` — out of scope per Area I; SKILL.md gets a code example, not a helper.

## 14. Anti-patterns Claude must avoid

1. **`df.selection.setAll(true)` to "clear" the selection.** That selects every row — the *opposite* of what filter's clear does. Use `df.selection.setAll(false)` or `grokky.clearSelection(df)`. **(The headline footgun.)**
2. **Confusing selection polarity with filter polarity in the same script.** Side-by-side:
   ```ts
   df.filter.setAll(true);     // clear filter (show all rows)
   df.selection.setAll(false); // clear selection (none selected)
   ```
   The skill must lead with this exact pair.
3. **Looping `df.rows` to select.** Use `df.selection.init(i => …)` for predicates or `selection.set(i, true, false)` + `fireChanged()` for index lists. RowList iteration is slow and explicitly warned against (`row.ts:191-192`).
4. **Setting `df.currentRowIdx` when you meant to select.** Current row is a single-row pointer (focused row). Selection is a set. They're independent — `SYNC_TYPE.CURRENT_ROW_TO_SELECTION` is opt-in (`const.ts:45`).
5. **`df.clone(df.selection)` with empty selection.** Returns a zero-row DF, no error, no toast. If the operation is user-facing, guard with `df.selection.anyTrue` first (matches chem package practice, `chem/src/utils/reaction-enumeration/enumerator-app.ts:640-643`).
6. **Per-bit `set` in a loop without `notify=false`.** N writes → N events → laggy UI. Always `set(i, x, false)` in the loop and `fireChanged()` at the end (`samples/data-frame/bitset/select-rows-fast.js:11-14`).
7. **`df.selection.init(...)` to "add to" the selection.** `init` zeros the buffer first (`bit-set.ts:174-181`). To extend, OR a freshly-built `DG.BitSet.create(df.rowCount, pred)` onto `df.selection` — or use `selectRows(df, pred, {mode: 'add'})`.
8. **Subscribing to `df.onSelectionChanged` inside a `datagrok-exec` block.** The closure dies with the block — sub goes with it (no leak, but also no effect after the block ends). For long-lived selection-driven reactions, the caller needs widget code with `viewer.subs[]`. Document, don't try to fix in-block.
9. **Assuming `getSelectedIndexes()` is fresh after a direct `getBuffer()` write.** It's cached and invalidates only on the public mutators (`bit-set.ts:188`, same warning as topic 2). Use the public API.
10. **Combining selections across DataFrames of different row counts.** `df.selection.and(otherDf.selection)` only makes sense if rows are aligned 1:1. The Dart side may not validate — silent garbage. If the user wants cross-table sync by key columns, that's `grok.data.linkTables(...)` (Area I).
11. **`df.rows.select(rowPred)` for perf-critical code.** Iterates `RowList` (slow). Same warning as topic 2's `df.rows.filter` (`row.ts:191-192`). Prefer `df.selection.init(i => …)` reading columns directly.

## 15. Open questions

1. **`selectRows` mode default.** Keep `'replace'` (matches current behavior, principle of least surprise). Worth a side-comment in the JSDoc that this is also the safest default for naive callers.
2. **`selectByIndices` vs `selectRows` with array input.** Two helpers for the same thing? `selectByIndices` improves discoverability ("how do I select these specific rows?") but is redundant with `selectRows(df, idx)`. Lean: keep both, with `selectByIndices` documented as the "explicit shape" and `selectRows` as the "polymorphic". Each links to the other in JSDoc.
3. **`Int32Array` typing.** `ArrayLike<number>` accepts `number[]` and `Int32Array` and `Uint32Array`. Verified by spec — both typed arrays satisfy `ArrayLike<number>`. Acceptable lint surface? Alternative: `readonly number[] | ArrayLike<number>` for clarity, but `ArrayLike<number>` is sufficient and shorter.
4. **`setCurrentRow` overload taking a `Row`?** `df.currentRow = row` works (`data-frame.ts:297`); but in `datagrok-exec` callers will almost always have an index. Number-only signature is cleaner. If a Row pops up, `df.currentRowIdx = row.idx` is one keystroke.
5. **`selectionFromFilter` polarity asymmetry.** `df.selection.copyFrom(df.filter)` literally copies bit-for-bit, which works because both masks use `true`-as-positive in their respective domains ("visible" maps to "selected"). The names line up. Worth a one-line comment in JSDoc to reassure callers who think about the polarity flip from clearing semantics.
6. **`describeSelection.selectedIndexes` truncation threshold.** 100 indices is arbitrary. The agent context is large (Opus 4.7 1M), so truncation is for stdout readability rather than token economics. Could be configurable via opts, but probably not worth the surface area for a debug helper.
7. **`selectAll` and visible rows.** Should `selectAll` select all rows in the DF, or all *currently filtered* rows? Ambiguous. The literal `setAll(true)` selects every row regardless of filter — which is what `selection.setAll(true)` does at the BitSet level. The "select all visible" intent maps to `selectionFromFilter`. Worth a JSDoc disambiguation: "selects every row in the DF, including filtered-out ones. For 'select all visible', use `selectionFromFilter`."
8. **`linkTables` exposure.** Topic 3 doesn't include a wrapper for `linkTables`, but should the SKILL.md (not helper.ts) include a runnable example for "make table B show only what's selected in table A"? Probably yes — discoverability — but the wrapper would just be `grok.data.linkTables(...)`, no value added.
9. **Selection BitSet identity stability.** `df.selection` getter constructs a fresh `new BitSet(...)` each call (`data-frame.ts:161-163`). Subscribing on `df.selection.onChanged` should still work because the underlying Dart object is the same — but if the agent caches the BitSet wrapper and the user later does `df.selection.setAll(false)`, both wrappers point at the same Dart bitset. Confirm with a quick experiment before relying on it in eval prompts.
10. **`mode: 'replace'` with index list — empty array clears.** Current `selectRows` with `[]` does `setAll(false)` + zero set-bit calls + `fireChanged()` → empty selection. Intuitive, but is that documented? Probably yes — the SKILL.md examples should include the empty-array case explicitly.
