# Row filtering API research — for `datagrok-filtering` skill

Citations use `<absolute-path>:<line>`. Paths abbreviated where unambiguous:
- `js-api/` = `/Users/aleksashka/Desktop/datagrok/reddata/public/js-api/src/`
- `samples/` = `/Users/aleksashka/Desktop/datagrok/reddata/public/packages/ApiSamples/scripts/`
- `chem/` = `/Users/aleksashka/Desktop/datagrok/reddata/public/packages/Chem/src/`
- `grokky/` = `/Users/aleksashka/Desktop/datagrok/reddata/public/packages/Grokky/src/`

## 1. Summary table

| # | Area | Key entrypoints | Canonical idiom | Wrapper-worthiness |
|---|------|-----------------|-----------------|--------------------|
| A | `BitSet` API | `setAll`, `set/get`, `init`, `invert`, `and/or/xor/andNot`, `copyFrom`, `clone`, `findNext/Prev`, `getSelectedIndexes`, `trueCount/falseCount/length`, `onChanged`, `fireChanged` (`js-api/dataframe/bit-set.ts:18-222`) | `df.filter.init(i => col.get(i) > 10)`; clear with `df.filter.setAll(true)`; invert with `df.filter.invert()` | **High** — teach the surface, not wrap; semantics + clear/invert are footguns |
| B | Clearing / inverting / combining | `df.filter.setAll(true)`, `df.filter.invert()`, `df.filter.and/or/andNot(other)` (all in-place, return `this`) (`bit-set.ts:97-134`) | clear = `setAll(true)`; AND another mask = `df.filter.and(bs)` | **Medium** — one-liner helpers `clearFilter(df)`, `invertFilter(df)` |
| C | Predicate-based filtering | `df.filter.init(i => …)` (`bit-set.ts:170`); `df.rows.filter(row => …)` (`js-api/dataframe/row.ts:297`) | `df.filter.init` is the canonical, fastest form | **High** — wrap as `filterByPredicate(df, i => …)` to encode the gotcha that "true = passes" |
| D | Event lifecycle | `onRowsFiltering` (run filters), `onRowsFiltered` (filters done), `onFilterChanged` = `filter.onChanged` (`js-api/dataframe/data-frame.ts:472-493`); `rows.requestFilter()` (`row.ts:309`) | viewer: subscribe `onRowsFiltering`, edit `df.filter` in-place; outside: subscribe `onRowsFiltered` to read result; trigger by calling `df.rows.requestFilter()` | **Low** — teach the difference; one-liner `onFilterDone(df, fn)` if requested |
| E | Existing `grokky.filterRows` | `view.getFiltersGroup().updateOrAdd({type, column, …})` (`grokky-api.ts:20-31`) | passes substructure string straight as `molBlock` → broken for SMILES/SMARTS | **High** — must fix substructure path |
| F | Substructure filter — robust | `DG.chem.isMolBlock` (`js-api/chem.ts:68`), `DG.chem.isSmarts` (`chem.ts:904`), `DG.chem.convert` (`chem.ts:888`), `grok.chem.searchSubstructure(col, pattern)` (`chem.ts:766`) | Either (a) convert input → molblock then push through `updateOrAdd`; or (b) call `searchSubstructure` and write BitSet onto `df.filter` directly | **High** — UI path vs BitSet path are both valid; default to UI |
| G | Filter UI (FilterGroup) | `view.getFiltersGroup({createDefaultFilters?})` returns `FilterGroup` (`js-api/views/view.ts:404`); `fg.add/updateOrAdd/setEnabled/setActive/remove/getStates/filters` (`js-api/viewer.ts:510-553`); `view.filters()` (deprecated) (`view.ts:491`) | `view.getFiltersGroup({createDefaultFilters:false}).updateOrAdd({type, column, …})` | **Medium** — `filterRows` already wraps the common case |
| H | Multi-condition / combining | `df.filter.init(i => predA(i) && predB(i))` or build separate BitSets and `df.filter.and(bs).and(bs2)` | predicate AND is one-liner; programmatic AND across pre-built BitSets uses `BitSet.and` | **Medium** — fold into `filterByPredicate` |
| I | Show-only-filtered vs drop-rows | "show only" = mutate `df.filter`; "drop" = `df.clone(df.filter)` (`js-api/dataframe/data-frame.ts:290`) or `df.rows.removeWhereIdx(i => !df.filter.get(i))` (`row.ts:238`) | very different intent; the skill must distinguish | **High** — explicit helpers `keepRowsView(view, pred)` vs `dropRows(df, pred)` |
| J | Subscribe to filter events | `df.onFilterChanged`, `df.onRowsFiltered`, `df.onRowsFiltering` are `Observable<any>`; subs in datagrok-exec leak unless caller pushes onto `view.subs` | inside `datagrok-exec` block subs die with the block; meaningful only in widget/viewer code | **Low** — note in skill, no helper |

## 2. The BitSet API (Area A)

`df.filter` is a `BitSet` constructed once in the DataFrame ctor (`js-api/dataframe/data-frame.ts:51`, `:62`). One bit per row, 0-indexed.

**Crucial convention: `true` = "passes filter", `false` = "filtered out".**
- Clearing the filter (show all rows) → `df.filter.setAll(true)` (`bit-set.ts:132`).
- Filtering everything out → `df.filter.setAll(false)`.

Full surface (`bit-set.ts`, all methods):

| Method | Line | Notes |
|---|---|---|
| `length` (getter) | `:66` | row count of the parent DF |
| `trueCount`, `falseCount` | `:71-77` | cached |
| `anyTrue`, `anyFalse` | `:81-84` | trivial booleans |
| `version` | `:87` | bumps on each in-place mutation |
| `clone()` | `:91` | returns a new BitSet |
| `invert(notify=true)` | `:97` | in-place; returns `this` |
| `and(other, notify=true)` | `:104` | in-place; returns `this` |
| `or(other, notify=true)` | `:111` | in-place; returns `this` |
| `xor(other, notify=true)` | `:118` | in-place; returns `this` |
| `andNot(other, notify=true)` | `:125` | in-place; returns `this` |
| `setAll(x, notify=true)` | `:132` | in-place; returns `this` |
| `findNext(i, x)` / `findPrev(i, x)` | `:139-149` | -1 if not found |
| `get(i)` / `set(i, x, notify=true)` | `:152-163` | direct bit access |
| `init(predicate, notify=true)` | `:170` | bulk write via fn(i) — buffer-direct; fastest path |
| `getSelectedIndexes()` | `:188` | `Int32Array` of set bits; result is cached |
| `copyFrom(other, notify=true)` | `:193` | in-place; returns `this` |
| `fireChanged()` | `:198` | manual notify |
| `onChanged` (getter) | `:203` | `Observable<any>` |
| `similarityTo(other, metric)` | `:211` | Tanimoto etc. |
| `setAll`, `init`, `and`, `or`, `xor`, `andNot`, `invert`, `copyFrom` all take an optional `notify` flag | — | pass `false` for batch then `fireChanged()` |

`BitSet.create(length, f?)` (`bit-set.ts:47`) returns a fresh standalone BitSet — useful for building masks not yet attached to a DF.

**Selection uses the same BitSet API.** `df.selection` is also a `BitSet`, conventions identical, so any reflex Claude learns here transfers to topic 3. Per the existing `grokky.selectRows` (`grokky-api.ts:38-48`) the same `BitSet | (i)=>boolean | number[]` shape works for both.

**Gotchas**:
- `setAll(true)` clears the filter (counter-intuitive — Claude will reach for `setAll(false)`).
- Bitwise ops mutate `this` and return `this` for chaining; they do **not** return a new BitSet.
- `init` *replaces* — it zeros the buffer first (`bit-set.ts:175`). To merge with existing filter, use `and` with a freshly-built BitSet.
- `getSelectedIndexes()` is cached; if you mutate via direct buffer write (`getBuffer()`), the cache is stale.

## 3. Clearing / inverting / combining (Area B)

Three idioms Claude will need constantly:

| Intent | Code |
|---|---|
| clear all filters (show every row) | `df.filter.setAll(true)` |
| invert current filter | `df.filter.invert()` |
| AND another mask onto current filter | `df.filter.and(otherBitSet)` |
| OR | `df.filter.or(otherBitSet)` |
| subtract (filter out these rows from current) | `df.filter.andNot(otherBitSet)` |
| make a fresh mask | `DG.BitSet.create(df.rowCount, i => …)` |

All `BitSet` mutating ops return `this`, so chaining works: `df.filter.setAll(true).and(catA).and(catB)`.

**The "clear" idiom is the audit's #1 gap.** The skill must lead with it.

In samples: `samples/data-frame/filtering/collaborative-filtering.js:5-9` uses `bs.init(i => …)` inside `onRowsFiltering`. Substructure samples (`samples/domains/chem/substructure-search.js:6`, `substructure-search-smarts.js:5`, `substructure-search-smarts-molfile-als.js:28`) use `t.filter.copyFrom(bs)` to overwrite the whole filter from a freshly-computed BitSet. `samples/performance/dataframe-access.js:16` shows `table.filter.getSelectedIndexes()` is the fast iteration path.

## 4. Predicate-based filtering (Area C)

Three ways to filter rows by a predicate; only two are recommended:

1. **`df.filter.init(i => …)`** (`bit-set.ts:170`) — buffer-direct, single notification, ~10× faster than `df.rows.filter`. **This is the canonical form.**
2. **`df.rows.filter(row => …)`** (`row.ts:297`) — uses `_applyPredicate` (`row.ts:280`), which calls `bitset.set(idx, …, false)` per row then `fireChanged`. Cleaner row-shape but slower; explicitly flagged as not for performance-critical paths in the source (`row.ts:191-192`).
3. **`for (let i = 0; i < n; i++) df.filter.set(i, p(i), false); df.filter.fireChanged();`** — same as `init` semantically; just write `init`.

Sample (`samples/data-frame/row-matching/select-rows.js:4-5`):

```js
demog.rows.select((row) => row.sex === 'M');
demog.rows.filter((row) => row.age > 42);
```

Pattern matcher: `df.rows.match('age > 50').select()` / `.filter()` / `.highlight()` (`row.ts:128-130`, `samples/data-frame/row-matching/patterns.js:6`). This is a `RowMatcher` — string DSL parsed by Dart side. Useful but less general than a JS predicate.

**Anti-pattern: looping `df.rows`** for read access — same warning as topic 1 (`row.ts:191-192`). For predicate filtering, *only* `df.filter.init` reads columns directly.

## 5. Filter event lifecycle (Area D)

Three events on DataFrame:

| Event | Fires when | Use |
|---|---|---|
| `onRowsFiltering` (`data-frame.ts:475`) | filter system is rebuilding the mask; a viewer's chance to AND its contribution onto `df.filter` | implementing a custom filter that participates in collaborative filtering |
| `onRowsFiltered` (`data-frame.ts:472`) | filter pass complete; mask is now final | observers that act on the filtered result (count, recompute viewer) |
| `onFilterChanged` (`data-frame.ts:493`) | `df.filter.onChanged` — fires on any BitSet mutation incl. direct writes | most general; equivalent to listening on the BitSet itself |

`df.rows.requestFilter()` (`row.ts:309`): asks the filter system to re-run. Use after **changing the inputs** to a custom filter (state outside `df.filter` itself). If you wrote directly to `df.filter` you don't need it — the `onChanged` fires automatically.

Collaborative filtering pattern from `samples/data-frame/filtering/collaborative-filtering.js:5-9`:

```js
view.dataFrame.onRowsFiltering.subscribe((_) => {
  let bs = view.dataFrame.filter;        // already AND-folded by previous filters
  let sex = view.dataFrame.col('sex');
  bs.init((i) => sex.get(i) === 'M');    // contribute our mask
});
```

The docstring on `requestFilter` (`row.ts:307-308`) makes the contract explicit: "Viewers that filter rows should subscribe to DataFrame.onRowsFiltering event. When filtering conditions are changed, viewers should call requestFilter()."

Custom-filter widgets (the `Filter` class, `js-api/widgets/filter.ts:15`) override `applyFilter()` (`:75`) and AND their contribution onto `dataFrame.filter`. Then call `dataFrame.rows.addFilterState(this.saveState())` (`row.ts:315`) so the state is round-trippable via layout save/load. Out of scope for `datagrok-exec` but useful context for "why does my filter get clobbered" questions.

**Subscription leaks**: inside a `datagrok-exec` block, the closure dies when the block finishes — subscriptions are GC'd along with it, so no leak. In viewer code, push to `viewer.subs[]` (`js-api/viewer.ts:398` "stream subscriptions - will be canceled when the viewer is detached"). The skill should note both.

## 6. Existing `grokky.filterRows` — current state and bugs (Area E)

Source (`grokky-api.ts:13-31`):

```ts
interface FilterCriteria {
  min?: number; max?: number;
  values?: string[];
  substructure?: string;       // molBlock or SMILES/SMARTS
}

function filterRows(view: DG.TableView, column: string, criteria: FilterCriteria): void {
  const fg = view.getFiltersGroup();
  if (criteria.min !== undefined || criteria.max !== undefined) {
    const col = view.dataFrame.col(column);
    fg.updateOrAdd({type: DG.FILTER_TYPE.HISTOGRAM, column,
      min: criteria.min ?? col?.min, max: criteria.max ?? col?.max} as any);
  } else if (criteria.values !== undefined)
    fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column, selected: criteria.values} as any);
  else if (criteria.substructure !== undefined)
    fg.updateOrAdd({type: DG.FILTER_TYPE.SUBSTRUCTURE, column, columnName: column, molBlock: criteria.substructure} as any);
}
```

`DG.FILTER_TYPE` (`js-api/const.ts:145-153`): `HISTOGRAM='histogram'`, `CATEGORICAL='categorical'`, `MULTI_VALUE='multi-value'`, `BOOL_COLUMNS='bool-columns'`, `FREE_TEXT='free-text'`, `COLUMN_FREE_TEXT='column-free-text'`, `SUBSTRUCTURE='Chem:substructureFilter'`.

**Bugs**:

1. **Substructure path silently matches zero rows** when input is SMILES or SMARTS. The `Chem:substructureFilter` widget reads `state.molBlock` (`chem/widgets/chem-substructure-filter.ts:476-481`) and feeds it to RDKit as a molblock. If the string isn't a molblock, RDKit returns no atoms → 0 matches with no error. Confirmed by the platform Chem package always converting before pushing the state: `chem/package.ts:1967-1976` does `convertMolNotation(molecule, Notation.Smiles, Notation.MolBlock)` when the column units are SMILES, then sets `molBlock: molblock`.

2. **Range filter with no `min`/`max`** falls into the histogram branch and uses `col?.min`/`col?.max` as both bounds — effectively no-op. Low impact (caller has to specifically omit both yet pass the criteria), but a guard `if (criteria.min === undefined && criteria.max === undefined)` is cleaner.

3. **`view.getFiltersGroup()` creates the default filter panel** as a side effect (`view.ts:404-405`, default `createDefaultFilters=true`). The Chem package consistently passes `{createDefaultFilters: false}` (`chem/package.ts:665`, `:1972`, `chem/tests/substructure-filter-tests.ts:417,426`, `chem/widgets/structural-alerts.ts:158`). Worth doing the same.

4. **No clear/invert/predicate paths.** Users say "clear the filter", "show only rows where X", "invert" — none of those map to `updateOrAdd`.

## 7. Substructure detection — recommended strategy (Area F)

Inputs we get: SMILES (`'c1ccccc1'`), SMARTS (`'[!#6]1[#6][#6][#6][#6][#6]1'`), molblock (multi-line, ends with `M  END`). The platform exposes:

- `DG.chem.isMolBlock(s)` (`js-api/chem.ts:68-70`) — pure JS check: `s.includes('M  END')`. Cheap.
- `DG.chem.isSmarts(s)` (`chem.ts:904-913`) — async-ish, calls `Chem:isSmarts` function. Requires the Chem package to be loaded.
- `DG.chem.checkSmiles(s)` (`chem.ts:896-902`) — same caveat.
- `DG.chem.convert(s, fromNotation, toNotation)` (`chem.ts:888-894`) — wraps `Chem:convertMolNotation`. Synchronous (`funcCall.callSync()`). Returns a string.
- `DG.chem.Notation` enum (`chem.ts:51-59`): `Smiles`, `CxSmiles`, `CxSmarts`, `Smarts`, `MolBlock`, `V3KMolBlock`.

How the Chem package itself routes (the canonical pattern, `chem/package.ts:1966-1977` and sketcher `chem.ts:291-300`):

```ts
if (substructure || isSmarts(s))
  setSmarts(s);
else if (isMolBlock(s))
  setMolFile(s);
else
  setSmiles(s);
```

For our wrapper, we need a molblock to pass to `FILTER_TYPE.SUBSTRUCTURE`. Two viable strategies:

**Strategy 1 — UI path (preferred, matches existing wrapper shape).** Detect → convert to molblock → `updateOrAdd`. Pseudocode:

```ts
function toMolblockForFilter(s: string): string {
  if (DG.chem.isMolBlock(s)) return s;
  // isSmarts is a function call — try SMARTS first because SMILES->molblock conversion
  // of a SMARTS string can succeed but lose meaning. The Chem package gates on isSmarts.
  if (DG.chem.isSmarts(s))  return DG.chem.convert(s, DG.chem.Notation.Smarts, DG.chem.Notation.MolBlock);
  return DG.chem.convert(s, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
}
```

Then pass the result as `molBlock` to `updateOrAdd`. The resulting filter shows up in the side panel, is editable, and persists in layouts.

**Strategy 2 — BitSet path (no UI).** Use `grok.chem.searchSubstructure(col, pattern, settings?)` (`chem.ts:766-775`) which accepts SMILES/SMARTS/molblock natively and returns a `Promise<BitSet>`. Then `df.filter.copyFrom(bs)`. Pattern from `samples/domains/chem/substructure-search.js:4-6`, `substructure-search-smarts.js:4-5`, and `substructure-search-smarts-molfile-als.js:5-28`:

```js
let bs = await grok.chem.searchSubstructure(t.col('smiles'), 'c1ccccc1');
t.filter.copyFrom(bs);
```

This is **simpler and more robust** (the chem package does all detection internally), but produces no UI widget — the user can't see or edit the filter, and combining with other filters requires the user to `and`/`copyFrom` manually.

**Recommendation**: default to Strategy 1 (UI path) to keep the contract of `filterRows` consistent (everything you do shows up in the filter panel). Expose Strategy 2 as a separate helper `filterSubstructure(view, column, query, {ui?: boolean})` with `ui` default `true`. When `ui: false`, the BitSet is written to `df.filter` directly (preferable for big tables where the filter widget is laggy).

`isSmarts` requires Chem package to be loaded. Add a guard: if `DG.chem.isSmarts` throws (Chem missing), fall back to the heuristic "treat as SMILES" with a `console.warn`.

For the molblock-failover when caller wants the molblock used directly, expose `{molBlockFailover}` (matches `searchSubstructure` settings, `chem.ts:766-768`).

## 8. Filter UI integration (Area G)

`TableView.getFiltersGroup({createDefaultFilters?})` (`view.ts:404-405`). Returns a `FilterGroup` (`js-api/viewer.ts:510-553`):

| Method | Use |
|---|---|
| `add(state)` | append a new filter |
| `updateOrAdd(state, requestFilter?)` | idempotent — update existing or add |
| `remove(filter)` | drop a filter |
| `getStates(columnName, filterType)` | inspect current state |
| `setEnabled(filter, active)` | toggle one |
| `setActive(active, notify)` | toggle the whole group |
| `setExpanded(filter, active)` | UI fold/unfold |
| `filters` (getter) | array of `Filter | Widget` |
| `getFilterSummary()` | HTML summary node |

`FilterState` is just `{type: FILTER_TYPE | string, column?: string, …filter-specific fields}` (`js-api/viewer.ts:503-506`). Per-type shapes seen in source:

| `type` | Required extra fields | Example |
|---|---|---|
| `HISTOGRAM` | `min`, `max` | range filter |
| `CATEGORICAL` | `selected: string[]` | category multi-select (note: `grokky-api.ts:28` uses `selected`, not `value`) |
| `MULTI_VALUE` | column, value parsing rules | for separator-delimited cells |
| `FREE_TEXT` / `COLUMN_FREE_TEXT` | text search; `samples/data-frame/filtering/update-filter-state.js:7-12` uses `type: 'text'` + `gridNames`+`value` for a free-text filter |
| `SUBSTRUCTURE` | `column`, `columnName`, `molBlock`, optionally `searchType`, `simCutOff`, `fp` (`chem/widgets/chem-substructure-filter.ts:476-481`) |

`view.filters(options?)` (`view.ts:491`) is **deprecated** — it's a shortcut to `addViewer(VIEWER.FILTERS, options)`. Prefer `view.getFiltersGroup()`.

There is no `view.getFiltersDataFrame()` in current source. (The audit findings mentioned it; I don't find it.) The closest API is `fg.getStates(columnName, type)` for inspecting current state, plus `fg.filters` to enumerate filter widgets.

## 9. Multi-condition idioms (Area H)

| Intent | Code |
|---|---|
| Lipinski-style AND in one pass | `df.filter.init(i => mw.get(i) < 500 && logP.get(i) < 5)` |
| AND with an already-computed BitSet (e.g. substructure) | `df.filter.and(subBs)` |
| OR | `df.filter.or(otherBs)` or fold into the predicate |
| Stack with UI filter | `fg.updateOrAdd({type: HISTOGRAM, column: 'MW', max: 500})` then `fg.updateOrAdd({type: HISTOGRAM, column: 'cLogP', max: 5})` — they collaborate via `onRowsFiltering` |
| Numeric range | `df.filter.init(i => col.get(i) >= min && col.get(i) <= max)` |

When mixing UI and programmatic, the rule is: UI filters write their contribution into `df.filter` during `onRowsFiltering`. If you then `df.filter.init(...)` outside that event, your write is overwritten on the next filter pass. **Programmatic mutations of `df.filter` only "stick" in two cases**:

1. The DF has no active UI filters → your write persists.
2. You subscribe to `onRowsFiltering` and apply your mask there.

This is the single biggest mental model gap for `df.filter` users.

## 10. "Show only" vs "drop rows" — distinct intents (Area I)

| User says | Operation | DF identity preserved? |
|---|---|---|
| "filter to X", "show only X" | mutate `df.filter` (BitSet) or `getFiltersGroup().updateOrAdd(...)` | yes; rows still in DF |
| "remove rows where X" | `df.rows.removeWhereIdx(i => !predicate(i))` (`row.ts:238`) or `df.rows.removeWhere(row => …)` (`row.ts:232`) | no; rows gone |
| "give me a DF of only the filtered rows" | `df.clone(df.filter)` (`data-frame.ts:290`) | yes; new DF returned, original untouched |
| "drop everything except selection" | `df.clone(df.selection)` (selection is a BitSet) | new DF |

`df.rows.removeWhereIdx` is destructive and irreversible. `df.clone(df.filter)` is the right shape when the user wants "the filtered rows" as a real DF (export, send to a script, etc.). The audit raises this distinction explicitly — the skill must teach it.

## 11. Subscribing to filter events (Area J)

```js
df.onFilterChanged.subscribe(_ => console.log('filter ->', df.filter.trueCount));
df.onRowsFiltered.subscribe(_ => /* mask is final */);
df.onRowsFiltering.subscribe(_ => /* mask is being computed; mutate df.filter here */);
```

In `datagrok-exec` blocks: the subscription closure dies with the block, no leak. The skill should note this matter-of-factly — Claude will sometimes ask whether to unsubscribe.

In viewer/widget code: push the subscription onto `viewer.subs[]` (`js-api/viewer.ts:398`) — it auto-disposes on detach. `DG.Widget` does the same with its `subs` array. The user memory already notes this ("DG.Widget auto-cleans subs").

`onMetadataChanged` filter example from `samples/data-frame/events/events.js:33-35` shows how to narrow with `rxjs.operators.filter`. Same shape applies to filter events when you want to react to specific column filters.

## 12. Proposed grokky helpers

Signatures only — no bodies. Rationale below each.

```ts
// 12.1 — fixed substructure path; ui flag selects FilterGroup-vs-BitSet strategy
function filterRows(view: DG.TableView, column: string, criteria: {
  min?: number;
  max?: number;
  values?: string[];
  substructure?: string;            // SMILES | SMARTS | molblock — auto-detected
  searchType?: 'substructure' | 'similar' | 'exact';  // forwards to chem filter state
  ui?: boolean;                      // default true — show in filter panel; false → BitSet write
}): Promise<void>;

// 12.2 — show every row
function clearFilter(df: DG.DataFrame): void;

// 12.3 — flip the current mask
function invertFilter(df: DG.DataFrame): void;

// 12.4 — write a fresh predicate mask in one shot; clobbers existing programmatic filter
function filterByPredicate(df: DG.DataFrame, predicate: (i: number) => boolean): void;

// 12.5 — robust substructure search; ui:true puts it in the panel, ui:false writes BitSet
async function filterSubstructure(
  view: DG.TableView,
  column: string,
  query: string,
  opts?: {ui?: boolean; molBlockFailover?: string; searchType?: string},
): Promise<void>;

// 12.6 — "give me a real DF of just the visible rows"
function filteredDf(df: DG.DataFrame, opts?: {withSelection?: boolean; cols?: string[]}): DG.DataFrame;

// 12.7 — destructive drop (distinct from filtering!)
function dropRows(df: DG.DataFrame, predicate: (i: number) => boolean): number;  // returns count removed

// 12.8 — JSON snapshot of the current filter state for debugging / reports
function describeFilter(df: DG.DataFrame): {
  rowCount: number;
  visibleCount: number;
  hiddenCount: number;
  uiFilters?: Array<{type: string; column?: string; summary?: string}>;
};
```

### Rationale

- **12.1 `filterRows` (fixed)** — backward-compatible signature plus `searchType` + `ui` flags. Internally routes through `toMolblockForFilter` (Strategy 1) when `ui !== false`. When `ui === false`, calls `filterSubstructure` under the hood. Returns `Promise<void>` because the substructure path is async (`searchSubstructure` and `isSmarts` both involve function calls). For range/categorical paths the promise resolves synchronously.
- **12.2 `clearFilter`** — pure one-liner (`df.filter.setAll(true)`). The reason it deserves a helper: Claude will reach for `setAll(false)` and silently break things. Forcing the named call kills that footgun. Also removes UI filters if a `TableView` is passed (overload — but keep first arg `DG.DataFrame` for predictability).
- **12.3 `invertFilter`** — `df.filter.invert()` plus a note about UI-filter collisions (programmatic invert gets clobbered next time a UI filter applies). Document, don't try to fight it.
- **12.4 `filterByPredicate`** — wraps `df.filter.init(predicate)`. The value is teaching: "predicate returns `true` if row passes" — backwards from the `removeWhere` predicate which is "true if row goes away". Naming it differently from `df.rows.filter` makes the difference unmissable.
- **12.5 `filterSubstructure`** — the substructure-specific helper. Robustly accepts SMILES/SMARTS/molblock, routes via `DG.chem.isMolBlock`/`isSmarts` heuristics, and falls back to either `searchSubstructure` (BitSet) or `updateOrAdd` (UI). Worth a separate helper because chem-specific knobs (`searchType: 'similar'`, `simCutOff`, `fp`, `molBlockFailover`) only apply here.
- **12.6 `filteredDf`** — sugar over `df.clone(df.filter, opts?.cols ?? null, opts?.withSelection ?? false)`. Fixes the "I want the filtered rows" intent without ambiguity.
- **12.7 `dropRows`** — explicit destructive op. Returns the count removed for confirmation. Wraps `df.rows.removeWhereIdx(i => predicate(i))`. **Note**: predicate matches the predicate's natural reading: "drop where this is true". This is the **opposite polarity** from `filterByPredicate` ("keep where this is true"), which we should call out loudly in the skill.
- **12.8 `describeFilter`** — read-only inspection. Returns counts and a compact list of UI filters (via `fg.filters.map(f => f.saveState?.())` or `fg.getFilterSummary().textContent`). Useful for "what's filtered right now?" questions.

**Helpers I considered and rejected**:

- `andFilter(df, other)` / `orFilter(df, other)` — too thin; `df.filter.and(other)` is already a one-liner.
- `getFilteredIndexes(df)` — `df.filter.getSelectedIndexes()` is already terse and discoverable via `describeFilter`.
- `onFilterChanged(df, cb)` subscription wrapper — sub leaks are a viewer-scope concern that doesn't apply in `datagrok-exec` blocks. Just document `df.onFilterChanged.subscribe(...)`.

## 13. Anti-patterns Claude must avoid

1. **`df.filter.setAll(false)` to "clear" the filter.** That hides every row. Use `df.filter.setAll(true)` or `grokky.clearFilter(df)`.
2. **Passing SMILES/SMARTS as `molBlock` to a `SUBSTRUCTURE` filter state.** The current `grokky.filterRows` does this (`grokky-api.ts:30`); it silently matches zero rows. Convert via `DG.chem.convert(...)` or use `grokky.filterSubstructure(...)`.
3. **Looping `df.rows` to build a filter.** Use `df.filter.init(i => …)` — reads columns directly and writes the buffer in one pass (`bit-set.ts:170`). The source explicitly warns against `RowList` in perf paths (`row.ts:191-192`).
4. **Confusing `df.filter` polarity.** `true` = passes. `false` = hidden. (`bit-set.ts:132` — `setAll(true)` is "all visible".)
5. **`df.clone(df.filter)` to "filter the DF".** Functional but wasteful — `clone` allocates a new DF. If the user wants the rows hidden, mutate `df.filter`. If they want a real new DF (export, downstream pipeline), then yes, clone.
6. **`df.rows.removeWhereIdx(...)` to "filter rows".** Destructive — rows are gone. Confirm with the user that they want deletion, not filtering.
7. **Writing to `df.filter` outside `onRowsFiltering` when UI filters exist.** Your write will be overwritten on the next filter cycle. Either subscribe to `onRowsFiltering`, or disable the UI filter group first (`fg.setActive(false)`).
8. **`view.filters(...)`** — deprecated (`view.ts:488-493`). Use `view.getFiltersGroup({createDefaultFilters: false}).updateOrAdd(...)`.
9. **Subscribing to `df.onFilterChanged` from inside a viewer without pushing onto `viewer.subs[]`.** Leak. (Not an issue inside `datagrok-exec` blocks — the closure is GC'd.)
10. **Forgetting `await` on `grokky.filterRows({substructure: ...})`.** The fixed version is async because conversion goes through Dart-backed `Chem:convertMolNotation`.
11. **Calling `view.getFiltersGroup()` without `{createDefaultFilters: false}`** when you only want to add one filter — it'll create default histogram/categorical filters for every column on first call (`view.ts:404-405`). Matches Chem package practice (`chem/package.ts:665`, etc.).
12. **`df.filter.init(...)` to **add** a condition on top of existing.** `init` zeros the buffer first (`bit-set.ts:175`). To intersect with current, build a separate BitSet via `DG.BitSet.create(df.rowCount, i => …)` and `df.filter.and(that)`.

## 14. Open questions

1. **Substructure detection robustness without Chem loaded.** `DG.chem.isSmarts` and `DG.chem.convert` route through `Func.find({package: 'Chem', ...})` (`chem.ts:889`, `:905`). If Chem isn't installed/loaded, they throw or return false-y. Should `filterSubstructure` detect Chem availability and fall back gracefully (e.g., "treat as SMILES, pass molblock-failover blank, rely on chem-substructure-filter widget if it loads later")? Proposal: yes, with `console.warn` and a `chemAvailable: boolean` field on the returned descriptor.
2. **`ui: true` default for `filterSubstructure`.** Confirm. The audit findings imply the user wants robust handling first, UI second. I lean `ui: true` because that matches existing `filterRows` contract (everything shows up in the panel), but for very large tables the widget is slow and `ui: false` (BitSet write) is materially better. Maybe a heuristic — if `df.rowCount > 50_000`, default `ui: false` and `console.info` the choice?
3. **`clearFilter` and UI filter group.** Should `grokky.clearFilter(df_or_view)` also disable/remove UI filters, or just clear the BitSet? Proposal: if given a `TableView`, also call `fg.setActive(false)` (preserves the widgets but turns them off). If given a `DataFrame`, just `setAll(true)`. Open: should we expose `removeAllFilters(view)` separately for "delete the widgets too"?
4. **Polarity inversion between `filterByPredicate` and `dropRows`.** Two helpers with opposite predicate polarity living in the same skill is dangerous. Options: (a) document loudly with a side-by-side table; (b) rename `dropRows` to `dropRowsWhere` to telegraph the polarity. I prefer (a) + (b).
5. **Filter event sample — viewer scope.** Should the skill include a runnable example of `onRowsFiltering` subscription, given that the right place is a viewer/widget? Or just docs? Proposal: docs only, with a pointer to the `Filter` class for users actually building a custom filter.
6. **`describeFilter` UI-filter inspection.** `fg.filters` returns `Array<Filter | Widget>` (`js-api/viewer.ts:534`). For Dart-side widgets (`Widget`), reading state may not be straightforward — `getStates(columnName, type)` is column+type scoped, not "give me everything". Worth a discovery pass: does `fg.filters` actually yield enough to round-trip? If not, fall back to text summary from `fg.getFilterSummary()` (`viewer.ts:529`).
7. **`grokky.filterRows` is now async.** Backwards compat: callers that didn't `await` get a `Promise` returned, which they may pass around (harmless for void-return). Worth flagging in the skill — and possibly in a CHANGELOG entry for Grokky.
8. **`DG.STRUCTURE_FILTER_TYPE`** (`const.ts:759-762`: `Sketch`, `Categorical`) — when is this used vs `FILTER_TYPE.SUBSTRUCTURE`? Quick read: it's the column-level tag `.structure-filter-type` (`const.ts:325`) that picks which widget renders for a molecule column. Not a filter-state field. Worth a footnote so Claude doesn't confuse them.
