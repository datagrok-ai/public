# Code Analysis — PowerGrid & CddVaultLink

Scope: source under `public/packages/PowerGrid/src/` and `public/packages/CddVaultLink/src/`.
Reviewed on 2026-04-16. Findings grouped by severity; line numbers reference current `master`.

---

## PowerGrid

### Bugs (wrong behavior)

#### 1. `tags-cell-renderer.ts` — `getColor()` returns `undefined` on first use of a tag
`src/cell-types/tags-cell-renderer.ts:50-57`

```ts
const getColor = (tag: string) => {
  const colors = gridCell.gridColumn.temp['catColors'] ??= {};
  if (colors[tag] || (colors[tag] === 0))
    return colors[tag];

  const keys = Object.keys(colors);
  colors[tag] ??= DG.Color.getCategoricalColor(keys.length);   // <-- no `return`
};
```

The "new tag" branch assigns the color into the cache but never returns it. The caller does
`g.fillStyle = DG.Color.toHtml(color)` and `DG.Color.getContrastColor(color)` with `color === undefined`,
so the very first render of each new tag ends up with an invalid fill style (silent corruption in Canvas —
the previous fill is reused). Subsequent renders hit the cache hit branch and look fine, which is why this
never surfaced as a reproducible test failure.

Fix:

```ts
return colors[tag] ??= DG.Color.getCategoricalColor(Object.keys(colors).length);
```

#### 2. `sparklines-lines.ts` — off-by-one in hit-testing when there is only one column
`src/sparklines/sparklines-lines.ts:86-87`

```ts
const activeColumn = Math.floor((mousePoint.x - b.left + Math.sqrt(minDistance)) /
  b.width * (cols.length - 1 > 0 ? cols.length - 1 : 1));
```

When `cols.length === 1`, the formula multiplies by `1`, so `activeColumn` becomes a function of the mouse
x (anything from 0..∞) instead of always 0. With `b.width` tiny or mouse past the right edge the computed
index exceeds 0 and the tooltip is suppressed. For the 1-column case `activeColumn` must be clamped to 0.

#### 3. `shared.ts` — `getScaledNumber` divides by zero / yields NaN for Row normalization on single-column rows
`src/sparklines/shared.ts:181-191`

When `normalization === Row` and only one column is passed (or the row has one non-null value),
`Math.min(...rowValues)` === `Math.max(...rowValues)`, which the inline `normalize` short-circuits to `0`
only when `min === max` strictly. If a `NaN` sneaks in through `col.getNumber(row)` on an integer
null cell (`DG.INT_NULL`), `scaleValue(NaN)` flows through and poisons the rest of the calc. The caller
filters `isNone` in the renderers but `scaleSettings`/`getScaledNumber` still reads *all* cols in `cols` for
the row, so one missing-value column drags the whole row to NaN.

Add `if (!isFinite(v)) continue;` in the `for (const col of cols)` loop, or filter `isNone(row)` columns
before calling `getScaledNumber`.

#### 4. `piechart.ts` — does not `filter((c) => c != null)` in hit-test
`src/sparklines/piechart.ts:118`

```ts
const cols = gridCell.grid.dataFrame.columns.byNames(settings.columnNames);
```

All the other renderers (and the render path in the same file on line 224) filter `null` from
`byNames(...)` to handle deleted source columns. `onHit` doesn't, so hovering a sparkline whose source
column has just been deleted throws (`cols[activeColumn].getNumber(row)` on undefined). The autostart
subscription in `package.ts` is supposed to prune missing names from `columnNames`, but that only runs
on `onColumnsRemoved` — a column that is *not* present at hover time (e.g. a typo in a saved project
JSON, or a row renamed before the listener was wired) still reaches this path.

#### 5. `piechart.ts` — divide-by-zero in Angle mode when all values are 0 or missing
`src/sparklines/piechart.ts:279-295`

```ts
const sum = getColumnsSum(cols, row);
...
const endAngle = currentAngle + 2 * Math.PI * cols[i].getNumber(row) / sum;
```

If every value on the row is 0 or `isNone`, `sum === 0` and every `endAngle` is `NaN`. The canvas silently
fails to draw but the loop still calls `g.arc(..., NaN, NaN)`, which on some browsers flags state
corruption for subsequent draws. Guard with `if (sum === 0) return;`.

#### 6. `radar-chart.ts` — hit-test uses `gridCell.bounds.midX/midY` but render uses a `fitSquare()` box
`src/sparklines/radar-chart.ts:40-56`

The hit test computes `vectorX/Y` from `gridCell.bounds.midX/midY` and feeds them into
`getAxesPointCalculator(cols, box)` where `box` is `bounds.fitSquare().inflate(-2,-2)`. When the cell is
not square, the rendered polygon is centered in the square sub-box but the hit test is centered on the
full rectangle. Tooltips light up off-center from the visible marker. Use `box.midX`/`box.midY` (or keep a
single rect for both paths).

#### 7. `radar-chart.ts` — `onHit` wrap-around block is correct but needlessly duplicated
`src/sparklines/radar-chart.ts:51-53`

```ts
activeColumn = activeColumn > cols.length - 1 ? 0 : activeColumn;
valueForColumn = Math.floor(valueForColumn + maxAngleDistance) > cols.length - 1 ?
  cols.length - valueForColumn : valueForColumn;
```

**Not a bug** (earlier revision of this doc claimed otherwise — retracted). `cols.length - valueForColumn`
when `valueForColumn` is near `cols.length` yields a small *positive* number — the wrap-around distance
to axis 0 — and `Math.abs()` on line 62 makes the sign irrelevant. Hit tests near axis-0 do fire from
both sides.

Still a readability smell: `Math.floor(valueForColumn + maxAngleDistance)` is computed twice and the two
ternaries branch on the same condition but mutate different variables. Collapses to one `if`:

```ts
if (activeColumn > cols.length - 1) {
  activeColumn = 0;
  valueForColumn = cols.length - valueForColumn;
}
```

#### 8. `forms.ts` — `onMouseMove` schedules `setTimeout` without a delay and without clearing on leave
`src/forms/forms.ts:227-233`

```ts
onMouseMove(gridCell: DG.GridCell, e: MouseEvent) {
  const el = scene?.hitTest(e.offsetX, e.offsetY);
  if (el?.style?.tooltip)
    setTimeout(() => ui.tooltip.show(el.style!.tooltip!, e.x + 20, e.y - 20));
  else
    ui.tooltip.hide();
}
```

Every mouse move over the form schedules a new `setTimeout(...)` with no delay, and nothing cancels the
pending ones. Moving the mouse across the cell stacks up N `tooltip.show` calls that fire on the next
microtask, causing the tooltip to flicker. Use `setTimeout(fn, 200)` *and* store the handle so the
next `onMouseMove`/`onMouseLeave` can `clearTimeout`. The `onMouseEnter` in the same file already uses
`200ms`; adopt the same pattern here.

#### 9. `forms.ts` — `scene` is a module-level `let`, shared across every cell of every grid
`src/forms/forms.ts:48` + `219-228`

```ts
let scene: Scene;
...
onMouseEnter(...) { scene = FormCellRenderer.makeScene(gridCell); ... }
onMouseMove(gridCell, e) { const el = scene?.hitTest(e.offsetX, e.offsetY); ... }
```

A single module-level `scene` is written by `onMouseEnter` and read by `onMouseMove`. If two grids are
visible (e.g. a side-by-side layout) hovering one form updates `scene` and then stale `hitTest`s happen
against the wrong grid's coordinate space until the next enter event. Make it a `WeakMap<DG.Grid, Scene>`
keyed on the grid or on the cell column — or simply recompute in `onMouseMove` when bounds change.

#### 10. `image-cell-renderer.ts` — `getImage` starts loading, returns immediately; the LRU caches an unresolved image
`src/cell-types/image-cell-renderer.ts:27-45` + `56`

```ts
const img: Promise<HTMLImageElement> = ImageCellRenderer.images.getOrCreate(url, () => this.getImage(gridCell));
```

`getImage` creates `new Image()`, sets `src`, and returns it synchronously wrapped in a `Promise`. That
Promise resolves instantly (with the not-yet-loaded image) because `getImage` is `async` but `await`s
nothing on the HTTP branch (only on the `.exists`/`files.list` branch). The renderer then checks
`img.complete` before drawing, so correctness is fine — but the cache entry is an already-resolved
promise for an image whose `onload` fires later and calls `gridCell.grid.invalidate()`. On the second
render, the cache hit bypasses a fresh `invalidate` wiring. For URLs that fail to load
(404, CORS denied), the cache permanently holds a broken entry and the image never retries.

Two-part fix: (a) only cache on successful load (`image.onload = () => { cache.set(url, image); ... }`);
(b) cache the image element, not the promise — promises-in-LRU is a code smell here.

#### 11. `bar-cell-renderer.ts` — `BarCellSettingsProperties` declares `radius` as STRING
`src/cell-types/bar-cell-renderer.ts:11`

```ts
const BarCellSettingsProperties = {
  color: DG.Property.js('color', DG.TYPE.STRING),
  radius: DG.Property.js('radius', DG.TYPE.STRING),   // <-- should be FLOAT/INT
};
```

The interface declares `radius: number` but the property descriptor says STRING. Dead code today
(the value isn't used by `render`), but if anyone wires up the settings panel it will render as a text
input. Also: `BarCellSettingsProperties` is defined and never exported or referenced anywhere else —
candidate for deletion.

#### 12. `scatter-plot.ts` — re-enters `document.body.appendChild(plot.root)` on every render
`src/sparklines/scatter-plot.ts:52-80`

The render path appends the plot's root to `document.body`, resizes the canvas, then removes it again.
`plot.render(g)` is commented out ("version conflict"), so this does nothing visible but still thrashes
the DOM on every cell repaint. Either remove the scratch-draw block entirely (the whole renderer is
effectively a no-op) or leave a single-line comment saying so.

### Code smells / maintainability

#### 13. `piechart.ts` — `let minRadius` at module scope; mutated inside `render`
`src/sparklines/piechart.ts:19` + `226`

```ts
let minRadius: number;
...
render(...) {
  ...
  minRadius = Math.min(box.width, box.height) / 10;   // writes shared module state
  ...
}
```

Any other caller reading `minRadius` between two `render` invocations sees a value dependent on whichever
cell was painted last. Only `onHit` reads it — which happens in between renders, off the most-recent
value — so it's a de-facto global. Move to a local in `render` and pass it into the `renderSubsector`
call; `onHit` already computes its own radius.

#### 14. `shared.ts` — `getRenderColor` is a deeply nested ternary
`src/sparklines/shared.ts:345-351`

One expression, four branches, one line wrap. Replace with a `switch` on `settings.colorCode` or an early-
return chain. It's currently the hardest-to-read function in the file.

#### 15. `shared.ts` — `createBaseInputs` inner `getMinMaxProperties` / `getAdditionalColumns` could be flattened
`src/sparklines/shared.ts:255-308`

Three nested declaration-only helpers that are each called exactly once. Inlining them removes ~40 lines
and an indentation level.

#### 16. Backwards-compat scattered across `getSettings` in every sparkline
`sparklines-lines.ts:50-53`, `bar-chart.ts:27-29`, (similar in `piechart.ts`, `radar-chart.ts`)

Every renderer repeats:

```ts
//@ts-ignore: convert old format to new - backwards compatibility
if (settings.globalScale !== undefined && settings.globalScale !== null)
  settings.normalization = settings.globalScale ? NormalizationType.Global : NormalizationType.Column;
```

Pull into `shared.ts#migrateLegacySettings(settings)` and call once in `getSettingsBase` so a future
breaking change touches one place, not four.

#### 17. `sparklines-lines.ts` — `interface getPosConstants` / `function getPos` / `getPosConstants` variable all share the same name
`src/sparklines/sparklines-lines.ts:23-40` + `79-83`

Type, function, and local variable are all `getPos*`. The interface is `camelCase` when TypeScript
convention is `PascalCase`. Rename to `GetPosCtx` (interface) and `ctx` (local).

#### 18. `sparklines-lines.ts` — `renderSettings` recomputes `getSettings` via a five-line isSummary/settings dance
`src/sparklines/sparklines-lines.ts:167-169`

```ts
const settings: SparklineSettings = isSummarySettingsBase(gridColumn.settings) ? gridColumn.settings :
  gridColumn.settings[SparklineType.Sparkline] ??= getSettings(gridColumn);
```

The exact same shape appears in `bar-chart.ts`, `piechart.ts`, `radar-chart.ts`, `forms.ts`. Put it in
one helper (`resolveSettings<T>(gc, type, getter)`) — currently five copies of the same conditional.

#### 19. `radar-chart.ts` — hand-rolled `it.range()`
`src/sparklines/radar-chart.ts:10-12`

```ts
class it {
  static range = (n: number) => [...Array(n).keys()];
}
```

`it.range(cols.length).map(...)` is a way to write `for (let i = 0; i < cols.length; i++)`.  Replace with
a regular `for` loop (every other renderer uses loops) and delete the `class it` shim. `wu` is also
already imported in this package and has a more idiomatic `wu.count().take(n)` if a range iterator is
really wanted.

#### 20. `radar-chart.ts` — path built via `.map`, then rendered via chainable `setFillStyle().polygon().fill()`
`src/sparklines/radar-chart.ts:115-123`

Fine stylistically, but note `canvas-extensions.ts`'s `polygon()` stubs are patched onto
`CanvasRenderingContext2D.prototype` as a side effect of importing this file indirectly. Any consumer
that imports PowerGrid's code (via tests, via TypeScript typings, etc.) inherits these mutations.
Document this coupling in `canvas-extensions.ts` or move to a non-prototype helper `polygon(g, pts)`.

#### 21. `package.ts` autostart subscribes to events but calls `grid.detach()` from its own detach handler
`src/package.ts:368`

```ts
const gridDetachedSub = grid.onDetached.subscribe(() => grid.detach());
```

Detaching a grid inside its own `onDetached` handler is either a typo (meant to call
`colsRemovedSub.unsubscribe()` / `colsRenamedSub.unsubscribe()`) or a re-entrancy bug that hasn't blown
up because `detach()` is idempotent on already-detached grids. Either way the intent is unclear. If the
goal is to clean up subscriptions, subscribe-and-push into the grid's own subscription bag via
`grid.sub(sub)` (which is already done two lines below) — no manual `onDetached` needed.

#### 22. `package.ts` — indentation inconsistency in the `PackageFunctions` class
`src/package.ts:132-145` + `159-161`

Two decorators are indented with 4 spaces while the rest use 2. ESLint config is 2-space but these lines
slipped through (they sit inside `@grok.decorators.func(...)` blocks that aren't picked up by the
`no-mixed-spaces-and-tabs` rule). Fix indentation.

#### 23. `package.ts:344` — `summaryCols[summaryCols.length] = gridCol` instead of `.push(gridCol)`

```ts
summaryCols[summaryCols.length] = gridCol;
```

Same behavior as `.push()`, slower in V8 if the engine can't specialize, and harder to read. Use `push`.

#### 24. `forms.ts:246` — `value: settings.showColumnNames ?? 'Auto'` typed against a wider `string[]`

The choice input accepts `string[]` but the enum has three string-literal members. The cast
`value as ColumnNamesVisibility` papers over it. Define `items: ColumnNamesVisibility[]` and let
`ui.input.choice<ColumnNamesVisibility>(...)` carry the type.

#### 25. `confidence-interval-cell-renderer.ts` — dead code
`src/cell-types/confidence-interval-cell-renderer.ts:6` + `320-329` + `261-304`

`subscribedGrids` is declared but unused; `renderHeader` / `getHeaderScale` are defined but never called
(the only call site is commented out inside `render`). Either wire them up or delete — ~60 lines.

#### 26. `forms.ts:152-153` — `setTimeout(fn, 200)` vs `setTimeout(fn)` inconsistency vs flicker fix

See #8. Same file has the correct 200ms delay in `onMouseEnter` (line 222) but not in `onMouseMove`.

#### 27. `canvas-extensions.ts:19` — spurious semicolon after `declare global { ... };`

```ts
declare global {
  ...
};   // <-- block statement with stray semicolon
```

Not a compile error but ESLint's `@typescript-eslint/no-extra-semi` will catch it.

#### 28. Duplicate `getSparklinesContextPanel` copy in every renderer

Each sparkline renderer declares `async getContextValue(gridCell: DG.GridCell): Promise<any>`. These are
identical modulo the `getSettings` call. Move to a shared mixin / base class, or bake into
`getSparklinesContextPanel(gridCell, cellType)`.

#### 29. `multi-choice-cell-renderer.ts:50` — `DG.Color.lightGray` vs `'lightgray'` vs `'#D0D0D0'` inconsistency

The codebase uses three different spellings for the same gray:
- `DG.Color.lightGray` (`multi-choice-cell-renderer.ts`, `bar-chart.ts`, `radar-chart.ts`)
- `'lightgrey'` (`sparklines-lines.ts`, `scene.ts`)
- `'#D0D0D0'` / `'#b9b9b9'` / `'#dcdcdc'` (`radar-chart.ts`, `stars-cell-renderer.ts`)

Per `CLAUDE.md` > "Datagrok CSS Conventions", prefer the design tokens / `DG.Color.*` constants.

#### 30. `stars-cell-renderer.ts:60` — redundant `!gridCell.cell.isNone()` check

```ts
if (starSize < MIN_STAR_SIZE || (!gridCell.cell.isNone() && (value < 0 || value > maxStars))) {
```

The early-return on line 54 already handles `isNone`. The `!isNone` inside the `||` is always true by the
time we reach line 60 — drop it.

#### 31. `hyperlink-cell-renderer.ts` — prints `link <url>` as plain text
`src/cell-types/hyperlink-cell-renderer.ts:21`

Renders `"link " + value` rather than the value itself styled as a link (no underline, no blue color,
no click handler). This is probably still a placeholder but it's been shipped — either style it properly
(blue text + `onClick` opening the URL) or mark as `@deprecated`.

#### 32. `svg-cell-renderer.ts:47` — cache eagerly stores `null` on first render, but retries happen
`src/cell-types/svg-cell-renderer.ts:40-54`

The cache flow is: on first render set `null`, kick off load, on load set the img. But the `has(cacheKey)`
+ `get(cacheKey)` dance skips reloading the image on subsequent renders even if the old `set(cacheKey, null)`
is still in place (for an SVG that is still decoding). So you can flash-of-no-content for a wide cell if
the user scrolls fast. Minor. More importantly: the cache key is `length + first 100 chars + last 50
chars`, which is not injective — two SVGs with the same prefix and suffix but different middles collide.
For a 200-char SVG cap that's fine; for long but near-identical SVGs (e.g. differently-colored variants)
it is not. Hash the string (djb2 / FNV-1a) in the key instead.

### Duplication / refactor targets (PowerGrid)

| Duplicated pattern                      | Locations                                                                     | Suggestion                                                                                                              |
|-----------------------------------------|-------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------|
| `getSettings` resolve-or-migrate shape  | all four sparkline renderers, `forms.ts`                                      | Extract `resolveSettings<T>(gc, type, defaults)` into `shared.ts`.                                                      |
| `cols = byNames(...).filter(c => c != null)`  | all four sparklines, `forms.ts`                                         | Put into `shared.ts#getVisibleColumns(gc, settings)`.                                                                   |
| `w < 20 \|\| h < 10` size gate          | sparklines-lines, bar-chart, piechart, radar-chart, image, binary-image       | Single `shouldRender(w, h, min={20,10})` helper or let each renderer define `minWidth`/`minHeight`.                     |
| `hasContextValue / getContextValue`     | all four sparklines                                                           | Base class (`SparklineCellRendererBase`) that every concrete renderer extends — 4 nearly identical implementations.     |
| `onMouseMove` tooltip show/hide         | all four sparklines                                                           | Same base class — the only thing that varies is `onHit`.                                                                |
| Legacy `globalScale` → `normalization`  | all four sparklines                                                           | `migrateLegacy(settings)` in `shared.ts`, called once inside `getSettingsBase`.                                         |

A `SparklineCellRendererBase<TSettings extends SummarySettingsBase>` with abstract `onHit` and
`renderInner` would remove ~150 lines of copy-paste across the sparkline folder.

---

## CddVaultLink

### Bugs

#### 1. `cdd-vault-api.ts` — `ApiResponse<T>` conflates two wire shapes
`src/cdd-vault-api.ts:9-18`

```ts
export interface ApiResponse<T> {
  count?: number;
  offset?: number;
  page_size?: number;
  objects?: T[];
  data?: T;
  ...
}
```

The CDD API returns two shapes: a bare JSON (`[{...}]`) for list endpoints and an envelope (`{count,
objects}`) for paged ones. `request()` always returns `{data}` wrapping the whole parsed JSON into the
`data` field — so for `queryVaults()` `response.data` is `Vault[]` (and the top-level `count/offset/objects`
properties on `ApiResponse` are *never* populated). The `objects` and `count` fields on `ApiResponse` are
dead — they look like a paging layer that never activated. Callers that write `res.data?.objects` get it
right (for paged) because `data = MoleculesQueryResult = { count, objects }`. Remove `count`, `offset`,
`page_size`, and `objects` from `ApiResponse` — they lie about the shape.

#### 2. `cdd-vault-api.ts` — mutates caller-provided params in `*Async` wrappers
`src/cdd-vault-api.ts:466-471` (and analogues in `queryCollectionsAsync`, `queryBatchesAsync`,
`queryMoleculesAsync`, `queryReadoutRowsAsync`)

```ts
export async function queryProtocolsAsync(vaultId: number, params: ProtocolQueryParams): Promise<...> {
  params.async = true;                                // <-- mutates caller's object
  const paramsStr = paramsStringFromObj(params);
  return request('GET', `/api/v1/vaults/${vaultId}/protocols${paramsStr}`);
}
```

Each `*Async` wrapper writes `params.async = true` into the caller's object. If the caller (e.g.
`getVaultStats` in `package.ts`) reuses the `{only_ids: true}` literal for multiple calls, the shared
object will carry `async: true` into later calls. Today each call creates a fresh literal so the bug is
latent, but the moment someone hoists the options object for cache-by-reference this bites. Clone:
`const q = {...params, async: true};`.

#### 3. `utils.ts` — `adjustIdColumnWidth` uses a fixed `DG.delay(100)` as a "render finished" proxy
`src/utils/utils.ts:784-789`

```ts
async function adjustIdColumnWidth(tv: DG.TableView) {
  await DG.delay(100);
  const idCol = tv.grid.col('id');
  if (idCol)
    idCol.width = 100;
}
```

A blind 100ms sleep waits for the grid to finish initial layout. On a slow render this is too short;
faster machines waste 100ms per tab switch. Subscribe to `tv.grid.onAfterRender` / the table view's
`onTableAttached` instead, or query `tv.grid.col('id')` inside a `requestAnimationFrame` loop with a cap.

#### 4. `utils.ts` — `attachLoadAllRibbon` cancellation checks the wrong token
`src/utils/utils.ts:495-530`

```ts
const ribbonViewToken = tv;           // line 497 — captured at call time (`createCDDTableView`)
...
if (openedView !== ribbonViewToken) { ... }   // line 509
```

`ribbonViewToken` is the `TableView` object, but `openedView` can be reassigned to a *different* TableView
by a sibling operation (e.g. another ribbon's "Load all" finishing after the user navigated). The check
compares TableView references, not the `viewToken` sentinel used correctly inside `createCDDTableView`
(line 471). If the user goes Molecules → Collections → back-to-Molecules before the async resolves, the
second Molecules view *is* a new TableView and the ribbon from the first one silently overwrites the
live view's dataframe. Use a local monotonic counter or compare against the module-level token created
per-open.

#### 5. `utils.ts` — `createNestedCDDNode` leaks the `openedView` pointer on error
`src/utils/utils.ts:385-389`

```ts
if (openedView) {
  ui.setUpdateIndicator(openedView.root, false);
  ui.empty(openedView.root);
  openedView.root.append(ui.divText(`Error`));
}
```

After an error, the view shows "Error" but `openedView` is still the preview view. The next tab click
closes it and replaces it, which is fine, but there's no `openedView = null` reset, so
`attachLoadAllRibbon`'s `openedView !== ribbonViewToken` check (bug #4) continues to compare against the
now-stale errored view. Set `openedView = null` after displaying the error UI, or better: keep the view
and reset the preview state.

#### 6. `utils.ts:655-665` — `onPathClick` subscription that fires arbitrarily many times
`src/utils/utils.ts:650-680`

`setBreadcrumbsInViewName` is called on every tab change and *each call* subscribes to
`breadcrumbs.onPathClick`. Nothing unsubscribes. If a user clicks between tabs 20 times and then clicks a
breadcrumb, the click handler fires 20 times, each trying to set `tree.currentItem`. Takes care of itself
since `currentItem =` is idempotent, but it does make the UI react 20× as much as needed.

Fix: store the subscription inside `openedView.onClose` (or pass the subscription into `usedView.temp`
and unsubscribe in a `usedView.closed` hook).

#### 7. `utils.ts` — `createObjectViewer` mixes `isDictionary` checks with direct property access
`src/utils/utils.ts:213-219`

```ts
function getPaneName(item: any, index: number, parentName: string): string {
  if (item.name) return item.name;          // <-- throws if item is null/undefined/primitive
  if (item.id) return item.id.toString();
```

Callers (line 232) already filter `isDictionary(item) || Array.isArray(item)` before calling, so in
practice `item` is always non-primitive. But the helper is exported-able and the contract is unclear —
at minimum, put a `if (!item || typeof item !== 'object') return 'Item ${index}';` guard at the top.

#### 8. `utils.ts` — `createCDDContextPanel` assumes `molecule_fields`/`stoichiometry`/`udfs`/`batch_fields` are always objects
`src/utils/utils.ts:708-711`

```ts
} else if (key === 'molecule_fields' || key === 'udfs' || key === 'stoichiometry' || key === 'batch_fields') {
  const acc = ui.accordion(key);
  acc.addPane(key, () => ui.tableFromMap((obj as any)[key] as any, true));
```

If the key exists but has a primitive/null value (which CDD *does* return for some fields),
`ui.tableFromMap(null, true)` will error. Guard with `typeof obj[key] === 'object' && obj[key] !== null`.

#### 9. `search-function-editor.ts` — `reset()` writes to `DG.InputBase` instances that were replaced by `init()`
`src/search-function-editor.ts:95-104` + `167`

`init()` rebuilds `this.protocolListChoice` with fresh items once protocols are fetched. `reset()` does
`this.protocolListChoice.value = ''` — but between construction and `init()` resolving, there's a window
where `this.protocolListChoice` is the empty stub. `reset()` awaits `initComplete`, which protects that
case, but the broader pattern (re-assigning inputs inside `init()`) is fragile — any other caller that
holds a reference to the pre-init input is stale. Pre-build with *all* items empty, then call
`.addOptions()` or reassign `items` after the fetch. Don't recreate the input.

#### 10. `package.ts:219` — hardcodes `vaults[0]` for context panel without user override
`src/package.ts:217-219`

```ts
//looking for molecule in the first vault
const vaultId = vaults[0].id;
const cddMols = await queryMolecules(vaultId, {structure: molecule, structure_search_type: 'exact'});
```

The context panel only checks the first vault. A user with multiple vaults will never see matches in
vault #2. The UI doesn't even indicate *which* vault was queried. At minimum, the panel should query all
vaults and aggregate, or expose a dropdown.

#### 11. `package.ts:533-576` — `cDDVaultSearch2` is an 80-line param-copying boilerplate that could be 5 lines
`src/package.ts:474-590`

```ts
if (molecules) params.molecules = molecules;
if (names) params.names = names;
params.include_original_structures = include_original_structures;
params.only_ids = only_ids;
... (26 more lines of the same)
```

All of this can be replaced by:

```ts
const params: MoleculeQueryParams = {};
for (const [k, v] of Object.entries({molecules, names, created_before, ...})) {
  if (v !== undefined && v !== null && v !== '')
    (params as any)[k] = v;
}
```

…or, given the explicit allow-list is actually useful for typing, keep individual assignments but drop
the `if` guard and let `paramsStringFromObj` filter undefineds (it already does, line 185 of `utils.ts`).

Also note the function is named `cDDVaultSearch2` with no deprecation marker on `cDDVaultSearch`. Is `2`
the keeper or the sandbox? The `@grok.decorators.func` annotation on it (`'name': 'CDD Vault search 2'`)
leaks an internal numbering into the UI — rename.

### Code smells / maintainability

#### 12. `utils.ts:168-180` — `reorderColumns` doesn't handle aliased names

The `firstColumns` array is `['id', 'name', 'smiles']`. After `createLinksFromIds`/`createMoleculeIdLinks`
replaces the `id` column, the replacement still has name `'id'`, so reordering works. But for
`createBatchesDfFromObjects` the link column becomes `molecule_id`, not `id`, and `reorderColumns` won't
pull it to the front. Consider: `const priority = ['id', 'molecule_id', 'name', 'smiles'];`.

#### 13. `utils.ts:784` — `adjustIdColumnWidth` mutates grid state after the function resolved; no way to tell it to stop

Related to #3 and #4 above. If a user switches tabs during the 100ms delay, the function still fires
`idCol.width = 100` on whichever grid now happens to have an `id` column. Same `viewToken` guard used for
`createCDDTableView` should apply here.

#### 14. `package.ts` — `CDDVaultSearchEditor` is declared as `@grok.decorators.editor()` but the TODO in the code says "is not used at the moment"
`src/package.ts:233-251`

Dead code still wired to a `'editor'` meta tag. Either hook it up (it's a reasonable customization
surface) or delete both the function and the meta reference on `cDDVaultSearch` and `cDDVaultSearchAsync`.

#### 15. `package.ts:127-128` — `createCDDTableView` called with a 9-arg positional signature

```ts
createCDDTableView([PROTOCOLS_TAB, item.name], 'Waiting for molecules',
  'CDDVaultLink:cDDVaultSearch', {...protocolSearchParams, page_size: PREVIEW_ROW_NUM},
  'CDDVaultLink:cDDVaultSearchAsync', protocolSearchParams,
  vault, treeNode);
```

9 positional args, including two pairs of `(name, params)`, with a boolean (`addFilters`) at the end.
This exact signature is repeated 5 times with subtle variations. A small `{sync: {name, params}, async:
{name, params} | null, ...}` options object would make each call site self-documenting and remove an
entire class of "which vault went where?" bugs.

#### 16. `utils.ts:14` — hardcoded path `'apps/Cddvaultlink'`

The app is registered as `CDD Vault` with decorator `'name': 'CDD Vault'` in `package.ts`. The URL path
derives from the decorator name lowercased/kebabed. Hardcoding `'apps/Cddvaultlink'` here couples this
module to the decorator in a way that breaks silently if either is renamed. Prefer deriving via
`_package.url` or a shared constant exported from `package.ts`.

#### 17. `utils.ts:316-324` — `handleInitialURL` parses the path manually

```ts
const currentTabs = url.pathname.includes(`${CDD_VAULT_APP_PATH}`) ?
  url.pathname.replace(`${CDD_VAULT_APP_PATH}`, ``).replace(/^\/+/, '')
    .split('/').map((it) => decodeURIComponent(it)) : [];
```

`URL.pathname` decoded against `CDD_VAULT_APP_PATH` that has a lowercased app name. Works today;
brittle because `createPath` (line 430) encodes as `encodeURIComponent` only the *vault*/tab parts.
Build a tiny `{vault, view, subView} = parseCddPath(url)` function and dual-test it against
`createPath()` output.

#### 18. `utils.ts:335` — `currentView ? vault.expanded = true : vault.root.click();`

Using a ternary for side effects is lint-hostile (`no-unused-expressions`). Convert to `if/else`.

#### 19. `utils.ts:376` — `let items: any[] | null = _initialItems;`

`_initialItems` arg is always passed as `null` by every caller. Dead parameter — remove it (the four
`createNestedCDDNode` call sites in `package.ts` all pass `null`/typed-null).

#### 20. `utils.ts:28-111` — three closely-related export helpers (`getAsyncResults`, `getAsyncResultsAsDf`, `runAsyncExport`, `runAsyncExportAsDf`) that CLAUDE.md says callers should not mix

The CLAUDE.md note says "do not call those three primitives directly from new code unless you need
fine-grained control." Make `getAsyncResults` / `getAsyncResultsAsDf` / `getExportId` **non-exported**
module functions (they already aren't used outside this file). Keeping them exported invites the exact
anti-pattern the docs warn against.

#### 21. `cdd-vault-api.ts:7` — `let apiKey = ''` at module scope

Module-level mutable state means two packages that load this file get independent `apiKey`s — not a
real concern under webpack's single-instance guarantee, but worth wrapping into a `getApiKey()` helper
so the lazy-init + 401-invalidate logic can be unit tested.

#### 22. `search-function-editor.ts:7` — `LAST_SEARCH_KEY = 'CDDVaultLink.lastSearch'` string

Matches `CLAUDE.md` convention, fine. But combined with the subkey `String(this.vaultId)` this is the
only place in the package using `grok.userSettings`. Document the shape in `cdd-vault-api.ts` alongside
`CDD_HOST` or move to a dedicated `preferences.ts` so future persistent preferences share a naming scheme.

#### 23. `package.ts:474-590` — repeated `@grok.decorators.param({...})` with near-identical `'options': {'category': '...', 'nullable': true, 'description': '...'}`

Over 20 parameters, most sharing `nullable: true` and a category. A couple of ideas:
- Extract a `const nullableStringParam = (category, description) => ({...})` helper. `@grok.decorators`
  accepts plain objects, so this is purely syntactic.
- Split into sub-functions by `category` — the 7 batch-field params, 6 structure params, etc. belong
  together.

#### 24. `package.ts:525-576` — identical construction pattern as `runSearch` helper above it

`runSearch` (lines 32-78) already parameterizes the sync/async switch. `cDDVaultSearch2` is a second
entry point that doesn't go through `runSearch` and duplicates the `queryMolecules → error → DataFrame`
plumbing. Decide which is the canonical path and delete or merge the other.

### Duplication / refactor targets (CddVaultLink)

| Duplicated pattern                   | Locations                                             | Suggestion                                                                                      |
|--------------------------------------|-------------------------------------------------------|-------------------------------------------------------------------------------------------------|
| `queryXxxAsync` boilerplate          | `queryProtocolsAsync`, `queryCollectionsAsync`, etc.  | `const asyncQuery = (vaultId, path, params) => request('GET', `${path}${paramsStringFromObj({...params, async: true})}`);` |
| `grok.shell.error(e?.message ?? e)`  | 8+ sites                                              | `toast(e)` helper.                                                                              |
| `getCol` / safe `byName`             | utils.ts + forms/confidence renderers in PowerGrid    | Put into `@datagrok-libraries/utils`.                                                           |
| Manual `paramsStringFromObj`         | 5 Async wrappers                                      | Pass `params` directly to `request` with a `GET` body-as-query variant.                         |

---

## Summary

**PowerGrid** — the core visual components (sparklines, forms) work but share ~150 lines of copy-paste
boilerplate across four nearly-identical renderer files that should have a base class. Real bugs
concentrated in `tags-cell-renderer` (silent color-cache bug), `piechart`/`radar-chart` hit-testing,
and the mutable module-level `scene` in `forms.ts`.

**CddVaultLink** — cleaner architecture (one entry point, well-documented CLAUDE.md), but the
`ApiResponse<T>` type is lying about the wire shape and a few spots (`*Async` wrappers, `cDDVaultSearch2`,
`adjustIdColumnWidth`) have patterns that would benefit from either deletion or extraction. The
view-token cancellation logic in `utils.ts#attachLoadAllRibbon` has a subtle race that should be fixed
before multi-tab scenarios get more complex.

Highest-leverage refactors (both packages):
1. Extract `SparklineCellRendererBase<T>` in PowerGrid — removes ~150 lines and centralizes the
   `getSettings`/`onHit`/`onMouseMove`/`getContextValue`/`renderSettings` quintet.
2. Clean up `ApiResponse<T>` in CddVaultLink — the shape lie costs readers 10+ minutes of tracing.
3. Replace the `setTimeout`-without-handle patterns in `forms.ts` and the `DG.delay(100)` in
   `utils.ts`.
