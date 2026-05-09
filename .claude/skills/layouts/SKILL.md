---
name: layouts
description: Save, restore, ship, and reapply Datagrok TableView layouts from a package using the JS API
---

# layouts

## When to use

You are scripting layout reuse from a package or app — capturing the
current view's viewer arrangement to ship with a dashboard, restoring a
`.layout` file alongside a CSV, persisting per-campaign layouts on the
gallery, or auto-applying a saved layout to a fresh table. Triggers:
"ship a layout with this dataset", "save the dashboard to the gallery",
"restore the user's last layout".

## Prerequisites

- A package scaffold (`grok create <Name>`); commands run from the
  package root, code lives under `src/`.
- `datagrok-api` imports — the article omits these:
  ```typescript
  import * as grok from 'datagrok-api/grok';
  import * as DG   from 'datagrok-api/dg';
  ```
- A `DG.TableView` (`grok.shell.addTableView(df)` or `grok.shell.tv`).
  Layout APIs live on the generic `View` base but the docstring is
  explicit: "Only applicable to certain views, such as `TableView`"
  (`DG-FACT-218`). Calling on a non-table `ViewBase` is a silent no-op.
- For server-gallery persistence: an authenticated `grok.dapi`
  session (default in any package running inside the platform).

## Steps

1. **Capture the current view's layout — TableView only.**
   `view.saveLayout()` returns a `DG.ViewLayout`. Article omits the
   optional flag `saveWithData` (default `false`); pass `true` only
   when you intend to ship a self-contained `.layout` carrying the
   underlying frame inline (`DG-FACT-218`, `DG-FACT-DRIFT-086`):
   ```typescript
   const view = grok.shell.addTableView(grok.data.demo.demog());
   const layout = view.saveLayout();                       // portable; data lives elsewhere
   const bundled = view.saveLayout({saveWithData: true}); // ships table data inline
   ```
   Expected: `layout` is a non-null `DG.ViewLayout`. Guard with
   `if (view instanceof DG.TableView)` if the caller's `view` may be
   a non-table `ViewBase`.

2. **Reapply / roll back via `loadLayout`.**
   `loadLayout(layout, pickupColumnTags?)` accepts a second flag the
   article doesn't mention (default `false`); set `true` only when
   the layout's stored column tags should be copied onto the
   destination frame (`DG-FACT-218`, `DG-FACT-DRIFT-086`):
   ```typescript
   const layout = view.saveLayout();
   view.addViewer('Histogram', {value: 'age'});
   view.grid.columns.rowHeader.width = 100;
   view.loadLayout(layout);             // rolls back both edits
   view.loadLayout(layout, true);       // also overwrites column tags
   ```
   Expected: viewer added in line 2 disappears; grid header width
   reverts. To clear without a saved layout: `view.resetLayout()`
   (`js-api/src/views/view.ts:610` — leaves only the grid). Pattern:
   `packages/HitTriage/src/app/hit-design-app.ts:687-689`.

3. **Round-trip through JSON: `fromJson` vs `fromViewState`.**
   `DG.ViewLayout.fromJson(json)` parses the FULL document (entity
   metadata + `viewStateMap` + `columns`) — what a `.layout` file or
   `layout.toJson()` produces. `DG.ViewLayout.fromViewState(state)`
   rebuilds from the BARE inner view-state JSON only — no entity
   identity, no column-matching contract (`DG-FACT-219`). Use the
   bare form when persisting layouts inside your own structures
   (HitTriage stores `template.layoutViewState: string` per campaign):
   ```typescript
   const json = layout.toJson();                       // full doc
   const restored = DG.ViewLayout.fromJson(json);

   const state = layout.viewState;                     // viewStateMap JSON only
   view.loadLayout(DG.ViewLayout.fromViewState(state)); // no column metadata
   ```
   `JSON.parse(layout.toJson())` exposes top-level keys `#type`,
   `viewStateMap`, `columns`, `name`, `friendlyName`, `author`,
   `createdOn` (`DG-FACT-221`). `ViewLayout` also carries
   `getUserDataValue(k)` / `setUserDataValue(k, v)` for arbitrary
   string metadata — article omits these (`DG-FACT-219`).

4. **Ship a `.layout` file with the package.**
   Drop the `.layout` next to its `.csv` under `<Package>/files/`,
   read with `_package.files.readAsText(...)`, then `fromJson` →
   `loadLayout`. Production sites precede the apply with
   `await DG.delay(100)` because Dart-side viewer construction is
   async even though `loadLayout` is synchronous on the JS side
   (`DG-FACT-223`):
   ```typescript
   const tv = grok.shell.addTableView(await grok.data.files.openTable(
     `${_package.webRoot}files/demo_files/demo_smiles.csv`));
   const s = await _package.files.readAsText('demo_files/Overview_demo.layout');
   await DG.delay(100);
   tv.loadLayout(DG.ViewLayout.fromJson(s));
   ```
   Pattern: `packages/Chem/src/demo/demo.ts:64-68, 181-184, 194-196`.

5. **Persist to the server gallery via `grok.dapi.layouts`.**
   `grok.dapi.layouts: LayoutsDataSource` extends `HttpDataSource<ViewLayout>`,
   inheriting `save`, `find`, `delete`, `filter`, `list`, `count`
   (`DG-FACT-220`). The article only mentions `list()` and
   `getApplicable()` — `save(layout)` is the JS equivalent of
   `View | Layout | Save to Gallery`. Always call `l.newId()` before
   saving a layout reconstructed from JSON, otherwise the server
   rejects the duplicate entity id:
   ```typescript
   const l = view.saveLayout();
   await grok.dapi.layouts.save(l);                              // direct from view

   const fromFile = DG.ViewLayout.fromJson(layoutString);
   fromFile.newId();                                             // mandatory
   await grok.dapi.layouts.save(fromFile);

   const found = await grok.dapi.layouts.filter(`friendlyName = "${l.friendlyName}"`).list();
   ```
   Pattern: `packages/ApiTests/src/dapi/layouts.ts:246-249`.

6. **Auto-apply on a new dataset via `getApplicable`.**
   `grok.dapi.layouts.getApplicable(df)` returns layouts whose column
   contract binds to `df`. The article numbers three matching rules
   but they are SUFFICIENT conditions, not sequential AND
   (`DG-FACT-222`): a column binds when (a) name + type both match,
   OR (b) both columns carry the same `layout-id` tag, OR (c) both
   carry the same `quality` semantic-type tag. To force a stable
   match across datasets where column names differ, set
   `c.tags[DG.TAGS.LAYOUT_ID] = '<shared-id>'` on both source and
   target columns:
   ```typescript
   const df = grok.data.demo.demog();
   const view = grok.shell.addTableView(df);
   df.col('age')!.tags[DG.TAGS.LAYOUT_ID] = 'demog-age';     // stable identity
   const layouts = await grok.dapi.layouts.getApplicable(df);
   if (layouts.length) view.loadLayout(layouts[0]);
   ```
   Pattern: `packages/ApiSamples/scripts/dapi/applicable-layouts.js`,
   `packages/UITests/src/views/layouts.ts:86,130,159`.

## Common failure modes

- **`view.saveLayout()` returns null / throws.** Called on a
  non-`TableView` (`DG-FACT-218`). Guard with
  `if (view instanceof DG.TableView)` before the call.
- **`grok.dapi.layouts.save(layout)` errors with a duplicate-entity
  id.** Layout was built via `DG.ViewLayout.fromJson(...)`, which
  preserves the source entity id. Call `layout.newId()` before `save`
  (`packages/ApiTests/src/dapi/layouts.ts:248`).
- **Freshly saved layout missing from `getApplicable()` in the same
  session.** Server index lag — keep the in-memory `ViewLayout`
  reference and pass it to `loadLayout` directly instead of round-
  tripping through `getApplicable` (`grid-run.md:31,64` in
  `packages/UsageAnalysis/files/TestTrack/Viewers/`).
- **`getApplicable(df)` returns `[]` for a layout that visually
  should match.** Column names + types differ AND no `layout-id` /
  `quality` tag bridges them (`DG-FACT-222`). Set
  `col.tags[DG.TAGS.LAYOUT_ID] = '<shared-id>'` on both columns and
  re-save the source layout.
- **`loadLayout(layout)` silently overwrote the destination frame's
  column tags.** You hit `pickupColumnTags = true`. Re-load with the
  second arg explicitly omitted (`DG-FACT-DRIFT-086`).
- **Layout applies but viewers paint blank for a moment.** Dart-side
  viewer construction is async; insert `await DG.delay(100)` before
  `loadLayout` (`packages/Chem/src/demo/demo.ts:67`, `DG-FACT-223`).

## Verification

- TypeScript build (`npm run build` or `grok check`) exits `0` with
  no `Property 'saveLayout' does not exist on type 'View'` errors —
  every call site is on a `TableView`.
- `JSON.parse(view.saveLayout().toJson())` exposes `#type`,
  `viewStateMap`, `columns`, `name`, `author` top-level keys
  (`DG-FACT-221`).
- After a save-modify-restore round trip, `[...view.viewers]` length
  matches the count captured before `saveLayout`.
- `(await grok.dapi.layouts.list()).some(l => l.id === saved.id)` is
  `true` after `grok.dapi.layouts.save(saved)`.

## See also

- Source articles:
  - `help/develop/how-to/views/layouts.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-218` … `DG-FACT-223` and drift `DG-FACT-DRIFT-086`.
- Reference packages:
  - `packages/Chem/src/demo/demo.ts:64-68, 181-184, 194-196` —
    file-shipped `.layout`: `readAsText` → `fromJson` → `loadLayout`.
  - `packages/HitTriage/src/app/hit-design-app.ts:647-674, 1075-1079` —
    dual storage (file path first, in-template `viewState` fallback);
    save-back via `view.saveLayout().viewState`.
  - `packages/ApiTests/src/dapi/layouts.ts:219-249` — `getApplicable`,
    `filter`, `save` after `newId()`, `delete`.
  - `packages/UITests/src/views/layouts.ts:86-93, 130-133, 159-162` —
    `getApplicable(df).sort(... createdOn ...)` cleanup pattern.
- Related skills: `manipulate-viewers` (defines the viewer attach /
  docking primitives whose configuration this layout captures);
  `user-settings-storage` (per-user storage instead of the gallery).
