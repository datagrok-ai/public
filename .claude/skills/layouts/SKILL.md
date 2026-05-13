---
name: layouts
version: 0.1.0
description: |
  Capture a `TableView`'s viewer arrangement (positions, properties,
  column tags) into a portable artifact your package can ship, persist
  to the server gallery, or auto-apply when a fresh dataset opens.
  Covers the JS API path — `view.saveLayout` / `loadLayout`,
  `DG.ViewLayout.fromJson` / `fromViewState`, `grok.dapi.layouts`, and
  the `layout-id` / semantic-type column tags the matching algorithm
  uses — that the help article describes only partially.
  Use when asked to "save and restore a multi-viewer arrangement",
  "ship a dashboard preset alongside a CSV", or "auto-apply a saved
  viewer setup when a matching table opens".
triggers:
  - save and restore a viewer arrangement
  - ship a dashboard preset with a dataset
  - reapply a captured viewer setup
  - persist a multi-viewer dashboard to the gallery
  - auto-apply a saved viewer arrangement on table open
  - round-trip viewer state through json
allowed-tools:
  - Read
  - Edit
  - Bash
harness-authored: true
---

# layouts

## When to use

You are scripting viewer-arrangement reuse from a package or app —
capturing the current view to ship next to a CSV, restoring a saved
arrangement on a fresh frame, persisting per-campaign presets on the
server gallery, or letting the platform auto-pick an applicable preset
for newly opened tables.

## Prerequisites

- A package scaffold (`grok create <Name>`); commands run from the
  package root, code lives under `src/`.
- `datagrok-api` imports — the article omits these:
  ```typescript
  import * as grok from 'datagrok-api/grok';
  import * as DG   from 'datagrok-api/dg';
  ```
- A `DG.TableView` (`grok.shell.addTableView(df)` or `grok.shell.tv`).
  `saveLayout` / `loadLayout` are defined on the generic `View` base
  but are no-ops on a plain `ViewBase` — only `TableView`-shaped views
  serialize meaningful state (`js-api/src/views/view.ts:297-307`).
- For server-gallery persistence: an authenticated `grok.dapi` session
  (default inside any package running in the platform).

## Steps

1. **Capture the current view's arrangement.**
   `view.saveLayout()` returns a `DG.ViewLayout`. The article omits the
   optional flag `saveWithData` (default `false`); pass `true` only
   when shipping a self-contained `.layout` that carries the underlying
   frame inline (`js-api/src/views/view.ts:303-307`):
   ```typescript
   const view = grok.shell.addTableView(grok.data.demo.demog());
   const portable = view.saveLayout();                    // data lives elsewhere
   const bundled  = view.saveLayout({saveWithData: true}); // inlines table data
   ```
   Expected: `portable` is a non-null `DG.ViewLayout`. Guard with
   `if (view instanceof DG.TableView)` if the caller's `view` may be a
   non-table `ViewBase`.

2. **Reapply / roll back via `loadLayout`.**
   `loadLayout(layout, pickupColumnTags?)` accepts a second flag the
   article doesn't mention (default `false`); set `true` only when the
   stored column tags should be copied onto the destination frame
   (`js-api/src/views/view.ts:296-299`):
   ```typescript
   const layout = view.saveLayout();
   view.addViewer('Histogram', {value: 'age'});
   view.loadLayout(layout);          // rolls back the added histogram
   view.loadLayout(layout, true);    // also overwrites column tags
   ```
   To clear without a saved layout, use `view.resetLayout()`
   (`js-api/src/views/view.ts:610` — leaves only the grid). Pattern:
   `packages/HitTriage/src/app/hit-design-app.ts:687-689`.

3. **Round-trip through JSON: `fromJson` vs `fromViewState`.**
   `DG.ViewLayout.fromJson(json)` parses the FULL document (entity
   metadata + `viewStateMap` + `columns`) — what a `.layout` file or
   `layout.toJson()` produces. `DG.ViewLayout.fromViewState(state)`
   rebuilds from the BARE inner `viewState` JSON only, no entity
   identity, no column-matching contract
   (`js-api/src/entities/view-layout.ts:22-36`). Prefer the bare form
   when persisting layouts inside your own structures (HitTriage stores
   `template.layoutViewState: string` per campaign):
   ```typescript
   const restored = DG.ViewLayout.fromJson(layout.toJson()); // full doc
   view.loadLayout(DG.ViewLayout.fromViewState(layout.viewState)); // state only
   ```
   `JSON.parse(layout.toJson())` exposes top-level keys `#type`,
   `viewStateMap`, `columns`, `name`, `friendlyName`, `author`,
   `createdOn`. `ViewLayout` also carries `getUserDataValue(k)` /
   `setUserDataValue(k, v)` (`js-api/src/entities/view-layout.ts:38-44`).

4. **Ship a `.layout` file with the package.**
   Drop the `.layout` next to its `.csv` under `<Package>/files/`,
   read with `_package.files.readAsText(...)`, then `fromJson` →
   `loadLayout`. Production sites precede the apply with
   `await DG.delay(100)` because Dart-side viewer construction is async
   even though `loadLayout` is synchronous on the JS side:
   ```typescript
   const tv = grok.shell.addTableView(await grok.data.files.openTable(
     `${_package.webRoot}files/demo_files/demo_smiles.csv`));
   const s = await _package.files.readAsText('demo_files/Overview_demo.layout');
   await DG.delay(100);
   tv.loadLayout(DG.ViewLayout.fromJson(s));
   ```
   Pattern: `packages/Chem/src/demo/demo.ts:64-68, 181-184, 194-196`.

5. **Persist to the server gallery via `grok.dapi.layouts`.**
   `grok.dapi.layouts` extends `HttpDataSource<ViewLayout>`, inheriting
   `save`, `find`, `delete`, `filter`, `list`, `count`
   (`js-api/src/dapi.ts:111-113, 645-660`). The article only names
   `list()` and `getApplicable()` — `save(layout)` is the JS equivalent
   of `View | Layout | Save to Gallery`. Always call `l.newId()` before
   saving a layout reconstructed from JSON, otherwise the server
   rejects the duplicate entity id:
   ```typescript
   await grok.dapi.layouts.save(view.saveLayout());       // direct from view

   const fromFile = DG.ViewLayout.fromJson(layoutString);
   fromFile.newId();                                      // mandatory
   await grok.dapi.layouts.save(fromFile);
   ```
   Pattern: `packages/ApiTests/src/dapi/layouts.ts:246-249`.

6. **Auto-apply on a new dataset via `getApplicable`.**
   `grok.dapi.layouts.getApplicable(df)` returns layouts whose column
   contract binds to `df`. The article numbers three matching rules
   but they are SUFFICIENT conditions, not sequential AND: a column
   binds when (a) name + type both match, OR (b) both columns carry
   the same `layout-id` tag (`js-api/src/const.ts:283`), OR (c) both
   carry the same `quality` semantic-type tag. To force a stable match
   across datasets where column names differ, set
   `c.tags[DG.TAGS.LAYOUT_ID] = '<shared-id>'` on both source and
   target columns before saving:
   ```typescript
   const df = grok.data.demo.demog();
   const view = grok.shell.addTableView(df);
   df.col('age')!.tags[DG.TAGS.LAYOUT_ID] = 'demog-age';
   const layouts = await grok.dapi.layouts.getApplicable(df);
   if (layouts.length) view.loadLayout(layouts[0]);
   ```
   Pattern: `packages/UITests/src/views/layouts.ts:86,130,159`.

## Common failure modes

- **`view.saveLayout()` returns null / throws.** Called on a
  non-`TableView`. Guard with `if (view instanceof DG.TableView)`
  before the call (`js-api/src/views/view.ts:305-307`).
- **`grok.dapi.layouts.save(layout)` errors with a duplicate-entity
  id.** Layout was built via `DG.ViewLayout.fromJson(...)`, which
  preserves the source entity id. Call `layout.newId()` before `save`
  (`packages/ApiTests/src/dapi/layouts.ts:248`).
- **`getApplicable(df)` returns `[]` for a frame that visually should
  match.** Column names + types differ AND no `layout-id` / semantic
  type bridges them. Set `col.tags[DG.TAGS.LAYOUT_ID] = '<shared-id>'`
  on both columns and re-save the source.
- **`loadLayout(layout)` silently overwrote the destination frame's
  column tags.** You hit `pickupColumnTags = true`. Re-load with the
  second arg explicitly omitted.
- **Layout applies but viewers paint blank for a moment.** Dart-side
  viewer construction is async; insert `await DG.delay(100)` before
  `loadLayout` (`packages/Chem/src/demo/demo.ts:67`).

## Verification

- TypeScript build (`npm run build` or `grok check`) exits `0` with no
  `Property 'saveLayout' does not exist on type 'View'` errors —
  every call site is on a `TableView`.
- `JSON.parse(view.saveLayout().toJson())` exposes `#type`,
  `viewStateMap`, `columns`, `name`, `author` top-level keys.
- `(await grok.dapi.layouts.list()).some(l => l.id === saved.id)` is
  `true` after `grok.dapi.layouts.save(saved)`.

## See also

- Source articles:
  - `help/develop/how-to/views/layouts.md`
- Reference packages:
  - `packages/Chem/src/demo/demo.ts:64-68, 181-184` —
    file-shipped `.layout`: `readAsText` → `fromJson` → `loadLayout`.
  - `packages/HitTriage/src/app/hit-design-app.ts:647-674` — dual
    storage (file path first, in-template `viewState` fallback).
  - `packages/ApiTests/src/dapi/layouts.ts:219-249` — `getApplicable`,
    `filter`, `save` after `newId()`, `delete`.
  - `packages/UITests/src/views/layouts.ts:86-93, 130-133` —
    `getApplicable(df)` cleanup pattern.
- Related skills: `manipulate-viewers` (defines the viewer attach /
  docking primitives whose configuration this captures);
  `user-settings-storage` (per-user storage instead of the gallery).
