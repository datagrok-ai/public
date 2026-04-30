---
name: ui
description: UI building guidelines for Datagrok TypeScript components, viewers, and drag-and-drop
when-to-use: When creating or modifying UI components, viewers, dialogs, file viewers, layouts, grids, or drag-and-drop
effort: low
---

# Datagrok UI Building Guidelines

Rules and patterns for building UI in TypeScript packages and libraries. These are
authoritative — follow them unless the user explicitly overrides.

## Tab Controls

- Always use `addPane` with a lazy `getContent` callback — never pass pre-built elements
- This avoids creating heavy components (grids, viewers, DataFrames) for tabs the user may never open

```typescript
// Good
const tabs = ui.tabControl();
tabs.addPane('Sheet 1', () => {
  const df = buildDataFrame();
  const grid = DG.Viewer.grid(df);
  return grid.root;
});

// Bad — all tabs built eagerly
const tabs = ui.tabControl({
  'Sheet 1': buildExpensiveGrid(),
  'Sheet 2': buildExpensiveGrid(),
});
```

## Grids and Viewers

- Prefer `DG.Viewer.grid(df)` for embedding a grid inside a composite layout
- Prefer `DG.TableView.create(df, false)` when the grid IS the entire view
- The second argument (`false`) prevents the table from being added to the workspace

## Layouts

- Use `ui.splitH` / `ui.splitV` for resizable split panels
- Use `ui.divV` / `ui.divH` for simple stacking without resize handles
- Set `flex: 1` on the element that should fill remaining space
- For tree + content layouts, use `ui.splitH([tree.root, contentPanel])`

## Dialogs and Inputs

- Use `ui.dialog()` for modal interactions
- Use `ui.input.choice()`, `ui.input.int()`, `ui.input.bool()`, etc. for typed inputs - full set of input functions is in `js-api/ui.ts`.
Prefer `onValueChanged` in the options object over `.onChanged.subscribe()`:
```typescript
  ui.input.bool('Debug', {value: DG.Test.isInDebug, onValueChanged: (v) => DG.Test.isInDebug = v});
```
Use `ui.form([...inputs])` to render a labeled list of inputs inside a dialog
- For property panels, prefer `DG.JsViewer` properties (`this.string(...)`, `this.int(...)`)
  which automatically appear in the context panel

## Toggle Settings in Popup Menus

For toggleable settings in a `DG.Menu.popup()`, use `menu.items()` with `isChecked` — **never** use text-prefix hacks like `` `${flag ? '✓ ' : ''}Label` ``:

```typescript
const toggles = [
  {label: 'Debug', get: () => DG.Test.isInDebug, set: (v: boolean) => { DG.Test.isInDebug = v; }},
  {label: 'Benchmark', get: () => DG.Test.isInBenchmark, set: (v: boolean) => { DG.Test.isInBenchmark = v; }},
];
const menu = DG.Menu.popup();
menu.closeOnClick = false;
const refresh = () => {
  menu.clear();
  menu.items(toggles, (t) => { t.set(!t.get()); refresh(); }, {isChecked: (t) => t.get()});
};
refresh();
menu.show();
```

## Performance

- Debounce resize and selection handlers: `DG.debounce(observable, 50)`
- For large datasets, prefer canvas-based rendering over DOM elements
- Avoid re-creating viewers on every data change — update in place when possible

## Accordion

- Use `ui.accordion()` with lazy `getContent` callbacks (same principle as tab controls)

```typescript
const acc = ui.accordion();
acc.addPane('Details', () => buildDetailsPanel());
acc.addPane('Statistics', () => buildStatsPanel());
```

## Drag and Drop

- `ui.makeDroppable(el, IDragAndDropOptions<T>)` — receive entities dragged from the browse tree, grid, or other sources.
  - `acceptDrop(obj)` — fast predicate for showing the zone.
  - `doDrop(args)` — handle the drop. `args` is a `DragDropArgs<T>` with `dragObject`, `dragSource`, `dragObjectType`, `copying` (Ctrl/Cmd), `link` (Alt), `handled`.
  - Rich hooks: `acceptDrag`, `onBeginDrag`, `onEndDrag`, `onMouseEnter/Over/Leave/Out`, `dropSuggestion`, `makeDropZone`, `dropZoneRectTransformation`, `dropIndication`.
- `ui.makeDraggable(el, {getDragObject, getDragCaption})` — make your own UI a drag source.
- Sample: `ApiSamples/scripts/ui/interactivity/drag-and-drop.js`.

## Common Anti-Patterns

- Do not use raw `document.createElement` when `ui.*` helpers exist
- Do not set `innerHTML` with user data — use `ui.divText()` or `textContent`
- Do not use `style.width = '100%'` on tables — let them size to content
- Do not build all content eagerly in multi-pane layouts (tabs, accordions)
