---
name: ui
description: UI building guidelines for Datagrok TypeScript plugins and libraries. Use when creating or reviewing UI components, viewers, dialogs, file viewers, or layouts.
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

## File Viewers

- When the result is a single DataFrame, return a `DG.TableView` directly — do not
  wrap a grid in a custom `DG.View`
- Transfer all source metadata to tags (table-level and column-level) so it survives
  round-trips through the platform

```typescript
// Good — file viewer for a single-table format
static async previewFoo(file: DG.FileInfo): Promise<DG.View> {
  const bytes = await file.readAsBytes();
  const data = await parseFoo(bytes);
  const df = toDataFrame(data);
  const view = DG.TableView.create(df, false);
  view.name = file.name;
  return view;
}
```

## Metadata Preservation

- Set `df.setTag('source.format', '...')` to record where the data came from
- Preserve column-level metadata (description, format, categories) as tags
- Use a namespace prefix for format-specific tags (e.g., `minitab.version`, `prism.sheetId`)
- Use the standard `format` tag for display formatting (e.g., `#.00` for 2 decimal places)
- Use the standard `description` tag for column descriptions

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
- Use `ui.input.choice()`, `ui.input.int()`, etc. for typed inputs
- For property panels, prefer `DG.JsViewer` properties (`this.string(...)`, `this.int(...)`)
  which automatically appear in the context panel

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

## Common Anti-Patterns

- Do not use raw `document.createElement` when `ui.*` helpers exist
- Do not set `innerHTML` with user data — use `ui.divText()` or `textContent`
- Do not use `style.width = '100%'` on tables — let them size to content
- Do not build all content eagerly in multi-pane layouts (tabs, accordions)
