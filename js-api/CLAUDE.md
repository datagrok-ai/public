# CLAUDE.md

## Overview

This is `datagrok-api`, the JavaScript/TypeScript API for the Datagrok platform - a data analytics and visualization 
platform. The API provides TypeScript bindings that communicate with a Dart backend via an interop layer.

## Build Commands

```bash
npm install
npm run build           # tsc && webpack
npm run build-ts        # TypeScript + ESLint fix only
npm run build-js-api    # Clean, compile, bundle
```

## Architecture

### Entry Points (Three Main Namespaces)

The API is organized into three main namespaces that packages import:

- **`grok`** (`grok.ts`) - High-level APIs: shell access, functions, events, data operations, server API, user settings, AI
- **`ui`** (`ui.ts`) - UI components: elements, dialogs, inputs, menus, accordions, forms, viewers
- **`dg`** (`dg.ts`) - Re-exports all types, constants, and classes from the entire API

### Core Modules (in `src/`)

| Module          | Purpose                                                                        |
|-----------------|--------------------------------------------------------------------------------|
| `dataframe.ts`  | Barrel export for `src/dataframe/` - DataFrame, Column, Row, BitSet            |
| `entities.ts`   | Barrel export for `src/entities/` - User, Group, Project, Package, etc.        |
| `widgets.ts`    | Barrel export for `src/widgets/` - Widget, Dialog, Menu, InputBase, etc.       |
| `viewer.ts`     | Viewer base class and built-in viewers (ScatterPlot, Histogram, etc.)          |
| `grid.ts`       | Grid viewer, cell rendering, GridCell, GridColumn                              |
| `dapi.ts`       | Server API (Dapi class) - HTTP data sources for entities, queries, files, etc. |
| `shell.ts`      | Shell class - access to views, tables, windows, settings, current objects      |
| `functions.ts`  | FuncCall, function registration, parameter handling                            |
| `events.ts`     | Event system - rxjs Observable-based platform events                           |
| `views/view.ts` | ViewBase, View, TableView - application views                                  |
| `const.ts`      | Enums and constants (TYPE, COLUMN_TYPE, AGG, JOIN_TYPE, VIEWER, etc.)          |

### Module Directories

Large modules are split into focused sub-modules for maintainability. The main files (`dataframe.ts`, `entities.ts`, `widgets.ts`) are barrel exports that re-export everything from their directories, preserving backward compatibility.

#### `src/dataframe/` - Core Data Structures
| File                | Contents                                                    |
|---------------------|-------------------------------------------------------------|
| `types.ts`          | Type aliases: RowPredicate, Comparer, ColumnId, etc.        |
| `qnum.ts`           | Qnum class for qualified numbers (with comparison operators)|
| `bit-set.ts`        | BitSet class for efficient boolean arrays                   |
| `stats.ts`          | Stats class, GroupByBuilder for aggregations                |
| `column.ts`         | Column and typed variants (FloatColumn, DateTimeColumn)     |
| `column-list.ts`    | ColumnList collection class                                 |
| `column-helpers.ts` | ColumnMetaHelper, ColumnColorHelper, ColumnMarkerHelper     |
| `row.ts`            | Row, Cell, RowList, RowGroup, RowMatcher, ValueMatcher      |
| `data-frame.ts`     | DataFrame class and helper classes                          |
| `formula-helpers.ts`| DataFrameFormulaLinesHelper, DataFrameAnnotationRegionsHelper|

#### `src/entities/` - Platform Entities
| File                | Contents                                                    |
|---------------------|-------------------------------------------------------------|
| `types.ts`          | PropertyGetter, PropertySetter, DataConnectionProperties    |
| `entity.ts`         | Base Entity class                                           |
| `user.ts`           | User, UserSession, Group (kept together for circular deps)  |
| `property.ts`       | IProperty, Property, EntityProperty                         |
| `func.ts`           | Func, Script, ScriptEnvironment                             |
| `data-connection.ts`| DataConnection, DataQuery, TableQuery, DataJob, Credentials |
| `table-info.ts`     | TableInfo, ColumnInfo, FileInfo                             |
| `logging.ts`        | LogEventType, LogEvent, LogEventParameter                   |
| `schema.ts`         | Schema, EntityType, HistoryEntry                            |
| `project.ts`        | Project, ProjectOpenOptions                                 |
| `view-layout.ts`    | ViewLayout, ViewInfo                                        |
| `reports.ts`        | UserReport, UserReportsRule, UserNotification               |
| `misc.ts`           | Model, Notebook, Package, DockerContainer, ProgressIndicator|
| `search-provider.ts`| SearchProvider types                                        |

#### `src/widgets/` - UI Widgets
| File                | Contents                                                    |
|---------------------|-------------------------------------------------------------|
| `types.ts`          | RangeSliderStyle, SliderOptions, IMenu* interfaces          |
| `base.ts`           | Widget, DartWidget, DartWrapper, ObjectPropertyBag          |
| `filter.ts`         | Filter abstract base class                                  |
| `func-call-editor.ts`| FuncCallEditor abstract class                              |
| `inputs-base.ts`    | InputBase, JsInputBase                                      |
| `code-editor.ts`    | CodeEditor                                                  |
| `inputs.ts`         | DateInput, ChoiceInput, TypeAhead, CodeInput                |
| `markdown-input.ts` | MarkdownInput (separate for Quill lazy loading)             |
| `forms.ts`          | Dialog, InputForm                                           |
| `containers.ts`     | Accordion, AccordionPane, TabControl, TabPane, ToolboxPage  |
| `menu.ts`           | Menu, Balloon                                               |
| `tree.ts`           | TagEditor, TagElement, TreeViewNode, TreeViewGroup          |
| `data-widgets.ts`   | RangeSlider, HtmlTable, ColumnComboBox, Legend, PropertyGrid|
| `specialized.ts`    | FilesWidget, FunctionsWidget, Favorites, VisualDbQueryEditor|
| `progress.ts`       | TaskBarProgressIndicator                                    |

### Dart-JavaScript Interop

The API wraps a Dart backend. Key patterns:

- **`api` object**: `window` cast as `IDartApi` — all `grok_*` functions are Dart handlers
- **`.dart` property**: every wrapper class holds the Dart handle — always pass `x.dart` to `api.grok_*`, never the wrapper itself
- **`toJs(dart)`**: wrap Dart results on the JS side for `reg`-based (sync) calls — `ra*` async calls do it automatically
- **`toDart(x)`**: extracts `.dart`, calls `.toDart()`, converts `dayjs`→DateTime, plain `{}`→Map

```typescript
// sync — needs toJs()
get columns(): ColumnList { return toJs(api.grok_DataFrame_Columns(this.dart)); }

// async — toJs() already applied
async find(id: string): Promise<Entity> { return toJs(await api.grok_MyClass_Find(this.dart, id)); }
```

### Generated API Files (in `src/api/`)

Files ending in `.api.g.ts` are auto-generated from Dart code:
- `grok_api.g.ts` - IDartApi interface (main Dart function bindings)
- `grok_shared.api.g.ts` - Shared types and enums
- `d4.api.g.ts`, `ddt.api.g.ts` - Additional generated types

Do not manually edit `.g.ts` files - they are regenerated from the Dart codebase.

### Webpack Configuration

Two build targets in `webpack.config.js`:
1. **Node target** - Outputs `datagrok.js` (CommonJS) for Node.js/server-side use
2. **Browser target** - Outputs `js-api.js` to `../../core/client/xamgle/web/js/api` for browser use

## Quick Lookups

For JS API method questions, check these files first:

| Looking for...                  | Check first                        |
|---------------------------------|------------------------------------|
| Expression/formula evaluation   | `src/functions.ts` (eval, call, scriptSync, Context) |
| Function calls, FuncCall        | `src/functions.ts`                 |
| DataFrame operations            | `src/dataframe/data-frame.ts`      |
| Column operations               | `src/dataframe/column.ts`          |
| UI components, dialogs, inputs  | `src/widgets/`                     |
| Filter panel                    | `src/viewer.ts` (`FilterGroup`)    |
| Server/HTTP API                 | `src/dapi.ts`                      |
| Viewers                         | `src/viewer.ts`                    |
| Grid                            | `src/grid.ts` · [usage guide](src/grid.md) |
| Events                          | `src/events.ts`                    |
| Shell (views, tables, windows)  | `src/shell.ts`                     |
| Constants and enums             | `src/const.ts`                     |

## Filter Panel (`FilterGroup`)

Access via `tv.getFiltersGroup()` on a `TableView`. Use `{ createDefaultFilters: false }` in plugins to avoid adding unwanted default filters.

```typescript
const fg = tv.getFiltersGroup({ createDefaultFilters: false });
fg.updateOrAdd({ type: DG.FILTER_TYPE.HISTOGRAM, column: 'age', min: 20, max: 60 });
fg.updateOrAdd({ type: DG.FILTER_TYPE.CATEGORICAL, column: 'sex' }, false); // false = defer filtering
await DG.delay(200); // wait for the filter to be applied
```

- `fg.filters` — mixed array: built-in filters are **opaque Dart handles**, custom JS filters are `DG.Filter` instances — don't introspect built-in ones directly, use `fg.getStates(colName, filterType)` instead
- `fg.setEnabled(f, false)` — disable one filter; `fg.setActive(false)` — disable the whole group
- `column` is the canonical state field (`columnName` is a backwards-compat alias)
- Use `DG.delay` for timing

## Menu (`DG.Menu`)

Source: `src/widgets/menu.ts`. Option interfaces: `src/widgets/types.ts`.

```typescript
DG.Menu.popup()
  .item('Action', () => grok.shell.info('clicked'))
  .separator()
  .group('Sub').item('Inner', () => {}).endGroup()
  .items(['A', 'B'], (s) => use(s), { radioGroup: 'grp', isChecked: (s) => s === cur })
  .show();
```

- All methods return `Menu` (fluent chain) — read `menu.ts` for the full API
- Viewer context menu: `viewer.onContextMenu.subscribe((menu) => { menu.item(...); })`
- Samples: `ApiSamples/scripts/ui/components/popup-menu.js`, `menu-customization.js`, `menu-advanced.js`

## Server API Usage (grok.dapi)

Plugin code accesses the server via `grok.dapi` (instance of `Dapi` class from `src/dapi.ts`).

Key sub-objects: `users`, `groups`, `connections`, `queries`, `tables`, `projects`, `scripts`,
`packages`, `files`, `docker`, `permissions`, `layouts`, `views`, `functions`, `log`, `spaces`.

Most sub-objects extend `HttpDataSource<T>` providing: `list()`, `find(id)`, `save(entity)`,
`delete(entity)`, `filter(query)`, `order(field)`, `page(n)`, `by(pageSize)`.

For external HTTP requests from plugins, `grok.dapi.fetchProxy(url, params)` proxies through the
server to avoid CORS. Raw `fetch()` should never be used in plugin code.

See samples: `packages/ApiSamples/scripts/dapi/`

## Canonical code samples

See [API usage samples](../packages/ApiSamples/scripts) 

Each sample can be executed with `eval` in plain JavaScript within the running Datagrok in the browser.

When new important functionality is created, a sample needs to be added as well.