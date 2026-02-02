# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is `datagrok-api`, the JavaScript/TypeScript API for the Datagrok platform - a data analytics and visualization 
platform. The API provides TypeScript bindings that communicate with a Dart backend via an interop layer.

## Build Commands

```bash
npm run build           # Compile TypeScript and bundle with webpack
npm run build-ts        # TypeScript compile only + ESLint fix
npm run build-js-api    # Clean src, compile TypeScript, bundle webpack
npm run build-docker    # Build for Docker deployment
npm run build-all       # Build chem-meta library dependency first, then build
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

The API mostly wraps a Dart backend (in `../../core/client`). Key patterns:

- **`api` object**: Global interface to Dart functions (`IDartApi` from `src/api/grok_api.g.ts`)
- **`toJs(dart)`**: Convert Dart objects to JavaScript (`wrappers.ts`)
- **`toDart(js)`**: Convert JavaScript objects to Dart (`wrappers.ts`)
- **`.dart` property**: Most API classes hold a reference to their Dart counterpart

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

## Code Style

- 2-space indentation
- Single quotes for strings
- Semicolons required
- Windows line endings (CRLF)
- TypeScript strict mode with `strictNullChecks: true`
- File naming convention is kebab-case: `my-component.ts`

## Tests

API should be tested by running the [ApiTests](../packages/ApiTests) tests that 
are executed within puppeteer with the running Datagrok. Whenever a change is made
to JS API, tests should be passed.

```
grok test
```

## Canonical code samples

See [API usage samples](../packages/ApiSamples/scripts) 

Each sample can be executed with `eval` in plain JavaScript within the running Datagrok in the browser.

When new important functionality is created, a sample needs to be added as well.