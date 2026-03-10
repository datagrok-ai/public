# CLAUDE.md - Flow Package

**IMPORTANT: Always update this file after adding new features, nodes, or changing architecture.**

## Overview

Flow (FuncFlow) is an interactive visual function chain designer for Datagrok. It uses LiteGraph.js to let users compose Datagrok functions into executable JavaScript scripts via a node-based graph editor.

## Architecture

```
src/
├── package.ts                    # Entry: @funcflowApp, @fileViewer decorators
├── funcflow-view.ts              # Main view: 2-panel layout + native context panel, ribbon, status bar
├── canvas/
│   ├── canvas-controller.ts      # LGraphCanvas wrapper, type validation, theme
│   └── graph-manager.ts          # LGraph state wrapper, change notifications
├── nodes/
│   ├── node-factory.ts           # Node registration, DG.Func discovery
│   ├── input-nodes.ts            # Source nodes (become //input: lines)
│   ├── output-nodes.ts           # Sink nodes (become //output: lines)
│   ├── func-node.ts              # Dynamic node generator for DG.Func
│   ├── utility-nodes.ts          # Utility + constant nodes
│   └── comparison-nodes.ts       # Comparison operator nodes
├── compiler/
│   ├── graph-compiler.ts         # Graph → CompiledStep[]
│   ├── script-emitter.ts         # CompiledStep[] → JavaScript source
│   ├── validator.ts              # Pre-compilation validation (cycles, required inputs, duplicates)
│   ├── topological-sort.ts       # Kahn's algorithm
│   └── graph-utils.ts            # Safe LGraph node accessor
├── panel/
│   ├── function-browser.ts       # Left sidebar: searchable function catalog
│   └── property-panel.ts         # Context panel content: node properties editor with tooltips
├── serialization/
│   ├── flow-schema.ts            # .ffjson type definitions
│   └── flow-serializer.ts        # Load/save/download
├── history/
│   └── undo-manager.ts           # Command pattern (infrastructure, not fully wired)
├── types/
│   ├── funcflow-node.ts          # FuncFlowNode interface extends LGraphNode
│   ├── type-map.ts               # DG type→slot type+color, type compatibility
│   └── litegraph-augment.d.ts    # TypeScript augmentations for LiteGraph
└── utils/
    └── dart-proxy-utils.ts       # Safe Dart proxy property access
```

## Data Flow

```
User double-clicks function in FunctionBrowser
  → CanvasController.addNodeAtCenter(nodeTypeName)
  → User connects nodes visually
  → Generate Script action:
      validateGraph() → compileGraph() → emitScript()
      → JavaScript source with //input: //output: annotations
```

## Node Types

### Input Nodes (src/nodes/input-nodes.ts)
Become `//input:` annotation lines in generated scripts.

| Node | DG Type | Qualifiers |
|------|---------|------------|
| Table Input | dataframe | - |
| Column Input | column | type, semType filters |
| Column List Input | column_list | type, semType filters |
| String Input | string | nullable, choices, caption, semType |
| Number Input | double | nullable, min, max, showSlider, caption |
| Int Input | int | nullable, min, max, showSlider, caption |
| Boolean Input | bool | nullable, caption |
| DateTime Input | datetime | nullable, caption |
| File Input | file | nullable, caption |
| Map Input | map | nullable, caption |
| Dynamic Input | dynamic | nullable, caption |
| String List Input | string_list | - |
| Blob Input | blob | nullable, caption |

### Output Nodes (src/nodes/output-nodes.ts)
Become `//output:` annotation lines.

| Node | Notes |
|------|-------|
| Table Output | Fixed dataframe type |
| Value Output | Configurable type: string, int, double, bool, dataframe, column, column_list, object, dynamic, list, view, viewer, widget, graphics, grid_cell_renderer, filter, map, datetime, blob, funccall |

### Utility Nodes (src/nodes/utility-nodes.ts)
Generate inline code in script body.

| Node | Code Generated |
|------|----------------|
| Select Column | `df.col('name')` |
| Select Columns | `[df.col('a'), df.col('b')]` |
| Add Table View | `grok.shell.addTableView(df)` |
| Log | `console.log(value)` |
| Info | `grok.shell.info(msg)` |
| Warning | `grok.shell.warning(msg)` |
| ToString | `(value).toString()` |
| FromJSON | `JSON.parse(json)` |
| ToJSON | `JSON.stringify(value)` |
| Constants | String, Int, Double, Boolean, List literals |

### Comparison Nodes (src/nodes/comparison-nodes.ts)
==, !=, >, >=, <, <=, Contains, StartsWith, EndsWith, IsNull

### Function Nodes (src/nodes/func-node.ts)
Dynamically created for each DG.Func. Generate `await grok.functions.call(...)`.

## Function Filtering (src/nodes/node-factory.ts)

Functions are filtered during registration using `EXCLUDED_TAGS` and `EXCLUDED_ROLES` constants. Functions with any excluded tag or role are skipped. Edit these arrays in `node-factory.ts` to control which DG functions appear.

## Script Generation

The compiler pipeline:
1. **Topological sort** - Kahn's algorithm orders nodes
2. **Compile** - Each node becomes a CompiledStep with inputs/outputs resolved
3. **Emit** - Steps become JavaScript lines with `//input:`, `//output:` headers

### Validation (`validator.ts`)
Runs before compilation. Checks:
- Empty graph (warning)
- Cycles via topological sort (error)
- Column input without Table input (error)
- **Required (non-nullable) func inputs**: unconnected inputs with no stored value are errors. `0` and `false` are valid values; only `undefined`/`null`/empty string = missing. Nullable check: `param.nullable || param.options.optional || param.options.nullable`
- Disconnected nodes (warning)
- Output nodes with no incoming connection (warning)
- Empty or invalid parameter names (error)
- Duplicate parameter names (error)

Input annotation qualifiers are generated from node widget values:
- `{type: numerical}` from Column type filter
- `{semType: Molecule}` from Column semType filter or String semType combo
- `{nullable: true}` from nullable toggle
- `{caption: Name}` from caption field
- `{choices: ["a","b"]}` from choices field
- `{min: 0; max: 100}` from min/max fields
- `{showSlider: true}` from showSlider toggle (Number/Int only)

## File Format

`.ffjson` files store the full flow state including LiteGraph graph data and FuncFlow metadata.

## Build

```bash
cd packages/Flow
npm install
npm run build    # grok api && grok check --soft && webpack
```

## UI Architecture: Widgets vs Property Panel

**Nodes have NO inline widgets** (except ConstStringNode). All property editing happens in Datagrok's native context panel (right side) when a node is selected.

### Context Panel Integration
- `grok.shell.windows.showContextPanel = true` enables the native panel
- `grok.shell.o = propertyPanel.root` sets its content on node selection
- Layout is 2-panel: `[leftPanel (FunctionBrowser), canvasContainer]` — no custom right panel

### Property Panel (`property-panel.ts`)
- Reads `node.properties` and creates appropriate editors:
  - `string` properties → `<textarea>` (auto-resizing)
  - `number` properties → `<input type="number">` with step (1 for int, 0.1 for double)
  - `boolean` properties → `<input type="checkbox">`
  - enum properties → `<select>` dropdown
- **Exception**: `ConstStringNode` keeps its inline `text` widget for quick editing
- **Func nodes**: primitive input defaults stored as `_input_${name}` properties, edited in panel
- **Collapse icon**: `NODE_DEFAULT_BOXCOLOR = '#888'` (visible on light theme)
- **Font**: Roboto via Google Fonts import + LiteGraph canvas font overrides

### Tooltip System
- **PROP_TOOLTIPS**: Maps label names (e.g. 'Param Name', 'Nullable') to tooltip strings for input/output node properties
- **UTILITY_PROP_TOOLTIPS**: Maps `node.title → property name → tooltip` for utility/constant nodes (e.g. List → value → "Comma-separated list of values")
- **buildFuncInputTooltip(param)**: Builds rich tooltips for DG.Func input parameters from `DG.Property` metadata: description, type, default value, nullable status
- Tooltips are bound to both **labels** and **input elements** via `ui.tooltip.bind()`
- All editor helpers (`createTextarea`, `createNumberInput`, `createToggle`, `createCombo`) accept optional `inputTooltip` parameter

### Function Browser Tooltips
- All built-in nodes (inputs, outputs, utilities, constants, comparisons) have descriptive tooltips
- Format: `"<description>. Double-click to add"` — e.g. "Dataframe input parameter. Double-click to add"
- DG function nodes show their `func.description` plus package name

## Key Dependencies

- `litegraph.js` ^0.7.18 - Graph canvas library
- `datagrok-api` - Platform API (external, not bundled)

## Development Guidelines

1. **Adding new input nodes**: Add class in `input-nodes.ts` (properties only, no widgets), register in `registerInputNodes()`, add to `inputNodes` array with `desc` in `function-browser.ts`, handle in `buildInputLine()` in `script-emitter.ts` if special qualifiers needed. Property panel auto-discovers properties by key name.
2. **Adding new utility nodes**: Add class in `utility-nodes.ts` (properties only, no widgets), register in `registerUtilityNodes()`, add to `utilityNodes` array with `desc` in `function-browser.ts`, add case in `emitUtilityStep()` in `script-emitter.ts`. Add property tooltips to `UTILITY_PROP_TOOLTIPS` in `property-panel.ts`.
3. **Adding new types**: Add to `DG_TYPE_MAP` in `type-map.ts`, add compatibility rules to `COMPATIBLE_TYPES`
4. **After any change**: Update this CLAUDE.md file
