# CLAUDE.md - Flow Package

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
│   ├── comparison-nodes.ts       # Comparison operator nodes
│   └── breakpoint-node.ts        # Debug breakpoint node (pauses in debug mode)
├── compiler/
│   ├── graph-compiler.ts         # Graph → CompiledStep[]
│   ├── script-emitter.ts         # CompiledStep[] → JavaScript source (+ instrumented mode)
│   ├── validator.ts              # Pre-compilation validation (cycles, required inputs, duplicates)
│   ├── topological-sort.ts       # Kahn's algorithm
│   └── graph-utils.ts            # Safe LGraph node accessor
├── execution/
│   ├── execution-state.ts        # NodeExecStatus enum, NodeExecState, ExecutionState class
│   ├── execution-controller.ts   # Run lifecycle, event subscriptions, breakpoint control
│   ├── execution-visualizer.ts   # Maps execution state → node visual properties (colors, overlays)
│   └── value-inspector.ts        # Builds context panel section for runtime value display
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
User double-clicks function in FunctionBrowser (or drags DG.Func / FileInfo onto canvas)
  → CanvasController.addNodeAtCenter(nodeTypeName)
  → User connects nodes visually
  → Generate Script action:
      validateGraph() → compileGraph() → emitScript()
      → JavaScript source with //input: //output: annotations
```

## Drag-and-Drop

The canvas accepts drops from the Datagrok browse tree via `ui.makeDroppable()`:
- **DG.FileInfo** (files) → creates an OpenFile function node with the file path pre-filled
- **DG.Func** (queries, scripts, functions) → finds the matching registered node type and adds it

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

#### Node Display Name
Node headers use `func.friendlyName` (falling back to `func.name`), then split by `|` and take the last segment trimmed. E.g. `"Browse | CHEM | All ChEMBL structures"` → `"All ChEMBL structures"`. Helper: `getFuncDisplayName()` in `dart-proxy-utils.ts`.

#### Pass-Through Outputs
Every func node automatically gets **pass-through output slots** mirroring each input, named `inputName →` (with arrow suffix). These solve the execution ordering problem for mutating functions:

- **Problem**: If `addNewColumn(table)` and `doSomething(table)` both take the same table input, topological sort can't determine order since there's no edge between them.
- **Solution**: Connect `table` to `addNewColumn`'s input, then connect `addNewColumn`'s `table →` pass-through output to `doSomething`'s input. This creates a topological edge enforcing `addNewColumn` runs first.
- **Compiler behavior**: Pass-through outputs resolve to the same variable as the corresponding input — no new code is generated. `_passthroughCount` stores how many pass-through slots are at the start of the outputs array.
- **Visual layout**: Pass-through outputs come **first** (aligned with their corresponding input slots on the left), real outputs come **last** with arrow-shaped slots (`LiteGraph.ARROW_SHAPE`). The `→` suffix on pass-through names also distinguishes them.
- **Tooltips**: Hovering a pass-through slot shows "pass-through" label + hint to connect for execution ordering.

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

### Instrumented Mode (`EmitOptions`)

`emitScript()` accepts an optional `EmitOptions` parameter:

```typescript
interface EmitOptions {
  instrumented?: boolean;     // false = clean script, true = try/catch + events
  runId?: string;             // UUID for this run (required when instrumented)
  enableBreakpoints?: boolean; // emit breakpoint pause code (debug mode)
  haltOnError?: boolean;      // throw on first error (default true)
}
```

When `instrumented=true`, each step is wrapped in try/catch and fires custom events via `grok.events.fireCustomEvent()` on channel `funcflow.exec.{runId}`. Event types: `run-start`, `node-start`, `node-complete`, `node-error`, `breakpoint-hit`, `run-complete`. Output values are summarized (DataFrame → row/col count, Column → name + sample, etc.) to avoid large payloads.

**Variable hoisting**: When `wrapInstrumented()` encounters a code line starting with `let varName = ...`, it hoists the declaration (`let varName;`) before the `try` block and puts only the assignment (`varName = ...`) inside. This ensures variables are accessible to downstream nodes outside the try/catch scope.

## Execution Visualization

KNIME-inspired live execution feedback. The script runs in the same browser context and communicates back to the Flow view via custom events.

### Architecture
- **ExecutionController** (`execution-controller.ts`): Orchestrates runs. Subscribes to `funcflow.exec.{runId}` events, updates `ExecutionState`, drives `ExecutionVisualizer`.
- **ExecutionState** (`execution-state.ts`): Tracks per-node status (idle/running/completed/errored/stale) and runtime output summaries.
- **ExecutionVisualizer** (`execution-visualizer.ts`): Maps status → node visual properties (`boxcolor`, `bgcolor`, `onDrawForeground` overlay dot).
- **ValueInspector** (`value-inspector.ts`): Renders runtime output values in the context panel when a completed/errored node is selected.

### Visual States
| State | boxcolor | bgcolor | Overlay |
|---|---|---|---|
| Idle | `#78909c` | `#ffffff` | none |
| Running | `#1976d2` (blue) | `#e3f2fd` | pulsing blue dot (top-left, at collapse icon) |
| Completed | `#43a047` (green) | `#ffffff` | green dot with checkmark (top-left) |
| Errored | `#e53935` (red) | `#ffebee` | red dot with "!" (top-left) |
| Stale | `#9E9E9E` (gray) | `#f5f5f5` | dimmed dot |

The status dot is positioned at the top-left of the title bar (same area as the collapse icon) for consistent visibility. During the running state, a periodic animation timer forces canvas redraws to ensure smooth pulsing even without mouse movement.

### Execution Modes
- **Run**: Instrumented script with live visualization. Breakpoint nodes are skipped.
- **Debug**: Same as Run but breakpoint nodes pause execution via `await new Promise(...)` that resolves when the user clicks "Continue" (fires `funcflow.exec.{runId}.continue` event).
- **Run Script (Classic)**: Opens in Datagrok script editor with no instrumentation.

### Invalidation
Graph structural changes increment a version counter. If the graph changes after a run, all completed/errored nodes become **stale** (dimmed gray, still inspectable). Starting a new run resets all states.

### Breakpoint Node (`breakpoint-node.ts`)
Pass-through node (dynamic in → dynamic out) in the "Debug" category. In debug mode, emits code that fires `breakpoint-hit` event and awaits a `continue` event from the view. In normal run mode, the node is skipped entirely. The compiler treats Breakpoint specially: its output slots resolve to its input expression (pure pass-through), so no `breakpoint` variable is declared — downstream nodes reference the original input variable directly.

## File Format

`.ffjson` files store the full flow state including LiteGraph graph data and FuncFlow metadata.

## UI Architecture: Widgets vs Property Panel

**Nodes have NO inline widgets** (except ConstStringNode). All property editing happens in Datagrok's native context panel (right side) when a node is selected.

### Context Panel Integration
- `grok.shell.windows.showContextPanel = true` enables the native panel
- `grok.shell.o = propertyPanel.root` sets its content on node selection
- Layout is 2-panel: `[leftPanel (FunctionBrowser), canvasContainer]` — no custom right panel

### Property Panel (`property-panel.ts`)
- Uses `ui.accordion()` with `accordion.addPane()` for collapsible sections
- **Connections pane** is collapsed by default
- Reads `node.properties` and creates appropriate editors:
  - `string` properties → `<textarea>` (auto-resizing)
  - `number` properties → `<input type="number">` with step (1 for int, 0.1 for double)
  - `boolean` properties → `<input type="checkbox">`
  - enum properties → `<select>` dropdown
- **Exception**: `ConstStringNode` keeps its inline `text` widget for quick editing
- **Func nodes**: primitive input defaults stored as `_input_${name}` properties, edited in panel
- **Collapse icon**: `NODE_DEFAULT_BOXCOLOR = '#78909c'` (Spotfire-inspired blue-gray)
- **Font**: Roboto via Google Fonts import + LiteGraph canvas font overrides
- **Theme**: Spotfire-inspired with soft shadows, subtle dot-grid background, blue selection highlight, rounded nodes

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
4. **Modifying script emission**: All code generation (clean and instrumented) lives in `script-emitter.ts`. The `EmitOptions.instrumented` flag controls whether try/catch + event code is emitted. Add new step types to both the clean and instrumented paths.
5. **Adding execution visual states**: Modify `execution-visualizer.ts` `STATUS_COLORS` map. Node overlays are drawn via `onDrawForeground` hooks.
