# CLAUDE.md - Flow Package

## Overview

Flow (FuncFlow) is an interactive visual function chain designer for Datagrok. It uses **Rete.js v2** to let users compose Datagrok functions, queries, and scripts into executable JavaScript scripts via a node-based graph editor.

The renderer is React-based (`rete-react-plugin`), but React is contained to the canvas: the rest of the package — view, panels, ribbon, status bar — is plain TypeScript on top of the Datagrok UI helpers, exactly as KetcherSketcher mounts Ketcher inside a `ui.div`.

## Architecture

```
src/
├── package.ts                    # Entry: @grok.decorators.app + @grok.decorators.fileViewer
├── funcflow-view.ts              # DG.ViewBase host: 2-panel layout, ribbon, status bar
├── rete/
│   ├── sockets.ts                # TypedSocket extends ClassicPreset.Socket — DG type → compatibility
│   ├── scheme.ts                 # FlowNode / FlowConnection / FlowScheme
│   ├── node-component.tsx        # React Node + Socket components (rendered to DOM by ReactPlugin)
│   ├── flow-editor.ts            # NodeEditor + AreaPlugin + ConnectionPlugin + ReactPlugin wiring
│   ├── node-factory.ts           # Type registry: createNode(typeName), DG.Func discovery
│   └── nodes/
│       ├── input-nodes.ts        # 13 input types — //input: lines
│       ├── output-nodes.ts       # Table & Value (with auto-type detect)
│       ├── utility-nodes.ts      # 9 helpers + 5 constants (ConstString has inline InputControl)
│       ├── comparison-nodes.ts   # 10 ops
│       ├── breakpoint-node.ts    # Debug pause node
│       └── func-node.ts          # Dynamic node factory per DG.Func, builds pass-through outputs
├── compiler/
│   ├── topological-sort.ts       # Kahn's over editor.getConnections()
│   ├── graph-utils.ts            # Thin shim around FlowEditor
│   ├── graph-compiler.ts         # FlowEditor → CompiledStep[]
│   ├── script-emitter.ts         # CompiledStep[] → JS source (clean + instrumented modes)
│   └── validator.ts              # Pre-compilation checks
├── execution/
│   ├── execution-state.ts        # NodeExecStatus enum, NodeExecState (string node IDs)
│   ├── execution-visualizer.ts   # Sets node.dgStatus → CSS handles all visuals
│   ├── execution-controller.ts   # Run lifecycle, event subscriptions, breakpoints
│   ├── value-inspector.ts        # Context panel runtime-value section
│   └── output-preview.ts         # Bottom-docked output preview tabs
├── panel/
│   ├── function-browser.ts       # Left sidebar catalog
│   └── property-panel.ts         # Side-panel node properties editor
├── serialization/
│   ├── flow-schema.ts            # .ffjson v2 type definitions (Rete-native, no LiteGraph payload)
│   └── flow-serializer.ts        # serialize / deserialize / download
├── types/
│   └── type-map.ts               # DG type → slot color, role color, type compatibility
└── utils/
    └── dart-proxy-utils.ts       # Safe Dart proxy property access
```

## Package Entry (`package.ts`)

Uses `@grok.decorators` (TS decorator API):
- `@grok.decorators.app({name: 'Flow', tags: ['app']})` → `funcflowApp(path?)` returns `FuncFlowView`
- `@grok.decorators.fileViewer({fileViewer: 'ffjson'})` → `viewFuncFlow(file)` opens `.ffjson` files

Both entry points hide the toolbox/help and show the context panel.

## Rete Pipeline

```
NodeEditor (data layer — nodes & connections, signals via addPipe)
  └─ AreaPlugin (DOM container, NodeView/ConnectionView per element, transform)
        ├─ ConnectionPlugin (interactive connection drawing — pointerdown on socket DOM)
        │     └─ ClassicFlow with custom canMakeConnection (TypedSocket compatibility)
        └─ ReactPlugin (renders React components into element AreaPlugin provides)
              └─ ClassicPreset with custom Node + Socket components
```

The four plugins are composed in `FlowEditor` ([rete/flow-editor.ts](src/rete/flow-editor.ts)). All other code in the package consumes the `FlowEditor` public API and never imports Rete plugins directly.

### Why React inside the canvas?

`rete-react-plugin` is the official, well-maintained renderer. Writing a vanilla DOM renderer would mean reimplementing socket position-watching, connection SVG paths, mount/unmount lifecycle, and StrictMode handling. The React surface is contained to one component file ([rete/node-component.tsx](src/rete/node-component.tsx)) — no React state outside it, no JSX in the rest of the package.

## Data Flow

```
User double-clicks function in FunctionBrowser   (or drags DG.Func / FileInfo onto canvas)
  → FuncFlowView.addNodeByType(typeName)
  → createNode(typeName) — instantiates a registered FlowNode subclass
  → flow.addNodeAtCenter(node)
  → User connects nodes by dragging between sockets (ClassicFlow validates types)
  → Run / View Script:
      validateGraph(flow) → compileGraph(flow) → emitScript(flow, settings, options)
      → JavaScript source with //input: //output: annotations
```

## Drag-and-Drop

The canvas container accepts drops via `ui.makeDroppable()`:
- **`DG.FileInfo`** (files) → creates an `OpenFile` function node with `inputValues['fullPath']` pre-set.
- **`DG.Func`** (queries, scripts, functions) → finds the matching registered node type and adds it.

## Node Types

All nodes extend [FlowNode](src/rete/scheme.ts) which extends `ClassicPreset.Node`. Per-node metadata: `dgNodeType` (`'input'` / `'output'` / `'utility'` / `'func'`), `dgOutputType`, `dgFunc`, `dgFuncName`, `dgRole`, `passthroughCount`, `properties`, `inputValues`, `pos`, `dgTypeName`, `description` (KNIME-style annotation rendered under the title; for input/output nodes it also becomes the `[description]` suffix in the generated `//input:` / `//output:` line), optional `dgStatus`.

### Input Nodes ([rete/nodes/input-nodes.ts](src/rete/nodes/input-nodes.ts))

Become `//input:` annotation lines.

| Node | DG Type | Qualifiers (in `node.properties`) |
|------|---------|-----------------------------------|
| Table Input | `dataframe` | — |
| Column Input | `column` | typeFilter, semTypeFilter |
| Column List Input | `column_list` | typeFilter, semTypeFilter |
| String Input | `string` | nullable, choices, caption, semType |
| Number Input | `double` | nullable, min, max, showSlider, caption |
| Int Input | `int` | nullable, min, max, showSlider, caption |
| Boolean Input | `bool` | nullable, caption |
| DateTime Input | `datetime` | nullable, caption |
| File Input | `file` | nullable, caption |
| Map Input | `map` | nullable, caption |
| Dynamic Input | `dynamic` | nullable, caption |
| String List Input | `string_list` | — |
| Blob Input | `blob` (slot type `byte_array`) | nullable, caption |

### Output Nodes ([rete/nodes/output-nodes.ts](src/rete/nodes/output-nodes.ts))

Become `//output:` annotation lines.

| Node | Notes |
|------|-------|
| Table Output | Fixed `dataframe` type |
| Value Output | Configurable `outputType`. On connect, `FlowEditor.maybeAutoTypeValueOutput` copies the source slot's `dgType` into `properties.outputType` (skipping `dynamic` / `object`). |

### Utility Nodes ([rete/nodes/utility-nodes.ts](src/rete/nodes/utility-nodes.ts))

| Node | Generated Code |
|------|----------------|
| Select Column | `let v = df.col('name')` |
| Select Columns | `let v = [df.col('a'), df.col('b')]` |
| Add Table View | `let v = grok.shell.addTableView(df)` |
| Log | `console.log([label,] value)` |
| Info | `grok.shell.info(msg)` |
| Warning | `grok.shell.warning(msg)` |
| ToString | `(value).toString()` |
| FromJSON | `JSON.parse(json)` |
| ToJSON | `JSON.stringify(value)` |
| Constants (String / Int / Double / Boolean / List) | inline literals |

`ConstStringNode` is the only node with an inline widget — a `ClassicPreset.InputControl` for fast text editing. All other property editing happens in the side panel.

### Comparison Nodes ([rete/nodes/comparison-nodes.ts](src/rete/nodes/comparison-nodes.ts))

`==`, `!=`, `>`, `>=`, `<`, `<=`, `Contains`, `Starts With`, `Ends With`, `Is Null`.

### Function Nodes ([rete/nodes/func-node.ts](src/rete/nodes/func-node.ts))

`FuncNode` is built dynamically per `DG.Func`. Header uses `func.friendlyName` (split by `|`, last segment). Color comes from `ROLE_COLORS` map. Generates `await grok.functions.call('Pkg:funcName', {...})`.

#### Pass-Through Outputs

Every func node automatically gets **pass-through output slots** mirroring each input. They solve the execution-ordering problem for mutating functions:

- **Problem**: `addNewColumn(table)` and `doSomething(table)` both consume the same table — topological sort can't determine order without an edge between them.
- **Solution**: Connect `addNewColumn`'s pass-through output to `doSomething`'s input. Compiler treats the pass-through edge as an ordering edge only.
- **Encoding**: Pass-through output keys are `<inputName>__pt`. The visible label is just `→`. `FuncNode.passthroughInputName(key)` extracts the original input name. `node.passthroughCount` records how many pass-through slots are at the start of the outputs map.
- **Compile**: `graph-compiler.ts` resolves a pass-through output to the same expression as the corresponding input — no new variable is generated.
- **Visual**: dashed border on the socket, faded italic label.

## Function Filtering ([rete/node-factory.ts](src/rete/node-factory.ts))

Functions with no inputs *and* no outputs are skipped. Functions whose role appears in `EXCLUDED_ROLES` (or any tag in `EXCLUDED_TAGS`) are skipped.

`registerBuiltinNodes()` populates the `FACTORIES` map with all built-in types. `registerAllFunctions()` discovers DG functions via `DG.Func.find({})` and registers a per-func factory under name `DG Functions/<role>/<funcName>` (or `DG Functions/<role>/<pkg>:<funcName>` on collision).

`createNode(typeName)` looks up the factory and stamps `dgTypeName` on the new instance — this is what the serializer persists.

## Type System (`types/type-map.ts`)

- `DG_TYPE_MAP`: DG type string → `{slotType, color}`. The slot color is what the React Socket component fills the dot with.
- `ROLE_COLORS`: DG role → title-bar color (white body always).
- `areTypesCompatible(out, in)`: source-of-truth for connection validity. Used by `TypedSocket.isCompatibleWith`. Permissive for `dynamic` and `object`; explicit pairs for `int↔double↔num` and `list↔string_list`.

`TypedSocket` ([rete/sockets.ts](src/rete/sockets.ts)) is one-instance-per-DG-type (cached), so reference equality holds. Its `isCompatibleWith` method is consulted at connection-pick time by `ClassicFlow.canMakeConnection`, which rejects incompatible drops before they enter the editor's data layer.

## Script Generation

Pipeline ([compiler/](src/compiler)):

1. **Topological sort** — Kahn's algorithm over `editor.getConnections()`.
2. **Compile** — every node becomes a `CompiledStep` with `inputs: Map<key, expr>`, `outputs: Map<key, varName>`, `properties`, `inputValues`. Variable names: camelCase of node label + first real output; collisions deduplicated by suffix.
3. **Emit** — steps become JS lines; dataframe inputs first; `//input:` / `//output:` headers from properties + qualifiers.

### Validation ([compiler/validator.ts](src/compiler/validator.ts))

- Empty graph → warning
- Cycles → error
- Column input without Table input → error
- Disconnected non-input nodes → warning
- Output node with no incoming connection → warning
- Empty / invalid JS-identifier `paramName` → error
- Duplicate `paramName` across input/output nodes → error

### Input qualifier emission

Compiler reads `node.properties` and emits `{type:..; semType:..; nullable: true; caption: ..; choices: [..]; min:..; max:..; showSlider: true}`. Description appended as `[description]`. Default values inlined as `//input: type name = <value>`.

### Auto Semantic-Type Detection

After every function step, the emitter writes
`if (varName != null) await varName.meta.detectSemanticTypes();`
for each non-pass-through dataframe output (clean *and* instrumented modes).

### Instrumented Mode (`EmitOptions`)

```typescript
interface EmitOptions {
  instrumented?: boolean;     // false = clean script, true = try/catch + events
  runId?: string;             // UUID for this run (required when instrumented)
  enableBreakpoints?: boolean; // emit breakpoint pause code (debug mode)
  haltOnError?: boolean;      // throw on first error (default true)
}
```

Each step is wrapped in try/catch and fires `funcflow.exec.<runId>` events: `run-start`, `node-start`, `node-complete`, `node-error`, `breakpoint-hit`, `run-complete`. Output values are summarized by an inline `__ff_summarize(value, declaredType?)` (DataFrame → row/col + clone, Column → name + 5-element sample, graphics → raw, primitives, etc.).

**Variable hoisting**: when wrapping `let x = ...`, the declaration is hoisted before `try` and only the assignment goes inside, so downstream nodes can reference `x`.

**In-place mutating function support**: when a func node has dataframe input(s) but **zero real outputs**, the wrapper emits a synthetic output entry `'<inputName> (modified)': __ff_summarize(<inputExpr>, 'dataframe')` so the modified table is previewable.

## Execution Visualization

KNIME-inspired live feedback. The script runs in the same browser tab and communicates back via custom events.

### Architecture

- **ExecutionController** ([execution/execution-controller.ts](src/execution/execution-controller.ts)) — orchestrates runs. Subscribes to `funcflow.exec.<runId>`, updates `ExecutionState`, drives `ExecutionVisualizer`, pushes outputs into `OutputPreviewPanel`. Exposes callbacks `onBreakpointHit`, `onRunEnd`, `onNodeStateChanged`.
- **ExecutionState** ([execution/execution-state.ts](src/execution/execution-state.ts)) — per-node status tracking (`idle` / `running` / `completed` / `errored` / `stale`) keyed by string node IDs.
- **ExecutionVisualizer** ([execution/execution-visualizer.ts](src/execution/execution-visualizer.ts)) — sets `node.dgStatus` and calls `flow.updateNode(id)` to re-render. The React Node component reads `dgStatus` and writes it to a `data-status` attribute. CSS does the rest (status circle color, pulse animation, body tint).
- **ValueInspector** ([execution/value-inspector.ts](src/execution/value-inspector.ts)) — runtime-value section in the side panel. DataFrame summaries embed a full `DG.Viewer.grid` preview with "Add to workspace".
- **OutputPreviewPanel** ([execution/output-preview.ts](src/execution/output-preview.ts)) — bottom-docked tabs after a run completes. Classifies outputs (dataframe/viewer/widget/graphics/primitive) using `getOutputTypeHints()`.

### Visual States (CSS, in [css/funcflow.css](css/funcflow.css))

| State | `.ff-node-status` | `.ff-node` body |
|---|---|---|
| idle | white circle + gray outline | white |
| running | blue, `@keyframes ff-pulse` | light blue |
| completed | green + checkmark via `::after` | white |
| errored | red + `!` via `::after` | light red |
| stale | gray, 65% opacity | light gray |

The status circle is part of the title bar; CSS keyframes drive the pulse animation — no JS animation timer needed.

### Execution Modes

- **Run** — instrumented script with live visualization. Breakpoints skipped.
- **Debug** — same as Run but breakpoint nodes pause via `await new Promise(...)` resolved on `funcflow.exec.<runId>.continue`.
- **Run Script (Classic)** — clean (non-instrumented) script run via `DG.Script.create(script).prepare()`, outputs piped into the same `OutputPreviewPanel`.
- **View Script** — opens a dialog with the generated source; buttons: Copy / Export `.js` / Open in ScriptView / Run.

### Invalidation

`onGraphChanged()` increments a graph version. Completed/errored nodes become **stale** when the graph changes after a run.

## File Format

`.ffjson` v2 — Rete-native. A v1 → v2 migrator lives at `tools/migrate-ffjson-v1-to-v2.py` (run with `py tools/migrate-ffjson-v1-to-v2.py path/to/flow.ffjson`); the demo files in `files/` were converted with it. The runtime loader rejects anything other than `version: '2.0'`.

```typescript
{
  version: '2.0',
  name, description, author, created, modified,
  nodes: [{id, typeName, label, pos, properties, inputValues}],
  connections: [{id, source, sourceOutput, target, targetInput}],
  metadata: { settings: { scriptName, scriptDescription, tags } }
}
```

`serializeFlow(flow, settings)` produces the doc; `deserializeFlow(doc, flow)` clears and rebuilds. Unknown `typeName`s are skipped with a console warning. Connections referencing missing nodes are silently dropped.

## UI Architecture

### Layout

```
+--------------------------------------------------------------+
|  Rete canvas container                                       |
|  (AreaPlugin mounts here; FunctionBrowser lives in the       |
|   platform toolbox window — the view sets `this.toolbox =    |
|   functionBrowser.root` and turns `showToolbox` on)          |
+--------------------------------------------------------------+
|  Status bar: Nodes / Links / Validation                      |
+--------------------------------------------------------------+
```

Property panel goes into Datagrok's native context panel via `grok.shell.o = propertyPanel.root`. The function browser is the view's toolbox; toggling `grok.shell.windows.showToolbox` shows/hides it (ribbon icon `list-ul`).

A bottom-docked **Output panel** is *lazy*: never auto-opened. The first time the user clicks a node that has captured runtime values from a prior run, [`ExecutionController.showOutputsForNode`](src/execution/execution-controller.ts) creates the dock at the bottom; subsequent clicks update it in place. Closed when the graph changes or a new run starts.

### Property Panel ([panel/property-panel.ts](src/panel/property-panel.ts))

- Title row at top (editable label) + node-type badge.
- Accordion with type-specific panes:
  - Func nodes: **Function** (description, full name, role) + **Input Parameters** (per-input editor for primitives via `node.inputValues`; "connected only" label otherwise).
  - Input nodes: **Input Configuration** (paramName, description, defaultValue, nullable, caption, type/semType filters, choices, min/max, showSlider).
  - Output nodes: **Output Configuration** (paramName + outputType combo for ValueOutput).
  - Utility nodes: **Configuration** for non-underscore properties (bool/number/text auto-detected).
  - **Connections** pane (collapsed by default): Inputs / Pass-through / Outputs grouped, with connection status from `flow.isInputConnected()` / `flow.getConnections()`.
- Editor helpers (auto-resizing textarea, number, toggle, combo) all support optional tooltips.
- After each property change, calls `flow.updateNode(node.id)` to re-render visible state (label, etc.).

### Function Browser ([panel/function-browser.ts](src/panel/function-browser.ts))

- Search input + Group-by selector (`role` / `tags` / `package` / `output`). The `output` mode buckets functions by what they produce: **Data Sources** (output dataframe), **Transformations** (no declared output — typically in-place mutators), **Column Operations**, **Visualization** (viewer/view/widget/graphics), **Compute** (primitives), **Utilities** (everything else).
- Built-in sections: **Inputs**, **Outputs**, **Constants**, **Comparisons** (collapsed by default), **Utilities**, **Debug**. Tooltips on section headers.
- DG functions follow, alphabetically grouped per the chosen mode.
- Item tooltips: built-ins use `<description>. Double-click to add`; DG functions show `func.description (packageName)`.

### Theme (CSS)

Spotfire-inspired light theme:
- White nodes with colored title bars (per `dgNodeType` or per-role).
- Soft drop shadow (`box-shadow: 0 3px 10px rgba(0,0,0,0.1)`), rounded `border-radius: 8px`.
- Selection: blue 1.5px border + outer halo.
- Background: `#ebedf2` with subtle dot grid (programmatically generated PNG data URL, applied inline by `FlowEditor`).
- Sockets: 12px dot, type-colored fill, white border, 1px gray ring; hover scales 1.18×.
- Pass-through sockets: dashed white border, faded label.
- Connections: cubic Bézier SVG path, 2.5px wide. **Stroke color is set per-connection from the source slot's DG type** — `FlowEditor.styleConnectionElement` runs on every `rendered` signal and writes the color into the `<path>` attribute. This way every connection visually carries its data type.
- Connection states (`data-status="active"`/`"completed"`/`"errored"`/`"stale"`): `active` adds the marching-dashes `@keyframes ff-flow-march` animation to show data flowing through the edge; `errored` overrides the stroke to red; `stale` dims and dashes the path.

### Interaction

- **Collapse / expand**: click the status circle in the title bar (the same dot that shows execution state). Sets `node.collapsed`; the React component re-renders. Collapsed nodes still expose their socket DOM in a hidden absolute-positioned row at the title bar's left/right edges so existing connections keep their endpoints.
- **Context menu (right-click on a node)**: Collapse/Expand · Duplicate · Hide · Delete. Right-click on a connection: Delete connection. Built with `DG.Menu.popup().item(...).separator().show({causedBy: ev})` so the platform handles positioning, dismissal, and styling. Trigger comes from `area-plugin`'s `contextmenu` signal — `data.context` is `'root'` / a `FlowNode` / a `FlowConnection` and we branch on it.
- **Delete key (or Backspace)**: removes every selected node and all connections touching it. The handler is registered on `window.keydown` and skips events whose target is an `<input>`/`<textarea>`/`<select>` so typing in the property panel never deletes nodes.
- **Drag from toolbox**: `ui.makeDroppable` on the canvas container — accepts `DG.FileInfo` (creates an `OpenFile` node with `inputValues['fullPath']`) or `DG.Func` (looks up the registered factory by func name). Native HTML5 drag from the FunctionBrowser carries the typeName via the `FF_DRAG_MIME` data type; the canvas drop handler reads it and calls `addNodeAt` with the drop point.
- **Snap-to-grid + alignment guides**: `nodetranslate` is intercepted in `FlowEditor.wireEvents` — `computeSnap` returns either an alignment-snapped position (when an edge/center is within `alignThreshold` of another node's edge/center) or a 20px grid-snapped position. While dragging, dashed guide lines on a screen-space overlay (`.ff-guide-overlay`) mark the alignment axes; hidden on `nodedragged`.
- **Pan-to-new-node**: `addNodeAtCenter` calls `panToNode(id)` after one rAF, translating the AreaPlugin so the freshly-added node is at viewport center.
- **Drag-out suggestion menu**: pointerdown on `.ff-socket-row-output .ff-socket` records the source slot; if pointerup lands on empty canvas (not a node, not a socket) AND no connection was created in between, [`openSuggestionMenu`](src/rete/flow-editor.ts) opens a searchable popup of every node type whose first input is type-compatible with the dragged source (Value Output sorted first). Selecting an item creates the node at the drop point and auto-connects the source to its first compatible input. Compatible-types lookup is cached per node type — see `findNodeTypesAcceptingInput` / `_sampleInputTypesCache` in `node-factory.ts`.

## Key Dependencies

- `rete` ^2.0.6 — core data layer
- `rete-area-plugin` ^2.1.5 — DOM container, NodeView/ConnectionView, transform
- `rete-connection-plugin` ^2.0.5 — interactive connection drawing + ClassicFlow type validation
- `rete-react-plugin` ^2.1.0 — React renderer with Classic preset
- `rete-render-utils` ^2.0.3 — `getDOMSocketPosition`, `classicConnectionPath`
- `react` / `react-dom` ^18.3 — for the renderer (NOT used outside `node-component.tsx`)
- `styled-components` ^5.3 — peer dep of `rete-react-plugin`
- `datagrok-api` ^1.27.0 — platform API (external, not bundled)

`react`, `react-dom`, `styled-components`, and all `rete-*` packages are in `dependencies` (not `devDependencies`) — they're runtime-required.

## Development Guidelines

1. **Adding a new built-in node**: subclass `FlowNode` in the appropriate `rete/nodes/*.ts` file, set `dgNodeType` and slot definitions in the constructor, then register a factory in `node-factory.ts` `registerBuiltinNodes()`. Add it to the appropriate `*Nodes` array in `function-browser.ts`. If it generates code, add a case to `emitUtilityStep()` in `script-emitter.ts`. Property tooltips go in `UTILITY_PROP_TOOLTIPS` in `property-panel.ts`.
2. **Adding a new DG type**: add to `DG_TYPE_MAP` in `type-map.ts` (slot type + color); extend `COMPATIBLE_TYPES` if needed.
3. **Modifying script emission**: all code generation lives in `script-emitter.ts`. The `EmitOptions.instrumented` flag toggles the try/catch + event path. New step kinds need handling in both clean (`emitFuncStep` / `emitUtilityStep`) and instrumented (`emitFuncStepInstrumented` / `wrapInstrumented`) paths.
4. **Adding execution visual states**: extend `NodeExecStatus` in `execution-state.ts` and add a CSS rule `.ff-node-status[data-status="newstate"] { ... }`. Then update `ExecutionVisualizer` to set the status, no JS animation needed (CSS keyframes handle it).
5. **Adding output preview types**: extend the union in `OutputPreviewPanel` and add a case to `classifyOutput` / `buildPreview`. Make sure `getOutputTypeHints()` propagates the right declared type.
