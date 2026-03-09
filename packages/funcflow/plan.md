# FuncFlow: Interactive Function Chain Designer for Datagrok

## Motivation

Everything in Datagrok is a function (`DG.Func`). Every script, query, viewer, transformation, and application shares the same API: typed inputs (`Property[]`), typed outputs (`Property[]`), tags, metadata (`options`), and a unified execution model (`FuncCall`). This uniformity creates a unique opportunity to build a visual function chain designer — similar in spirit to KNIME, n8n, or Node-RED — but purpose-built for Datagrok's function ecosystem.

FuncFlow leverages the full metadata richness of Datagrok functions (types, semantic types, roles, tags, validation) to provide intelligent type-checking, auto-completion, and one-click script generation. The generated scripts are valid Datagrok scripts ready for immediate use.

---

## Technology Choice: LiteGraph.js

After researching React Flow, Rete.js, JointJS, Drawflow, and custom canvas implementations, **LiteGraph.js** is the recommended library:

| Requirement | LiteGraph.js |
|------------|-------------|
| Vanilla JS/TS (no React) | Yes — zero framework dependencies |
| Canvas-based rendering | Yes — HTML5 Canvas2D |
| Typed port connections | Yes — built-in type validation |
| Inline parameter editors | Yes — widget system (number, string, combo, slider, toggle) |
| Performance (100+ nodes) | Proven at scale (ComfyUI uses it) |
| JSON serialization | Built-in `graph.serialize()` / `graph.configure()` |
| License | MIT |
| Dependencies | Zero runtime dependencies |

LiteGraph provides four layers:
1. **LiteGraph** — global namespace and node type registry
2. **LGraph** — graph container (nodes + connections)
3. **LGraphCanvas** — rendering and interaction (Canvas2D)
4. **LiteGraph.Editor** — optional wrapper with sidebar/search

What LiteGraph lacks (we will implement): minimap overlay, undo/redo (Command Pattern), and Datagrok type mapping.

---

## Architecture Overview

```
FuncFlow Package
├── src/
│   ├── package.ts                    # Package entry, registers the app
│   ├── package.g.ts                  # Auto-generated (grok api)
│   ├── package-test.ts               # Test infrastructure
│   ├── funcflow-view.ts              # Main view (DG.ViewBase)
│   ├── canvas/
│   │   ├── graph-manager.ts          # LGraph wrapper, lifecycle management
│   │   ├── canvas-controller.ts      # LGraphCanvas setup, events, resize
│   │   └── minimap.ts                # Minimap overlay canvas
│   ├── nodes/
│   │   ├── node-factory.ts           # Creates LiteGraph nodes from DG.Func
│   │   ├── func-node.ts              # Base node wrapping a DG.Func
│   │   ├── input-nodes.ts            # Primitive input nodes (table, column, string, etc.)
│   │   ├── output-nodes.ts           # Output/result nodes
│   │   └── utility-nodes.ts          # Utility nodes (select column, etc.)
│   ├── panel/
│   │   ├── function-browser.ts       # Left sidebar: searchable function catalog
│   │   └── property-panel.ts         # Right panel: selected node properties
│   ├── compiler/
│   │   ├── graph-compiler.ts         # Compiles graph → DAG → script
│   │   ├── topological-sort.ts       # Topological ordering of nodes
│   │   ├── script-emitter.ts         # Generates Datagrok script text
│   │   └── validator.ts              # Graph validation with error messages
│   ├── serialization/
│   │   ├── flow-serializer.ts        # Save/load .funcflow.json files
│   │   └── flow-schema.ts            # JSON schema types
│   ├── types/
│   │   ├── type-map.ts               # DG.TYPE ↔ LiteGraph type mapping + colors
│   │   └── type-colors.ts            # Color palette for types
│   ├── history/
│   │   └── undo-manager.ts           # Command Pattern undo/redo
│   └── utils/
│       └── dart-proxy-utils.ts       # Safe access for Dart proxy objects (tags, options)
├── css/
│   └── funcflow.css                  # Custom styles
├── files/                            # Sample .funcflow.json files
├── package.json
├── webpack.config.js
└── tsconfig.json
```

---

## Detailed Implementation Plan

### Phase 1: Foundation — Canvas, Nodes, and Connections

#### Step 1.1: Install LiteGraph.js

```bash
cd packages/funcflow
npm install litegraph.js
```

Add to `webpack.config.js` (it should be bundled, not external):
```js
// Do NOT add litegraph.js to externals — it should be bundled
```

#### Step 1.2: Create `type-map.ts` — Datagrok ↔ LiteGraph Type Mapping

**File:** `src/types/type-map.ts`

Map Datagrok `TYPE` enum values to LiteGraph slot types with assigned colors:

| DG Type | LiteGraph Slot Type | Color | Hex |
|---------|-------------------|-------|-----|
| `TYPE.DATA_FRAME` | `"dataframe"` | Orange | `#E67E22` |
| `TYPE.COLUMN` | `"column"` | Blue | `#3498DB` |
| `TYPE.COLUMN_LIST` | `"column_list"` | Light Blue | `#5DADE2` |
| `TYPE.STRING` | `"string"` | Green | `#2ECC71` |
| `TYPE.INT` | `"int"` | Teal | `#1ABC9C` |
| `TYPE.FLOAT` | `"double"` | Cyan | `#00BCD4` |
| `TYPE.BOOL` | `"bool"` | Red | `#E74C3C` |
| `TYPE.DATE_TIME` | `"datetime"` | Purple | `#9B59B6` |
| `TYPE.STRING_LIST` | `"string_list"` | Lime | `#8BC34A` |
| `TYPE.OBJECT` | `"object"` | Gray | `#95A5A6` |
| `TYPE.DYNAMIC` | `"dynamic"` | White | `#ECF0F1` |
| `TYPE.MAP` | `"map"` | Amber | `#FFC107` |
| `TYPE.FUNC_CALL` | `"funccall"` | Pink | `#E91E63` |

The type map also defines compatibility rules (e.g., `int` connects to `double`, `dynamic` connects to anything).

Register all types with `LiteGraph.registerSlotType()` for color-coded ports.

#### Step 1.3: Create `func-node.ts` — DG.Func → LiteGraph Node

**File:** `src/nodes/func-node.ts`

Each `DG.Func` becomes a LiteGraph node type. For a function like:
```
//name: Sin
//input: double x
//output: double result
```

The node factory creates:
- **Title:** `Sin` (from `func.name`)
- **Input slots:** One input port `x` with type `"double"`
- **Output slots:** One output port `result` with type `"double"`
- **Color:** Based on `func.options['meta.role']` — see role color table below
- **Inline widgets:** For primitive inputs (`int`, `double`, `string`, `bool`), add LiteGraph widgets directly on the node

**Role Colors:**

| meta.role | Color | Hex |
|-----------|-------|-----|
| `app` | Indigo | `#3F51B5` |
| `panel` | Deep Purple | `#673AB7` |
| `viewer` | Blue | `#2196F3` |
| `transform` / `Transform` | Teal | `#009688` |
| `filter` | Amber | `#FFC107` |
| `converter` | Orange | `#FF9800` |
| `widget` | Pink | `#E91E63` |
| `cellRenderer` | Brown | `#795548` |
| `semTypeDetector` | Lime | `#CDDC39` |
| `fileViewer` / `fileExporter` | Cyan | `#00BCD4` |
| `editor` | Light Blue | `#03A9F4` |
| (no role / absent) | Gray | `#9E9E9E` |

**Key implementation details:**
- Access `func.options` and `func.tags` using `Object.entries()` (they may be Dart object proxies)
- Store the original `DG.Func` reference on the node: `this.dgFunc = func;`
- Store parameter metadata: `this.dgInputProps = func.inputs;` `this.dgOutputProps = func.outputs;`
- For primitive inputs, when the user edits a widget value, store it as a "hardcoded" value on that input. If a wire is connected, the widget is hidden and the wire value takes precedence.

#### Step 1.4: Create Input Nodes

**File:** `src/nodes/input-nodes.ts`

These are "source" nodes that provide external inputs to the chain. They will appear as script header parameters (`//input: ...`) in the generated script.

| Node Name | Type | Outputs | Widget |
|-----------|------|---------|--------|
| `Table Input` | Input | `table: dataframe` | Text field for parameter name |
| `Column Input` | Input | `column: column` | Text field for name + optional type filter |
| `String Input` | Input | `value: string` | Text field with default value |
| `Number Input` | Input | `value: double` | Number field with default value |
| `Int Input` | Input | `value: int` | Int field with default value |
| `Boolean Input` | Input | `value: bool` | Toggle with default value |
| `String List Input` | Input | `value: string_list` | Text area (comma-separated) |

Each input node stores:
- `paramName: string` — the parameter name for the generated script
- `defaultValue: any` — optional default value
- `description: string` — optional description

These are registered as LiteGraph node types under the `"Inputs"` category.

#### Step 1.5: Create Output Nodes

**File:** `src/nodes/output-nodes.ts`

| Node Name | Type | Inputs | Purpose |
|-----------|------|--------|---------|
| `Table Output` | Output | `table: dataframe` | Marks a result table |
| `Value Output` | Output | `value: dynamic` | Marks a result value |

Output nodes generate `//output:` lines in the script header.

#### Step 1.6: Create Utility Nodes

**File:** `src/nodes/utility-nodes.ts`

| Node Name | Inputs | Outputs | Description |
|-----------|--------|---------|-------------|
| `Select Column` | `table: dataframe`, widget: `colName: string` | `column: column` | Gets column from table by name |
| `Select Columns` | `table: dataframe`, widget: `colNames: string` | `columns: column_list` | Gets multiple columns |
| `Create Table` | (none — add rows/columns via widgets) | `table: dataframe` | Creates a table from hardcoded data |
| `Log` | `value: dynamic` | (none) | Logs value with `grok.shell.info()` |

These are registered under the `"Utilities"` category.

#### Step 1.7: Create `node-factory.ts` — Batch Registration

**File:** `src/nodes/node-factory.ts`

```typescript
export class NodeFactory {
  // Register all DG.Func as LiteGraph node types
  static async registerAllFunctions(): Promise<void> {
    const allFuncs = DG.Func.find({});
    for (const func of allFuncs) {
      const role = getRole(func);        // safely access meta.role
      const tags = getTags(func);        // safely access tags
      const category = role || 'Uncategorized';

      // Register as: "DG Functions/<role>/<funcName>"
      const typeName = `DG Functions/${category}/${func.name}`;
      registerFuncNode(typeName, func);
    }
  }

  // Register built-in input/output/utility nodes
  static registerBuiltinNodes(): void {
    registerInputNodes();   // Under "Inputs/"
    registerOutputNodes();  // Under "Outputs/"
    registerUtilityNodes(); // Under "Utilities/"
  }
}
```

**Safe Dart proxy access pattern:**
```typescript
function getRole(func: DG.Func): string | null {
  try {
    const opts = func.options;
    if (!opts) return null;
    // Use Object.entries or direct access — handle Dart proxy
    const entries = Object.entries(opts);
    for (const [key, val] of entries) {
      if (key === 'meta.role') return String(val);
    }
    return null;
  } catch { return null; }
}

function getTags(func: DG.Func): string[] {
  try {
    const tags = func.tags;
    if (!tags) return [];
    if (Array.isArray(tags)) return tags;
    return Object.values(tags).map(String);
  } catch { return []; }
}
```

---

### Phase 2: View, Panels, and Search

#### Step 2.1: Create `funcflow-view.ts` — Main View

**File:** `src/funcflow-view.ts`

```typescript
//name: FuncFlow
//description: Interactive function chain designer
//tags: app
//meta.role: app
export function funcflowApp(): DG.ViewBase {
  return new FuncFlowView();
}
```

The view extends `DG.ViewBase` with this layout:

```
┌──────────────────────────────────────────────────────────────────┐
│ [Ribbon Menu: File | Edit | View | Generate]  [Ribbon Buttons]  │
├──────────┬───────────────────────────────────────┬───────────────┤
│          │                                       │               │
│ Function │      Canvas (LGraphCanvas)            │  Properties   │
│ Browser  │                                       │  Panel        │
│          │      ┌──────┐     ┌──────┐            │               │
│ [Search] │      │ Node ├────►│ Node │            │  (selected    │
│          │      └──┬───┘     └──────┘            │   node's      │
│ ▸ Inputs │         │                             │   parameters) │
│ ▸ Outputs│      ┌──┴───┐                         │               │
│ ▸ Utils  │      │ Node │                         │               │
│ ▸ Transform│    └──────┘                         │               │
│ ▸ Viewer │                                       │               │
│ ▸ Panel  │                          [Minimap]    │               │
│ ...      │                                       │               │
├──────────┴───────────────────────────────────────┴───────────────┤
│ [Status Bar: node count | connection count | validation status] │
└──────────────────────────────────────────────────────────────────┘
```

**Layout implementation:**
- **Root:** `ui.splitH([leftPanel, centerCanvas, rightPanel])` with ratios `[0.2, 0.6, 0.2]`
- **Left panel** (`toolbox`): Function browser with accordion + search
- **Center:** Canvas element managed by `LGraphCanvas`
- **Right panel:** Property panel for selected node (via Datagrok accordion)
- **Ribbon menu:** File (New, Open, Save, Export), Edit (Undo, Redo, Delete), View (Zoom, Fit, Minimap), Generate (Script, Copy, Export .js)
- **Status bar:** Node count, connection count, validation warnings

**Canvas resize handling:**
```typescript
ui.onSizeChanged(this.root).subscribe(() => {
  this.canvas.width = this.canvasContainer.clientWidth;
  this.canvas.height = this.canvasContainer.clientHeight;
  this.graphCanvas.resize();
});
```

#### Step 2.2: Create `function-browser.ts` — Left Sidebar

**File:** `src/panel/function-browser.ts`

The function browser provides:
1. **Search bar** at the top — filters by name, description, tags, and meta.role
2. **Grouping mode** toggle — group by: `meta.role` (default), `tags`, `package`
3. **Tree view** — expandable categories, each containing function items
4. **Drag support** — drag a function from the browser onto the canvas to create a node

```typescript
class FunctionBrowser {
  private searchInput: HTMLInputElement;
  private treeContainer: HTMLElement;
  private allFuncs: DG.Func[];
  private groupBy: 'role' | 'tags' | 'package' = 'role';

  async init(): Promise<void> {
    this.allFuncs = DG.Func.find({});
    // Also add built-in nodes (Inputs, Outputs, Utilities)
    this.render();
  }

  private render(): void {
    const tree = ui.tree();
    const grouped = this.groupFunctions(this.filterBySearch());

    for (const [category, funcs] of Object.entries(grouped)) {
      const group = tree.group(category);
      for (const func of funcs) {
        const item = group.item(func.name);
        // Make draggable — on drop, create node on canvas
        ui.makeDraggable(item.root, {
          getDragObject: () => func,
          getDragCaption: () => func.name
        });
      }
    }

    this.treeContainer.innerHTML = '';
    this.treeContainer.appendChild(tree.root);
  }

  private filterBySearch(): DG.Func[] {
    const query = this.searchInput.value.toLowerCase();
    if (!query) return this.allFuncs;
    return this.allFuncs.filter(f => {
      return f.name.toLowerCase().includes(query)
        || (f.description || '').toLowerCase().includes(query)
        || getTags(f).some(t => t.toLowerCase().includes(query))
        || (getRole(f) || '').toLowerCase().includes(query);
    });
  }

  private groupFunctions(funcs: DG.Func[]): Record<string, DG.Func[]> {
    // Group by selected mode (role, tags, or package)
  }
}
```

#### Step 2.3: Create `property-panel.ts` — Right Sidebar

**File:** `src/panel/property-panel.ts`

When a node is selected on the canvas, show its properties:
- For **DG.Func nodes**: show all input parameters with appropriate Datagrok input controls
- For **Input nodes**: show parameter name, default value, description
- For **Output nodes**: show output name, type
- For **any node**: show node title (editable), description, connected inputs/outputs summary

Use Datagrok's `ui.accordion()` with categories for organizing parameters:
```typescript
const acc = ui.accordion();
acc.addPane('Parameters', () => this.renderParams(node));
acc.addPane('Connections', () => this.renderConnections(node));
acc.addPane('Info', () => this.renderInfo(node));
```

---

### Phase 3: Graph Compilation and Script Generation

#### Step 3.1: Create `topological-sort.ts`

**File:** `src/compiler/topological-sort.ts`

Performs topological sort on the node graph to determine execution order. Uses Kahn's algorithm:
1. Find all nodes with no incoming connections (source nodes)
2. Process them in order, removing edges and adding newly available nodes
3. Detect cycles (error if any nodes remain unprocessed)

Returns an ordered list of node IDs representing valid execution order.

#### Step 3.2: Create `validator.ts`

**File:** `src/compiler/validator.ts`

Validates the graph before compilation. Returns `ValidationResult[]` with severity and message:

| Validation Rule | Severity | Message |
|----------------|----------|---------|
| Column Input without Table Input | Error | `"Column input '${name}' requires a Table input connected upstream"` |
| Unconnected required input | Error | `"Required input '${param}' on node '${node}' is not connected and has no default value"` |
| Type mismatch on connection | Error | `"Cannot connect '${outType}' to '${inType}' — incompatible types"` |
| Cycle detected | Error | `"Cycle detected in graph: ${cycleNodes}"` |
| Multiple connections to single input | Error | `"Input '${param}' on '${node}' has multiple connections — only one allowed"` |
| Disconnected node (no inputs and no outputs connected) | Warning | `"Node '${name}' is disconnected from the graph"` |
| Empty graph | Warning | `"Graph is empty — nothing to generate"` |
| Output node without connection | Warning | `"Output node '${name}' has no incoming connection"` |

#### Step 3.3: Create `graph-compiler.ts`

**File:** `src/compiler/graph-compiler.ts`

Compiles the validated, topologically-sorted graph into executable steps:

```typescript
interface CompiledStep {
  nodeId: number;
  funcName: string;            // Full qualified name: "Package:FuncName"
  variableName: string;        // Variable name for the result
  inputs: Map<string, string>; // paramName → source expression (variable name or literal)
  outputs: Map<string, string>;// paramName → variable name
  isAsync: boolean;            // Whether to use await
}
```

**Variable naming strategy:**
- Each node's output gets a variable name derived from: `funcName_outputName` (camelCase, deduped with suffix number if needed)
- Example: `sin_result`, `pearson_corr`, `tableInput_table`

#### Step 3.4: Create `script-emitter.ts`

**File:** `src/compiler/script-emitter.ts`

Generates the final Datagrok script text from compiled steps.

**Generated script format example:**

For a graph like: `[Table Input] → [Sin (on column)] → [Table Output]`

```javascript
//name: MyFuncFlow
//description: Generated by FuncFlow
//language: javascript
//input: dataframe inputTable {caption: Input Table}
//output: dataframe result

let col = inputTable.col('values');
let sinResult = await grok.functions.call('Sin', {x: col});
result = sinResult;
```

**Script generation rules:**
1. **Header comments:** Input nodes → `//input:` lines. Output nodes → `//output:` lines.
2. **Parameter qualifiers:** Preserve type, caption, default value, optional flag from the input node configuration.
3. **Function calls:** Use `await grok.functions.call('Package:FuncName', { paramName: value, ... })` for each function node. Functions are always called async with `await`.
4. **Variable references:** When a node's output feeds into another node's input, use the variable name.
5. **Hardcoded values:** When a primitive input has a widget value (no wire connected), emit the literal value directly in the function call parameters.
6. **Column access:** For `Select Column` utility nodes, emit `table.col('columnName')`.

**Full example — more complex chain:**

```
[Table Input: "df"]
   ├──→ [Select Column: "age"] ──→ [Sin] ──→ [Value Output: "sinAge"]
   └──→ [Select Column: "weight"] ──→ [Pearson (with sinAge)] ──→ [Value Output: "corr"]
```

Generated script:
```javascript
//name: MyAnalysis
//language: javascript
//input: dataframe df
//output: double sinAge
//output: double corr

let age = df.col('age');
let weight = df.col('weight');
let sinAge = await grok.functions.call('Sin', {x: age});
let corr = await grok.functions.call('Pearson', {c1: sinAge, c2: weight});
```

---

### Phase 4: Save/Load Format

#### Step 4.1: Define JSON Schema

**File:** `src/serialization/flow-schema.ts`

```typescript
interface FuncFlowDocument {
  version: '1.0';
  name: string;
  description: string;
  author: string;
  created: string;           // ISO date
  modified: string;          // ISO date

  // LiteGraph's native serialization (contains all node/connection data)
  graph: LGraphSerializedData;

  // FuncFlow-specific metadata overlay
  metadata: {
    // Map of nodeId → additional metadata not in LiteGraph
    nodes: Record<number, FuncFlowNodeMeta>;
    // Global flow settings
    settings: FlowSettings;
  };
}

interface FuncFlowNodeMeta {
  dgFuncName: string;        // Full qualified function name (e.g., "Chem:SimilaritySearch")
  dgFuncType: string;        // "func" | "input" | "output" | "utility"
  paramName?: string;         // For input/output nodes: the script parameter name
  defaultValue?: any;         // For input nodes: default value
  description?: string;       // For input/output nodes: parameter description
  paramQualifiers?: Record<string, string>; // Additional qualifiers (caption, units, etc.)
}

interface FlowSettings {
  scriptName: string;
  scriptDescription: string;
  scriptLanguage: string;     // Always "javascript" for now
  tags: string[];
}

// LiteGraph's built-in serialization captures:
// - All node positions, sizes, types
// - All connections (origin_id, origin_slot, target_id, target_slot)
// - All widget values
// - Node properties
```

**File extension:** `.funcflow.json`

#### Step 4.2: Implement Save/Load

**File:** `src/serialization/flow-serializer.ts`

```typescript
class FlowSerializer {
  // Save graph to JSON
  static serialize(graph: LGraph, metadata: FlowMetadata): FuncFlowDocument {
    return {
      version: '1.0',
      name: metadata.scriptName,
      description: metadata.scriptDescription,
      author: grok.shell.user?.login ?? 'unknown',
      created: new Date().toISOString(),
      modified: new Date().toISOString(),
      graph: graph.serialize(),    // LiteGraph built-in
      metadata: buildMetadata(graph),
    };
  }

  // Load graph from JSON
  static deserialize(doc: FuncFlowDocument, graph: LGraph): void {
    graph.configure(doc.graph);    // LiteGraph built-in
    applyMetadata(graph, doc.metadata);
  }
}
```

Save/load triggers in the ribbon menu:
- **Save:** Opens a save dialog (file name input), then writes JSON via `DG.Utils.download(fileName, jsonString)`
- **Open:** File input dialog, reads JSON, deserializes into current graph

---

### Phase 5: Ribbon Menu and Export

#### Step 5.1: Ribbon Menu Structure

```
File
  ├── New Flow          (Ctrl+N)     — Clear canvas, start fresh
  ├── Open Flow...      (Ctrl+O)     — Load .funcflow.json file
  ├── Save Flow         (Ctrl+S)     — Save current flow as .funcflow.json
  └── Save Flow As...   (Ctrl+Shift+S)

Edit
  ├── Undo              (Ctrl+Z)
  ├── Redo              (Ctrl+Y)
  ├── Delete Selected   (Delete)
  ├── Select All        (Ctrl+A)
  └── Duplicate         (Ctrl+D)

View
  ├── Zoom to Fit       (Home)
  ├── Zoom In           (Ctrl+Plus)
  ├── Zoom Out          (Ctrl+Minus)
  ├── Toggle Minimap
  └── Toggle Grid

Generate
  ├── Generate Script          — Validate + compile + show in dialog
  ├── Copy Script to Clipboard — One-click copy
  ├── Export as .js File       — Download as .js file
  └── Validate Graph           — Run validation, show errors/warnings
```

#### Step 5.2: Ribbon Button Panels

```typescript
v.setRibbonPanels([
  [
    ui.iconFA('play', () => this.generateAndPreview(), 'Generate Script'),
    ui.iconFA('copy', () => this.copyScriptToClipboard(), 'Copy Script'),
    ui.iconFA('download', () => this.exportAsJs(), 'Export .js'),
  ],
  [
    ui.iconFA('undo', () => this.undoManager.undo(), 'Undo'),
    ui.iconFA('redo', () => this.undoManager.redo(), 'Redo'),
  ],
  [
    ui.iconFA('search-plus', () => this.zoomIn(), 'Zoom In'),
    ui.iconFA('search-minus', () => this.zoomOut(), 'Zoom Out'),
    ui.iconFA('compress-arrows-alt', () => this.zoomToFit(), 'Zoom to Fit'),
  ]
]);
```

#### Step 5.3: Script Preview Dialog

When "Generate Script" is clicked:
1. Run `validator.validate(graph)` — show errors if any
2. Run `graphCompiler.compile(graph)` → `scriptEmitter.emit(steps)`
3. Show dialog with code editor (read-only) displaying the generated script
4. Dialog has buttons: "Copy to Clipboard", "Export as .js", "Close"

```typescript
const dialog = ui.dialog({title: 'Generated Script'})
  .add(ui.input.code('', {value: generatedScript, language: 'javascript'}))
  .addButton('Copy to Clipboard', () => {
    navigator.clipboard.writeText(generatedScript);
    grok.shell.info('Script copied to clipboard');
  })
  .addButton('Export .js', () => {
    DG.Utils.download('funcflow-script.js', generatedScript);
  })
  .show({width: 700, height: 500});
```

---

### Phase 6: Undo/Redo and Minimap

#### Step 6.1: Undo Manager (Command Pattern)

**File:** `src/history/undo-manager.ts`

```typescript
interface Command {
  execute(): void;
  undo(): void;
  description: string;
}

class UndoManager {
  private undoStack: Command[] = [];
  private redoStack: Command[] = [];

  execute(command: Command): void {
    command.execute();
    this.undoStack.push(command);
    this.redoStack = [];  // Clear redo on new action
  }

  undo(): void { /* pop from undo, call undo(), push to redo */ }
  redo(): void { /* pop from redo, call execute(), push to undo */ }
}
```

Tracked commands:
- `AddNodeCommand` — add/remove node
- `RemoveNodeCommand` — remove/re-add node
- `MoveNodeCommand` — store old/new position (coalesce during drag)
- `ConnectCommand` — add/remove connection
- `DisconnectCommand` — remove/re-add connection
- `ChangeWidgetValueCommand` — old/new widget value

Hook into LiteGraph events:
```typescript
graph.onNodeAdded = (node) => { /* wrap in AddNodeCommand */ };
graph.onNodeRemoved = (node) => { /* wrap in RemoveNodeCommand */ };
graphCanvas.onNodeMoved = (node) => { /* wrap in MoveNodeCommand */ };
graph.onLinkAdded = (link) => { /* wrap in ConnectCommand */ };
graph.onLinkRemoved = (link) => { /* wrap in DisconnectCommand */ };
```

#### Step 6.2: Minimap

**File:** `src/canvas/minimap.ts`

A small overlay canvas (200×150px) in the bottom-right corner:
- Renders all nodes as small colored rectangles
- Shows current viewport as a semi-transparent white rectangle
- Click/drag on minimap to pan the main canvas
- Toggle visibility from View menu

---

### Phase 7: Connection Rules and Type Checking

#### Step 7.1: Connection Validation

LiteGraph supports connection validation via `LiteGraph.isValidConnection(typeA, typeB)` override.

Rules:
1. **Same type** → always valid
2. **`dynamic`** → connects to anything
3. **`int` → `double`** → valid (numeric widening)
4. **`column` → `dynamic`** → valid
5. **`dataframe` → `column`** → invalid (use Select Column node)
6. **Everything else** → must match exactly

Visual feedback: when dragging a connection, compatible ports glow, incompatible ones dim.

#### Step 7.2: Single-Input Rule

Override `LGraphCanvas.onConnectionCreated` to enforce: if an input port already has a connection and a new one is attempted, the old one is removed first (replacement behavior, not multi-connect).

Output ports can have unlimited connections (fan-out).

---

## Key API Patterns to Use

### Finding All Functions
```typescript
const allFuncs: DG.Func[] = DG.Func.find({});
```

### Calling a Function (generated script pattern)
```typescript
const result = await grok.functions.call('PackageName:FuncName', {
  param1: value1,
  param2: value2
});
```

### Getting Function Metadata
```typescript
const func = DG.Func.find({name: 'Sin'})[0];
const inputs: DG.Property[] = func.inputs;    // Input parameter definitions
const outputs: DG.Property[] = func.outputs;   // Output parameter definitions
const role = func.options['meta.role'];         // Function role
const tags = func.tags;                         // Tags (may be Dart proxy)
const desc = func.description;                  // Description
const pkg = func.package;                       // Package it belongs to
```

### Script Header Format (what we generate)
```javascript
//name: ScriptName
//description: Description
//language: javascript
//tags: funcflow, generated
//input: dataframe df {caption: Input Table}
//input: column col {type: numerical}
//input: string name = "default" {optional: true}
//output: dataframe result
//output: double value
```

### Input Parameter Qualifiers (for generated headers)
```
{caption: Display Name}
{type: numerical}                  // column type filter
{optional: true}                   // nullable parameter
{choices: ["a", "b", "c"]}        // string choices
{min: 0; max: 100}                // numeric range
{units: meters}                   // units
```

---

## Implementation Order

| Phase | Tasks | Est. Complexity |
|-------|-------|-----------------|
| **1** | Install LiteGraph, type map, func-node, input/output/utility nodes, node factory | High |
| **2** | FuncFlowView, function browser with search/grouping, property panel, canvas resize | High |
| **3** | Topological sort, validator, graph compiler, script emitter | High |
| **4** | JSON save/load format, serializer | Medium |
| **5** | Ribbon menu, export buttons, script preview dialog | Medium |
| **6** | Undo/redo manager, minimap | Medium |
| **7** | Connection validation, type checking, single-input rule | Low |

**Recommended implementation sequence:** Phase 1 → Phase 2 → Phase 7 → Phase 3 → Phase 4 → Phase 5 → Phase 6

Start with the canvas and nodes (Phase 1), then the view layout and browser (Phase 2), then connection rules (Phase 7) since they affect the graph structure, then compilation (Phase 3), then persistence (Phase 4), then UI polish (Phase 5 + 6).

---

## Files to Reference During Implementation

| Purpose | File Path |
|---------|-----------|
| DG.Func class | `js-api/src/entities/func.ts` |
| Property class | `js-api/src/entities/property.ts` |
| FuncCall class | `js-api/src/functions.ts` |
| TYPE enum | `js-api/src/const.ts` |
| View class | `js-api/src/views/view.ts` |
| ui.* functions | `js-api/src/ui.ts`, `js-api/src/widgets/` |
| Menu class | `js-api/src/widgets/menu.ts` |
| Accordion, Tree | `js-api/src/widgets/containers.ts`, `js-api/src/ui/tree-view.ts` |
| Function calling samples | `packages/ApiSamples/scripts/functions/calling.js` |
| Function discovery | `packages/ApiSamples/scripts/functions/discover-functions.js` |
| FuncCall usage | `packages/ApiSamples/scripts/functions/func-call.js` |
| Registration | `packages/ApiSamples/scripts/functions/register-function.js` |
| Ribbon example | `packages/ApiSamples/scripts/ui/views/ribbon.js` |
| Toolbox example | `packages/ApiSamples/scripts/ui/views/toolbox.js` |
| Drag and drop | `packages/ApiSamples/scripts/ui/interactivity/drag-and-drop.js` |
| Tree view | `packages/ApiSamples/scripts/ui/components/tree-view.js` |
| Script samples | `packages/ApiSamples/scripts/scripting/` |
| Canvas viewer example | `packages/GIS/src/gis-viewer.ts` |
