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
│   ├── graph-layout.ts           # Shared layered/banded layout (importer + Clean Layout ribbon)
│   └── nodes/
│       ├── input-nodes.ts        # 13 input types — //input: lines
│       ├── output-nodes.ts       # Table & Value (with auto-type detect)
│       ├── utility-nodes.ts      # 9 helpers + 5 constants (ConstString has inline InputControl)
│       ├── comparison-nodes.ts   # 10 ops
│       ├── breakpoint-node.ts    # Debug pause node
│       └── func-node.ts          # Dynamic node factory per DG.Func, builds pass-through outputs
├── compiler/
│   ├── topological-sort.ts       # Kahn's, component-by-component (top-y first), y/x-stable
│   ├── graph-utils.ts            # Thin shim around FlowEditor
│   ├── graph-compiler.ts         # FlowEditor → CompiledStep[]
│   ├── script-emitter.ts         # CompiledStep[] → JS source (clean + instrumented modes)
│   ├── creation-script-emitter.ts # FlowEditor → grok-language creation script (inverse of the importer)
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
├── import/
│   └── creation-script-importer.ts  # Table-creation script → flow graph (reverse of the compiler)
├── types/
│   └── type-map.ts               # DG type → slot color, role color, type compatibility
└── utils/
    └── dart-proxy-utils.ts       # Safe Dart proxy property access
```

## Package Entry (`package.ts`)

Uses `@grok.decorators` (TS decorator API):
- `@grok.decorators.app({name: 'Flow', tags: ['app']})` → `funcflowApp(path?)` returns `FuncFlowView`
- `@grok.decorators.fileViewer({fileViewer: 'ffjson'})` → `viewFuncFlow(file)` opens `.ffjson` files
- `@grok.decorators.func(...)` → `flowFromCreationScript(script)` builds a flow from a table-creation
  script, adds the view to the shell, and returns it

Both entry points hide the toolbox/help and show the context panel.

## Creation-Script Import ([import/creation-script-importer.ts](src/import/creation-script-importer.ts))

The reverse of the compiler: takes a **table-creation script** (the cascade of function calls Datagrok
records for reproducibly-created tables — `df.getTag(DG.Tags.CreationScript)`, used by project data
sync) and rebuilds it as a flow graph.

```
Mol1K = OpenFile("System:AppData/Chem/mol1K.csv") //{"timestamp": …}
Chem:addChemPropertiesColumns(Mol1K, "molecule", true, …)
AddNewColumn(Mol1K, "${HBA}+${HBD}+${LogP}", "sumOfSome", subscribeOnChanges = true)
```

- Each line parses via `grok.functions.parse(line, false)` into a `DG.FuncCall` (falls back to
  stripping the trailing `//{…}` metadata comment, quote-aware so URLs survive).
- **Assignments** are `SetVar`-shaped calls (two inputs: string `variableName` + non-null `value`);
  the value call's primary output is recorded in a variable table. The *first* variable is the
  script's result and gets wired to a Table/Value Output node at the end. Assignments do **not**
  advance their inputs (they produce a new variable, they don't mutate inputs).
- **Plain (bare) calls** become `FuncNode`s (via `ensureFuncNodeType` in `node-factory.ts`, which
  registers catalog-excluded funcs on the fly). They are treated as in-place mutators.
- **Inputs by type**: a primitive on an *editable* slot (`string`/`int`/`double`/`num`/`bool`) goes
  into `node.inputValues` (shown/edited in the panel). A primitive on any other slot — notably the
  `dynamic` `value` of `ResolveColumn` (a column name) — gets a **Constant node** wired in, because
  those slots aren't editable in the panel and need an explicit producer. `FuncCall` inputs connect
  recursively; `GetVar` inputs resolve through the variable table.
- **Column arguments** parse to `ResolveColumn(value:dynamic, parentTable:dataframe)`. The platform
  `ResolveColumn` misbehaves at runtime, *and* a Select Column node per column argument clutters the
  graph — so the importer **inlines the column name into `node.inputValues[paramName]`** (the
  `column`/`column_list` slot is editable in the panel, seeded by `FuncNode`). `ResolveColumnList` →
  a **comma-separated** value string. At compile time the value becomes `table.col('name')` /
  `[table.col('a'), …]` (see the compiler note below). These never advance variables.
  - **Which table** a column resolves against is stored explicitly in the func node's
    `properties['columnTables']` (`{columnParam → dataframeParam}`), seeded by `FuncNode` to the
    dataframe input sharing the param's numeric suffix (`keys2`→`table2`) else the first dataframe
    input (`defaultTableParam`). The panel shows a table-choice combo when the func has **≥2**
    dataframe inputs (unambiguous for one). The compiler reads this association (falling back to the
    same suffix/first heuristic for older saves).
  - **Marginal case**: a column input on a func with **no** dataframe input can't resolve a table —
    `FuncNode` doesn't seed it (stays connection-only) and the importer falls back to a real
    **Select Column / Select Columns** node (the older behavior, still wired by `parentTable` / the
    enclosing call's table). Arrays of plain primitives still become a **List** constant.
- **Table-name strings** parse to `ResolveTable(value)` calls — substituted with the **Select Table**
  utility (`grok.shell.tableByName(name)`, also broken platform-side), titled `table: <name>` so
  collapsed nodes stay readable.
- **Constant nodes** are titled after their value (`const: <value>`, via `constLabel` in
  utility-nodes.ts; empty → the type name). Because labels are user-editable, anything semantic must
  NOT read `node.label`: utility emission dispatches on `dgTypeName` (`utilityKind` in
  graph-compiler.ts), and the Value Output auto-typing checks `dgTypeName` too.
- **Outputs / variables**: there are **no Output nodes**. *Every* script variable
  (`BuiltGraph.outputVariables`, first-assignment order) feeds a real **`SetVar(variableName, value)`**
  func node (labeled `set: <name>`) — the single terminal per variable — wired from the variable's
  final ref. Running the flow registers each value in the context under its original name, so
  downstream consumers (a Select Table in a lower disjoint path, other scripts) can resolve it.
- **Inferred order edges**: a Select Table reads its table at runtime by name
  (`grok.shell.tableByName`) with no data edge back to the producer, so nothing *forces* the producer
  to run first. After all `SetVar`s are wired, `inferOrderEdges` matches each Select Table's
  `tableName` against the script's variable names — **normalized** (`normalizeName`: lowercase, drop
  non-alphanumerics) so the name↔friendlyName convention resolves (`Mol1KLocal` ↔ `"mol1K local"`) —
  and adds an **order edge** (exec-out → exec-in, see Execution-Ordering Edges) from that variable's
  `SetVar` to the Select Table. Creation scripts are linear/acyclic (a line references only
  already-created tables), so the edges never form a cycle. They are **excluded from the layout** (the
  `order` flag on `BuiltConnection`) — banding still relies on `orderedComponents`' name match — but
  they do constrain the topological sort, so reordering nodes can't break run-order.
- **Ordering**: after a *bare* call consumes a variable (a direct `GetVar`), the variable ref
  advances to that node's `<input>__pt` pass-through output. The next consumer connects there, so the
  topological sort reproduces the script's sequential line order (critical for in-place mutators like
  `addChemPropertiesColumns`), while the compiler still resolves the pass-through to the same
  expression — no spurious variables in the generated script.
- **Layout**: all imported nodes start **collapsed** (title bar only — expand per node as needed).
  The arrangement itself lives in [rete/graph-layout.ts](src/rete/graph-layout.ts) (`layoutGraph`), shared
  with the **Clean Layout** ribbon action so both produce the same result. Layered left-to-right with
  **one horizontal band per disjoint path** (weakly connected component):
  - Columns are **global** — shared x per layer, width = widest estimated node in that layer — so every
    edge points right and same-depth nodes line up across bands. The importer feeds its incrementally
    assigned layer map (`max(source layers)+1`); `FlowEditor.autoLayout` derives layers from the
    connection structure via `computeLayers` (longest-path).
  - Each component is a contiguous band; within a band/column, nodes order by predecessor barycenter
    and greedily stack (`max(nextFreeY, barycenter-h/2, bandTop)`) so chains read as straight lanes,
    branches fan out, and nothing overlaps.
  - Bands are stacked in **dependency order** (`orderedComponents`): a path that produces a table
    (its `SetVar` variable) is placed **above** the path that reads it through a `Select Table` node
    (matched by normalized name), with ties broken by node order. Because the execution topological
    sort ranks components by topmost-node `y`, this band order *is* the execution order — producers run
    before the consumers that read their tables.
  - `estimateNodeWidth`/`estimateNodeHeight` live in `graph-layout.ts` (re-exported from the importer
    for the layout-invariant tests: edges-point-right, no-overlap, producer-above-consumer).

The core is a **pure, synchronous, DOM-free** `buildCreationScriptGraph(script): BuiltGraph` — it
constructs `FlowNode` instances + connection records but touches no editor, so it is the unit-test
entry point. `applyGraphToEditor(graph, flow)` pushes it into a live editor;
`buildFlowFromCreationScript(flow, script)` does both.

UI: `File > Import Creation Script...` in the ribbon menu opens a dialog with a script textarea and a
"From table" picker prefilled from open tables that carry a creation script.
`FuncFlowView.loadFromCreationScript(script)` clears the canvas, builds, and zooms to fit.

## Creation-Script Emit ([compiler/creation-script-emitter.ts](src/compiler/creation-script-emitter.ts))

The **inverse** of the importer: `emitCreationScript(flow): {script, warnings}` compiles the graph back
into a grok-language creation script. Unlike [script-emitter.ts](src/compiler/script-emitter.ts) (pure,
DOM-free JS), it uses the platform's **own** serializer — `DG.Func.prepare(params).toString()` (Dart
`toConsole`/`valueToString`) — as the single source of truth, so it needs a **live backend** (tests run
on a local stand). The serialization contract (verified against the Dart source):

- `SetVar.prepare({variableName, value}).toString()` → `Name = <value>`; a `FuncCall` value recurses
  inline. `GetVar.prepare({variableName}).toString()` → the **bare** identifier — the *only* way to render
  a dataframe argument as a variable reference (a plain string or a real `DataFrame` both render as a
  *quoted name*).
- `column` arg → quoted column **name** (not JS `table.col(...)`); `column_list` → `[...]`; scalars/bools
  handled by the serializer; **optional params equal to their default are omitted** (so the importer's
  seeded `''` must be dropped — see `isEmptyLiteral` — or `sheetName = ""` leaks).

Graph → script mapping (walks `topologicalSort`, exec/order ports filtered by `isExecKey`):

- **Producer** = a func whose *real* (non-`__pt`) output is **consumed or anchored** → `X = f(...)`
  (assignment), downstream uses a bare `X`. The test is whether the real output is *used*, NOT whether the
  func declares one: `AddNewColumn` declares a real output but threads its table through a pass-through, so
  it stays a **bare call**.
- **In-place mutator / side-effecting call** = no consumed real output → bare `f(...)`; its `__pt` outputs
  forward the same variable (`forwardPassthroughs`).
- **SetVar node** = the variable-name anchor; `computeAnchors` walks its `value` back through the `__pt`
  chain (`walkToProducer`) to name the producer at the head, so an imported chain re-emits as
  `Mol1K = OpenFile(...)` / bare mutators / (no redundant `Mol1K = Mol1K`).
- **Select Table** → table-name string; **Select Column(s)** → column name(s); **Constants** → literals;
  **Input** node → bare `GetVar(paramName)` ref; **order edges** → ordering only, no line.
- **Warn & skip**: JS-only nodes (comparisons, ToString, FromJSON, ToJSON, Log/Info/Warning, Add Table
  View) have no grok-language form → a warning is collected and the node is skipped; consumers resolve to
  `skip`, so the rest still emits.

Round-trips faithfully: import → emit reproduces the chem/join example scripts (incl. qualified column
names, namespaced funcs, bare refs), and `emit → import → emit` is idempotent.

UI: the **Compile to Creation Script** ribbon action (`stream` icon, and `Script > Compile to Creation
Script...`) opens a dialog with the script, a copy button, and a warnings strip (`compileToCreationScript`
in [funcflow-view.ts](src/funcflow-view.ts)).

### Per-table split + Save (`creationScriptEditor`)

When the view edits *specific tables'* creation scripts — the `openCreationScriptFlowDialog`
(`meta.role: creationScriptEditor`) entry point passes `tableIds`, loaded into `DG.TableInfo[]` and
handed to `new FuncFlowView(tableInfos)` — the combined cascade must be split back into **one creation
script per table** and written to each.

- `emitCreationScriptsForTables(flow, tableVarNames)` (in [creation-script-emitter.ts](src/compiler/creation-script-emitter.ts))
  walks the graph **once** (shared with the single-script `emit()`), tagging every emitted line with its
  **owner** — the table variable it builds/mutates: a producer's anchor/assignment name, or, for a bare
  mutator, the variable resolved on its **first dataframe input** (`primaryTableOwner`). Lines are then
  bucketed by owner; a line whose owner isn't a requested table goes to `unassigned`. For an imported
  flow every line's owner is a real table variable, so the split is a clean partition.
- **Owner ↔ table** match key is the table's `.VariableName` tag (`tableInfo.tags[DG.Tags.VariableName]`,
  fallback `tableInfo.name`) — which equals the importer's `SetVar` variable names.
- Each line is stamped with a `//{"timestamp": <ms>}` comment from a **single global counter**
  (`Date.now()`, +1 ms per line) advanced in topological order across *all* tables — never per-table —
  so the value also encodes which table Datagrok builds first (an earlier-built table's lines carry
  smaller timestamps than a table that depends on it).
- UI: a primary **Save** ribbon button (`ui.bigButton`, only present when `tableInfos` is non-empty)
  opens `saveCreationScriptsDialog` — a horizontal `ui.tabControl`, one tab per table showing its
  standalone script, and a Save (OK) action that calls `TableInfo.saveCreationScript(script)` for each.
  An empty per-table script is legitimate (a table updated locally) and is shown plainly, not warned;
  only genuine emitter warnings (JS-only nodes, serialization failures) surface in the dialog.

## Tests ([src/tests/](src/tests))

`package-test.ts` imports the suites; run with `grok test --host localhost`. `test-utils.ts` provides
`makeEditor`/`destroyEditor` (a detached, off-screen `FlowEditor` whose data layer is populated
synchronously) and `BuiltGraph` query helpers (`nodesByFunc`, `sourceOf`, …).

| File | Category | Covers |
|---|---|---|
| `type-map-tests.ts` | Flow: type-map | `areTypesCompatible` matrix, `dgTypeToSlotType`, colors |
| `node-factory-tests.ts` | Flow: node-factory | `createNode`, registry, `ensureFuncNodeType` idempotency, pass-throughs |
| `compiler-tests.ts` | Flow: topological sort / script emitter / validator | order, cycles, emitted headers + body, instrumented mode, validation rules |
| `serializer-tests.ts` | Flow: serializer | serialize shape + round-trip topology, unknown-type skip |
| `minimap-tests.ts` | Flow: minimap | node rects + viewport drawn, `setMinimapCollapsed`, header-click collapse, hidden on empty canvas / shown once a node exists |
| `order-edge-tests.ts` | Flow: order edges | type isolation, exec ports on every node, order overrides `y` in the sort, sequenced-but-data-free emission, cycle detection, serialization round-trip |
| `layout-tests.ts` | Flow: layout | `computeLayers` (chain/diamond longest-path), `FlowEditor.autoLayout` (edges-point-right, no-overlap, producer-above-consumer in the editor) |
| `panel-tests.ts` | Flow: property panel | `stringChoiceOptions` (choices/nullable/current-preservation) + `propertyChoices` reading live func-input choices |
| `creation-script-import-tests.ts` | Flow: creation script import | exact `BuiltGraph` checks incl. the chem-properties example (column arg → Select Column wired to the table, pass-through ordering, output wiring), inferred order edges (friendly-name match, no-match, live-editor sort) + editor integration (emits `table.col(...)`, no `ResolveColumn`) |
| `function-browser-tests.ts` | Flow: function browser | exclusion list (no dev/test pkgs, `test*`, funccall wrappers), `categorizeFunc` placement (JoinTables→Combine, OpenFile→Data Sources, …), category order, `statusLabel` |
| `inspect-tests.ts` | Flow: inspect / slice | `sliceUpTo` (target + ancestors, excludes downstream/unrelated), `emitScript` `onlyNodeIds` filtering, `missingRequiredInputs` |
| `creation-script-emit-tests.ts` | Flow: creation script emit | round-trips (producer assignment, bare-call mutators in order, bare variable refs, join-by-name, friendly-name ref, no leaked optionals, full chem), `emit→import→emit` idempotency, warn-and-skip for JS-only nodes (needs a live backend) |

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
| Select Table | `let v = grok.shell.tableByName('name') ?? grok.shell.getVar('name') ?? …` (tries the exact, no-spaces, and lower-camel name variants) |
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

## Execution-Ordering ("order") Edges

A second, **data-free** connection type that expresses pure run-order — "node1 (and everything before it) must run before node2" — so users don't have to rely on vertical node position. Where pass-throughs are func-specific and tied to a data input, order edges work between **any** two nodes.

- **Ports**: `createNode` appends one **exec-in** (`__exec_in`, accepts many predecessors via `multipleConnections`) and one **exec-out** (`__exec_out`) of socket type `order` to *every* node — defined in [scheme.ts](src/rete/scheme.ts) (`EXEC_IN_KEY` / `EXEC_OUT_KEY` / `ORDER_SOCKET_TYPE` / `isExecKey`). They render as gray squares at the node's **top corners** (KNIME flow-variable style); the data-port rows exclude them.
- **Type isolation**: `areTypesCompatible` returns `false` whenever either side is `order` (checked *before* the `dynamic`/`object` wildcards), so a data port can never connect to an exec port and vice-versa.
- **Topological sort**: zero changes — the sort already works at the node level over `getConnections()`, so an order edge is just another dependency. Vertical position (`y`) degrades to a tiebreaker between genuinely-unordered nodes.
- **Compiler**: ignores exec ports — every port-iterating loop in `graph-compiler.ts` filters `isExecKey`, so order edges produce no variable, no data, and never leak `__exec*` into the output.
- **Layout**: **Clean Layout** (`FlowEditor.autoLayout`) **includes** order edges — it recomputes layers from scratch with `computeLayers`, so an order edge is just another forward dependency: the "after" node lands in a higher layer (further right) and explicit run-order shapes the layout left-to-right. The **importer** is the one exception: it reuses a *stale incremental* layer map (a `Select Table` was assigned layer 0 at creation, before its producer's `SetVar` existed), so feeding order edges to that map would draw a backward wire — the importer therefore **excludes** them (the `order` flag on `BuiltConnection`) and keeps its `orderedComponents` producer-above-consumer banding, with the order edge drawn as a diagonal between bands. (Click Clean Layout after import to re-flow it left-to-right.)
- **Visual**: static gray dotted edge (no data-flow march); `FlowEditor.tagConnectionElement` stamps `data-order="true"`, CSS does the rest.
- **Serialization**: free — order edges use the existing connection schema (`sourceOutput: __exec_out`, `targetInput: __exec_in`).

## Function Filtering ([rete/node-factory.ts](src/rete/node-factory.ts))

Functions with no inputs *and* no outputs are skipped. Functions whose role appears in `EXCLUDED_ROLES` (or any tag in `EXCLUDED_TAGS`) are skipped.

`registerBuiltinNodes()` populates the `FACTORIES` map with all built-in types. `registerAllFunctions()` discovers DG functions via `DG.Func.find({})` and registers a per-func factory under name `DG Functions/<role>/<funcName>` (or `DG Functions/<role>/<pkg>:<funcName>` on collision).

`createNode(typeName)` looks up the factory and stamps `dgTypeName` on the new instance — this is what the serializer persists.

## Type System (`types/type-map.ts`)

- `DG_TYPE_MAP`: DG type string → `{slotType, color}`. The slot color is what the React Socket component fills the dot with.
- `FUNC_NAME_COLORS`: per-function title-bar color, keyed by simple function name (case-insensitive). `getNodeColors(role, funcName)` checks this **before** role coloring, so specific functions can be pinned regardless of role (e.g. `SetVar` → red `#EF5350`, `GetVar` → light red). Add an entry to pin any function.
- `ROLE_COLORS`: DG role → title-bar color (white body always).
- `areTypesCompatible(out, in)`: source-of-truth for connection validity. Used by `TypedSocket.isCompatibleWith`. Permissive for `dynamic` and `object`; explicit pairs for `int↔double↔num` and `list↔string_list`.

`TypedSocket` ([rete/sockets.ts](src/rete/sockets.ts)) is one-instance-per-DG-type (cached), so reference equality holds. Its `isCompatibleWith` method is consulted at connection-pick time by `ClassicFlow.canMakeConnection`, which rejects incompatible drops before they enter the editor's data layer.

## Script Generation

Pipeline ([compiler/](src/compiler)):

1. **Topological sort** — Kahn's algorithm over `editor.getConnections()`, made deterministic and
   layout-aware: **disjoint subgraphs** (weakly connected components) run one after another, ranked by
   their topmost node's `(y, x)` — a path placed above another finishes completely before the lower
   one starts (lower paths may implicitly read what upper ones produced, e.g. a Select Table reading
   a table an upper path opened). **Within** a component, ready nodes are picked top-to-bottom
   (`y`, then `x`, then insertion order). `y` is only a *fallback* ordering — an explicit
   **order edge** (see Execution-Ordering Edges) makes "A before B" a real graph edge, so users aren't
   forced to encode run-order through vertical position.
2. **Compile** — every node becomes a `CompiledStep` with `inputs: Map<key, expr>`, `outputs: Map<key, varName>`, `properties`, `inputValues`. Variable names: camelCase of node label + first real output; collisions deduplicated by suffix. Func input resolution runs in two passes: ordinary inputs (connections, primitive `inputValues`) first, then unconnected **`column`/`column_list`** `inputValues` — these inline to `table.col('name')` / `[table.col(…), …]`, where the table comes from the node's `properties['columnTables']` association (so column args need no Select Column node; `tableExprForColumnParam`/`columnSelectionExpr`).
3. **Emit** — steps become JS lines; dataframe inputs first; `//input:` / `//output:` headers from properties + qualifiers.

This pipeline targets **JavaScript**. For the alternate **creation-script** target (grok-language cascade, via the platform's own serializer), see **Creation-Script Emit** above — it walks the same topological order but emits `prepare().toString()` lines instead of JS.

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

**SetVar preview**: `SetVar` declares no output, but the instrumented wrapper captures its incoming `value` as a synthetic output keyed by the variable name (`'<varName>': __ff_summarize(<valueExpr>)`), so clicking a SetVar node opens the docked output panel and renders the stored value by type (table → grid, column → sample, …) — same as any output-bearing node.

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

## Guide system (tutorials + how-to) — [`src/guide/`](src/guide)

Interactive, in-app onboarding. A non-invasive floating help button (`ff-help-fab`, bottom-left of
the canvas), a ribbon icon, and the Start panel's **“Take a 2-minute tour”** button all launch
guides: **tutorials** (multi-step) and **how-to questions** (single-answer). Each step highlights a
**concrete** UI element (a browser item, a canvas node, a context-panel input row, a ribbon icon —
never "the whole canvas") and **waits for the user's real action** before advancing — this is why
the [Test IDs](#test-ids-data-testid) matter: targets are addressed by `data-testid`/`data-*`.

The instruction card is **our own popup** (not `ui.hints.addHint`, whose injected ✕ collided with our
Exit link and whose placement never flips to the side with room). `computePlacement` (pure,
unit-tested) picks the best side and clamps fully on-screen; the card re-anchors on a timer and the
target gets a pulsing `.ff-guide-target` outline. No `scrollIntoView` (it used to jerk the whole UI).

| File | Role |
|---|---|
| `guide-model.ts` | Types (`Guide`/`GuideStep`/`GuideContext`/`GuideHost`; a step may set `highlights` for multiple pulse targets and `skipIf` to skip a satisfied prerequisite) + DOM-light helpers: `poll`, `waitForClick`, `untilClick`, `untilNodeType`/`untilFuncNode`, `untilMore{Nodes,Connections}`, `untilMoreCollapsed`(delta), `untilSectionExpanded`, `untilValueContains`/`untilValueMatches`(space-insensitive)/`untilValueNonEmpty`, `untilNodeSelected`/`untilNodeSelectedOfFunc`, `untilNodeRightOf`, `untilExists`; predicates `hasFuncNode`/`hasNodeType`/`nodeCount`; target resolvers `byTid`/`bySel`/`byNodeFunc`/`byNodeType`/`byBrowserFunc`/`byParam`(+`paramFieldSelector`)/`socketOf`; `copyToClipboard`/`prefillSearch`; and the pure `computePlacement`. |
| `guide-runner.ts` | `GuideRunner` — runs one guide at a time: highlight → own positioned card → `Promise.race(action, exit)` → re-anchor timer + cleanup; single ✕ exits; Skip; completion card. |
| `guide-content.ts` | The 4 tutorials (flagship `load-data-add-column` first — the Start-panel tour) + how-to questions (`TUTORIALS`, `QUESTIONS`). |
| `guide-launcher.ts` | `createHelpButton` + `openGuideMenu` (DG.Menu with Tutorials / How do I…? groups). |

The view ([funcflow-view.ts](src/funcflow-view.ts)) implements `GuideHost` (`getFlow`,
`showFunctionBrowser`, `anchorEl = helpButton`), mounts the button, and wires the Start-panel tour to
`TUTORIALS[0]`. Conditions, placement, and a live playthrough are tested in
[tests/guide-tests.ts](src/tests/guide-tests.ts). **Adding a guide**: append a `Guide` to `TUTORIALS`
or `QUESTIONS`; reuse the `until*` helpers and target a concrete element via `byTid`/`byNodeFunc`/
`byBrowserFunc`/`byParam`. Concrete targets (params, demo file) should be verified empirically (e.g.
`grok test --host localhost`) before being hard-coded.

## Test IDs (`data-testid`)

Every Flow UI surface carries a stable **`data-testid`** (the attribute Playwright's
`getByTestId()` reads by default) so the canvas, panels, ribbon, and dynamic lists can be
addressed deterministically from UI tests. Helper: [`utils/test-ids.ts`](src/utils/test-ids.ts) —
`tid(...parts)` builds an `ff-`-namespaced, slugified value; `setTid(el, ...parts)` stamps a plain
DOM element; React components pass `tid(...)` as a `data-testid` prop.

Convention: `ff-<area>-<thing>[-<dynamic-slug>]`. Dynamic parts (a function name, a param name, a
socket key, a node type) are slugified (`tidSlug`: lowercase, non-alphanumerics → `-`).

**Semantic identity attributes** (alongside the generic test-id, for finding a *specific* element by
its raw, un-slugified identity — the test-id `ff-node` is the same on every node):
- **Canvas node**: `data-node-id` (UUID), `data-node-type` (input/output/utility/func),
  `data-node-type-name` (registered type, e.g. `DG Functions/Data Sources/OpenFile`),
  `data-func` (qualified function name), `data-node-label` (current title).
- **Toolbox item / suggestion item**: `data-node-type-name`, `data-func`, `data-package`.
- **Browser section**: `data-section` (raw title).
- **Context-panel input row**: `data-param` (the input/param name).

| Surface | Test ID(s) |
|---|---|
| View shell | `ff-view`, `ff-root`, `ff-canvas`, `ff-statusbar`(+`-nodes`/`-links`/`-validation`) |
| Ribbon icons | `ff-ribbon-<action>` (`run`/`debug`/`continue`/`stop`/`view-script`/`save`/`open`/`undo`/`redo`/`layout`/`zoom-in`/`zoom-out`/`zoom-fit`/`toggle-browser`/`save-creation-scripts`) |
| Start panel | `ff-start-overlay`, `ff-start-panel`, `ff-start-bg`, `ff-start-template-<file>`, `ff-start-blank`, `ff-start-open`, `ff-start-import` |
| Function browser | `ff-browser`, `ff-browser-search`, `ff-browser-groupby`, `ff-browser-tree`, `ff-browser-section-<title>`, `ff-browser-item-<typeName>` |
| Canvas node | `ff-node` (+ `data-node-id`, `data-node-type`), `ff-node-title`/`-title-text`/`-status`/`-caret`/`-statusline`/`-hint`/`-description`/`-body`, `ff-exec-in`/`ff-exec-out`, `ff-socket-input-<key>`/`ff-socket-output-<key>`, `ff-socket-<type>` |
| Connections | `ff-connection`, `ff-edge-count`, `ff-waypoint-<i>` |
| Property panel | `ff-property-panel`, `ff-property-content`, `ff-property-title-row`, `ff-property-type-badge`, `ff-prop-input-<label/param>` |
| Editor extras | `ff-minimap`(+`-header`/`-toggle`), `ff-guide-overlay`, `ff-suggest-popup`/`-search`/`-list`, `ff-suggest-item-<typeName>`, `ff-hover-docs`, `ff-port-preview`, `ff-output-panel`, `ff-value-inspector`/`ff-value-previews`/`ff-preview-block-<name>` |

When adding a UI element, give it a `data-testid` via the helper. Tests live in
[tests/test-ids-tests.ts](src/tests/test-ids-tests.ts) (convention lock).

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
  - Func nodes: **Function** (description, full name, role) + **Input Parameters** (per-input editor for primitives via `node.inputValues`; a `string` input that declares `.choices` renders a **combo** instead of a text field — with a leading empty option when the property is `nullable` — via `propertyChoices`/`stringChoiceOptions`; `column` → a column-name text field, `column_list` → a comma-separated field, laid out side by side (≈70%/30%, `createColumnRow`) with a table-choice combo when the func has ≥2 dataframe inputs writing `properties['columnTables']`; "connected only" label otherwise).
  - Input nodes: **Input Configuration** (paramName, description, defaultValue, nullable, caption, type/semType filters, choices, min/max, showSlider).
  - Output nodes: **Output Configuration** (paramName + outputType combo for ValueOutput).
  - Utility nodes: **Configuration** for non-underscore properties (bool/number/text auto-detected).
  - **Connections** pane (collapsed by default): Inputs / Pass-through / Outputs grouped, with connection status from `flow.isInputConnected()` / `flow.getConnections()`.
- Editor helpers (auto-resizing textarea, number, toggle, combo) all support optional tooltips.
- After each property change, calls `flow.updateNode(node.id)` to re-render visible state (label, etc.).

### Function Browser ([panel/function-browser.ts](src/panel/function-browser.ts))

- Search input + Group-by selector. **Default mode is `category` ("what it does")**; other modes are `role` / `tags` / `package`.
- **`categorizeFunc(func, role)`** buckets by **input/output signature** (validated against the live catalog — see [docs/func-catalog-snapshot.md](docs/func-catalog-snapshot.md)). `FUNC_CATEGORIES` is also the **tree order** (Data Sources first):
  - **Data Sources** — outputs a table, **0 table inputs** (OpenFile, DB queries). A join is *not* a data source.
  - **Combine Tables** — **≥2 table inputs** (JoinTables, LinkTables).
  - **Transform Tables** — 1 table in → table out, or table in → no output (Aggregate, Unpivot, FilterRows).
  - **Column Operations** — outputs column/column_list (AddNewColumn, descriptors).
  - **Compute Values** — outputs a scalar. **Visualize** — viewer/view/widget/graphics (or role `viewer`). **Other** — the rest.
- Built-in sections (**Inputs / Outputs / Constants / Comparisons / Utilities / Debug**) are **all collapsed by default** — building-blocks reached for deliberately, not the first scan. The DG function categories below them are what leads.
- DG functions follow, ordered by `FUNC_CATEGORIES` in category mode (alphabetical in other modes); each category section is collapsed until clicked.
- Item tooltips: built-ins use `<description>. Double-click to add`; DG functions show `func.description (packageName)`.
- **Persistence**: the group-by mode and which sections are expanded are saved to `localStorage` (`funcflow.browser.v1`) and restored on open — keys are `b:<title>` for built-ins, `<mode>:<group>` for function groups. Search forces all open without touching the saved state (`isExpanded`/`setExpanded`).

### Catalog exclusions ([node-factory.ts](src/rete/node-factory.ts))

`shouldExcludeFunc(func, role, tags, pkgName)` drops ~⅔ of the raw catalog (firehose → ~800 usable nodes). Excluded when **any**: package ∈ `EXCLUDED_PACKAGES` (`Dbtests, ApiTests, UiTests, DevTools, Tutorials, ApiSamples, UsageAnalysis`); name is/starts-with `test`; role (comma-split) ∈ `EXCLUDED_ROLES` (now also `editor, cellEditor, panel, widgets, tooltip`); tag ∈ `EXCLUDED_TAGS`; or it takes a **`funccall` input** (command/dialog wrapper). Locked in by [tests/function-browser-tests.ts](src/tests/function-browser-tests.ts).

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

- **Collapse / expand**: click the **caret** (`▾`/`▸`) at the right of the title bar. The status dot is now **display-only** (U7 — so clicking the run indicator never folds the node out from under you). Sets `node.collapsed`; the React component re-renders. Collapsed nodes still expose their socket DOM in a hidden absolute-positioned row at the title bar's left/right edges so existing connections keep their endpoints.
- **Plain-language node status** (U5): below the title, completed/running/errored/stale nodes show a short line — `node.statusText`, set by `ExecutionVisualizer.statusLabel(status, detail)` (e.g. *Running… / Done · 1,204 × 8 / Error / Out of date*). The row/col detail comes from the captured `ValueSummary` via `summarizeOutputs` in the controller. Shown collapsed too. Idle → empty.
- **"Needs input" hint** (U6/U5b): when a node is idle/stale and a *required* input (a non-optional `dataframe`/`column`/`column_list`, listed in `FuncNode.requiredInputs`) is neither connected nor filled, an amber hint line ("Requires: table, molecules") replaces the status line and the status dot turns amber (`data-attention="true"`). Computed at render time by `missingRequiredInputs(node, isConnected)` ([scheme.ts](src/rete/scheme.ts)); `FuncFlowView.refreshNodeHints()` re-renders all nodes on `onGraphChanged` so it tracks wiring live.
- **Row counts on wires** (P0.2): after a run, each data connection is labelled at its midpoint with the count flowing through it (`_count` stuffed into the payload, rendered as `<text class="ff-edge-count">` by `FlowConnectionComponent`; set by `ExecutionController.labelOutgoingConnections`, cleared on edit/re-run via `FlowEditor.clearConnectionLabels`).
- **Inspect anywhere** (P1.1, the slice-compile keystone): right-click an output port → **"Run up to here & preview"** runs only the slice up to that node and opens its data. `sliceUpTo(flow, targetId)` ([graph-compiler.ts](src/compiler/graph-compiler.ts)) = the target + all transitive predecessors; `EmitOptions.onlyNodeIds` filters the emitted steps to that set (the set is closed under ancestors, so no surviving step references a dropped one); `ExecutionController.previewNodeData` runs it (skips global validation) and focuses the node on `run-complete` via `pendingPreviewNodeId`.
- **Context menu (right-click on a node)**: Collapse/Expand · Duplicate · Hide · Delete. Right-click on a connection: Delete connection. Built with `DG.Menu.popup().item(...).separator().show({causedBy: ev})` so the platform handles positioning, dismissal, and styling. Trigger comes from `area-plugin`'s `contextmenu` signal — `data.context` is `'root'` / a `FlowNode` / a `FlowConnection` and we branch on it.
- **Delete key (or Backspace)**: removes every selected node and all connections touching it. The handler is registered on `window.keydown` and skips events whose target is an `<input>`/`<textarea>`/`<select>` so typing in the property panel never deletes nodes.
- **Drag from toolbox**: `ui.makeDroppable` on the canvas container — accepts `DG.FileInfo` (creates an `OpenFile` node with `inputValues['fullPath']`) or `DG.Func` (looks up the registered factory by func name). Native HTML5 drag from the FunctionBrowser carries the typeName via the `FF_DRAG_MIME` data type; the canvas drop handler reads it and calls `addNodeAt` with the drop point.
- **Snap-to-grid + alignment guides**: `nodetranslate` is intercepted in `FlowEditor.wireEvents` — `computeSnap` returns either an alignment-snapped position (when an edge/center is within `alignThreshold` of another node's edge/center) or a 20px grid-snapped position. While dragging, dashed guide lines on a screen-space overlay (`.ff-guide-overlay`) mark the alignment axes; hidden on `nodedragged`.
- **Pan-to-new-node**: `addNodeAtCenter` calls `panToNode(id)` after one rAF, translating the AreaPlugin so the freshly-added node is at viewport center.
- **Overview minimap** (bottom-right): a screen-space SVG overlay built by `FlowEditor.installMinimap` (sibling of the guide overlay in `this.container`, *not* in the transformed canvas). Draws each node as a small rect (title color) plus the current viewport rectangle; click/drag pans (centers the viewport on the clicked graph point). Redraws are rAF-coalesced via `scheduleMinimapRedraw`, hooked into the area pipe (`translated`/`zoomed`/`render`/`nodetranslated`) and graph-change events. **Clicking anywhere on the header** minimizes it to the title bar (`data-collapsed`; the chevron is just an affordance). It's always present — there's no hide-entirely option. `FlowEditor.setMinimapCollapsed(bool)` sets the state programmatically; `FuncFlowView.setMinimapCollapsed(bool)` queues it until the (async) editor exists — used by `openCreationScriptFlowDialog`, which opens minimized in the dialog and expands on **Open In Editor**. The fit transform is stashed on the svg `data-*` so click→canvas mapping reuses it. rect-select and double-click-to-fit skip `.ff-minimap`.
- **Drag-out suggestion menu**: pointerdown on `.ff-socket-row-output .ff-socket` records the source slot; if pointerup lands on empty canvas (not a node, not a socket) AND no connection was created in between, [`openSuggestionMenu`](src/rete/flow-editor.ts) opens a searchable popup of every node type whose first input is type-compatible with the dragged source. **Ranked** (U9): Value Output first, then common next-step funcs (`COMMON_NEXT_FUNCS` — Join, AddNewColumn, Aggregate, Filter, …), then other built-ins, then DG funcs alphabetically. Selecting an item creates the node at the drop point and auto-connects the source to its first compatible input. Compatible-types lookup is cached per node type — see `findNodeTypesAcceptingInput` / `_sampleInputTypesCache` in `node-factory.ts`.
- **Start panel** (U1): an overlay over the empty canvas ([`buildStartPanel`](src/funcflow-view.ts)) — template cards (bundled `files/*.ffjson` via `_package.files.readAsText`), a Blank-canvas card, **Open a flow…** / **Import from a table…** buttons, and a discovery hint. `updateStartPanelVisibility()` shows it whenever `getNodeCount() === 0` (hooked into `onGraphChanged`, `initEditor`, `newFlow`); the ribbon **Flow › Templates…** re-opens it. Never face a blank page.

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
3. **Modifying script emission**: JS code generation lives in `script-emitter.ts`. The `EmitOptions.instrumented` flag toggles the try/catch + event path. New step kinds need handling in both clean (`emitFuncStep` / `emitUtilityStep`) and instrumented (`emitFuncStepInstrumented` / `wrapInstrumented`) paths. A new **utility** node that should also appear in **creation scripts** needs a branch in `emitUtility` of [creation-script-emitter.ts](src/compiler/creation-script-emitter.ts) (resolve its output to an `ArgValue`); otherwise it falls into warn-and-skip.
4. **Adding execution visual states**: extend `NodeExecStatus` in `execution-state.ts` and add a CSS rule `.ff-node-status[data-status="newstate"] { ... }`. Then update `ExecutionVisualizer` to set the status, no JS animation needed (CSS keyframes handle it).
5. **Adding output preview types**: extend the union in `OutputPreviewPanel` and add a case to `classifyOutput` / `buildPreview`. Make sure `getOutputTypeHints()` propagates the right declared type.
