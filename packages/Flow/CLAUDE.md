# CLAUDE.md - Flow Package

## Overview

Flow (FuncFlow) is an interactive visual function chain designer for Datagrok. It uses **Rete.js v2** to let users compose Datagrok functions, queries, and scripts into executable JavaScript scripts via a node-based graph editor.

The renderer is React-based (`rete-react-plugin`), but React is contained to the canvas: the rest of the package — view, panels, ribbon, status bar — is plain TypeScript on top of the Datagrok UI helpers, exactly as KetcherSketcher mounts Ketcher inside a `ui.div`.

## Rules — always follow

- **No page-global mutable state.** Several Flow views and editors can be alive in one page at
  the same time — file previews, Browse entity previews, the creation-script dialog, detached
  compile editors, multiple open tabs. Anything stored on `window`/`globalThis` is shared by all
  of them: a newer instance rebinding it, or a destroyed instance deleting it, silently breaks
  the others. (This exact bug once detached collapsed-node sockets and killed collapse carets —
  the old `window.__ff_editor` bridge.) Keep state on the owning instance and reach it through
  back-references (e.g. `FlowNode.editorBridge`). UI that belongs to a view must live **inside
  that view's root** (the bottom output panel is a splitter pane of the view — never a
  `grok.shell.dockManager` dock, which outlives and overlaps views).
  The one sanctioned global is `globalThis.__ffFlowLive`, the live-value registry the *emitted
  script* writes to (it runs outside any view). It is keyed by node id — unique per editor — and
  must only ever be cleared per-flow (`ExecutionController.clearLiveRegistry` deletes this flow's
  node ids, never the whole object).

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
│       ├── viewer-node.ts        # Manual viewer nodes (table.plot.fromType + setOptions)
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
│   └── output-preview.ts         # In-view bottom output panel (a pane of the view's splitter)
├── panel/
│   ├── function-browser.ts       # Left sidebar catalog
│   ├── property-panel.ts         # Side-panel node properties editor
│   └── column-picker.ts          # Column/columns selector menu seeded by the upstream table
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
- **Assignments** are `SetVar`-shaped calls carrying a string `variableName` + non-null `value`
  (`asAssignment` validates *those two inputs*, **not** the input count — the platform's `SetVarFunc`
  also has optional `outputName` / `outputIndex` for multi-output binding, so an assignment surfaces
  up to four inputs; gating on `=== 2` silently dropped every assignment);
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
- **SetVar node OR Output node** = the variable-name anchor (an Output's `paramName` IS its variable
  name — SetVar and Output are the same concept); `computeAnchors` walks the anchored value back through
  the `__pt` chain (`walkToProducer`) to name the producer at the head, so an imported chain re-emits as
  `Mol1K = OpenFile(...)` / bare mutators / (no redundant `Mol1K = Mol1K`), and a flow terminated by a
  Table Output emits `T = OpenFile(...)` with no intermediate variable.
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
| `type-map-tests.ts` | Flow: type-map | `areTypesCompatible` matrix, `dgTypeToSlotType`, colors, **`domainSection`/`domainCategory`/`isDomainOperation` chem/bio routing (operations only, not sources), disjoint package sets** |
| `node-factory-tests.ts` | Flow: node-factory | `createNode`, registry, `ensureFuncNodeType` idempotency, pass-throughs, suggestion-menu labels (friendly name + `funcCategory`, no "Uncategorized"), suggestion ranking (domain-in-play tier, exact-over-wildcard, used-func float), reverse suggestions (`findNodeTypesProducingOutput`: matching Input node leads, real-output-over-passthrough, domain boost) |
| `compiler-tests.ts` | Flow: topological sort / script emitter / validator | order, cycles, emitted headers + body, instrumented mode, validation rules, Select Table fail-fast guard, Output ⇄ SetVar same-contract emission |
| `serializer-tests.ts` | Flow: serializer | serialize shape + round-trip topology, unknown-type skip |
| `minimap-tests.ts` | Flow: minimap | node rects + viewport drawn, `setMinimapCollapsed`, header-click collapse, hidden on empty canvas / shown once a node exists |
| `order-edge-tests.ts` | Flow: order edges | type isolation, exec ports on every node, order overrides `y` in the sort, sequenced-but-data-free emission, cycle detection, serialization round-trip |
| `layout-tests.ts` | Flow: layout | `computeLayers` (chain/diamond longest-path), `FlowEditor.autoLayout` (edges-point-right, no-overlap, producer-above-consumer in the editor) |
| `panel-tests.ts` | Flow: property panel | `stringChoiceOptions` (choices/nullable/current-preservation) + `propertyChoices` reading live func-input choices |
| `creation-script-import-tests.ts` | Flow: creation script import | exact `BuiltGraph` checks incl. the chem-properties example (column arg → Select Column wired to the table, pass-through ordering, output wiring), inferred order edges (friendly-name match, no-match, live-editor sort) + editor integration (emits `table.col(...)`, no `ResolveColumn`) |
| `function-browser-tests.ts` | Flow: function browser | exclusion list (no dev/test pkgs, `test*`, funccall wrappers, **denylisted `nqName`s, machinery tags, `semantic_value`/filter-call, `meta.includeInFlow: false` opt-out; widgets KEPT**), `categorizeFunc` placement (JoinTables→Combine, OpenFile→Data Sources, **Chem/Bio operations→domain sections, chem/bio sources stay out, flow scripts→Workflows in every mode**), Cheminformatics section rendered, category order, `statusLabel`, queries grouped per-connection + kept out of the categories |
| `files-tree-tests.ts` | Flow: files tree | name-based test-ids on connection/folder/file rows; lazy expand loads + stamps a connection's files (Demo → demog.csv) |
| `execution-preview-tests.ts` | Flow: execution preview | widget/viewer outputs render their live `.root` and are renderable; context-panel meta names the kind (not `[object Object]`); a rootless widget is not renderable; panel state machine (hidden → expanded on first renderable output; minimize remembered — content updates never pop it up; `clear()` hides but keeps the preference; caret click toggles + fires `onStateChanged`, header body does not; disabled panels never show); same node + same state → no preview rebuild (new state / other node → rebuild); the panel is a pane of the view splitter, disabled in embedded views; **opening the panel minimizes the overview minimap** (one-shot, hidden → visible edge) |
| `viewer-tests.ts` | Flow: viewers | core viewer node types registered; a viewer node's table input / viewer output / type+specs; emits `plot.fromType` + `setOptions` (clean + instrumented); no table → no emission |
| `column-picker-tests.ts` | Flow: column picker | `column`/`column_list` inputs render a `ff-prop-pick-columns-<param>` icon; the request resolves the right dataframe input per column (JoinTables `keys1`→`table1`, `keys2`→`table2`); **viewer** axis options and **Select Column(s)** utilities get the picker too (resolving their `table` input); the request carries the icon as its menu `anchor`; `buildColumnMatchFilter` gates by `semType` / `columnTypeFilter`; no icon without an `onPickColumns` handler |
| `connect-interaction-tests.ts` | Flow: connect interaction | `soleCompatibleInput` (drop-on-node decision): one compatible free input → its key, several/zero → null, taken inputs excluded; `soleCompatibleOutput` (reverse drop): real output wins over passthroughs, sole-passthrough fallback, ambiguous/incompatible → null |
| `rerun-node-tests.ts` | Flow: rerun node | full run stashes each node's live outputs (`__ff_stash`, incl. passthrough table); a single-node re-run resolves external connections to `_ffLive(...)` (no `//input:`, no upstream call); a viewer re-plots the passthrough table; `canRerunNode` gating (required inputs + captured connection values; input nodes excluded) |
| `string-list-tests.ts` | Flow: string-list inputs, Flow: plain-list inputs | `string_list` inputs are seeded editable, render a text field, compile to a trimmed JS array, and omit when empty (`isStringListType`/`stringListToArrayLiteral` covered in `type-map-tests`); plain `list` inputs (incl. `list<string>` params) are seeded as arrays, render a DG List input via `forProperty`, compile to a JSON array literal, and omit when empty |
| `func-editor-tests.ts` | Flow: func editor | `shouldUseFunctionEditor` routing (allowlist / `editor:` meta / plain majority), `editorValueToPanelValue` conversions (column→name, lists, dataframe rejected), `tableParamForColumn` ladder, `applyEditorResult` (connected inputs win, column→name write-back), the full launcher ladder live (gate refuses unconnected → captured table opens the function's own dialog → close → write-back), **the dialog-hijack race** (a concurrent same-function call while the dialog is open — the autorun scenario — must execute uncanceled and must NOT resolve the round-trip; validated by disabling `ignoreEvent` → the concurrent call gets canceled), the header icon gating (needs a live backend) |
| `inspect-tests.ts` | Flow: inspect / slice | `sliceUpTo` (target + ancestors, excludes downstream/unrelated), `emitScript` `onlyNodeIds` filtering, `missingRequiredInputs` |
| `invalidation-tests.ts` | Flow: invalidation, Flow: autorun, Flow: in-place isolation | **In-place isolation**: instrumented emission snapshot-clones dataframe inputs (`__ff_clone`; the snapshot feeds the call, pass-through, and stash) while clean emission doesn't; a real Select Table → AddNewColumn run leaves the shell table + upstream capture pristine, previews show at-node state, and an autorun slice re-run is idempotent (validated by disabling `cloneDataframeInputs` → the shell table gets mutated). `sliceDownFrom` (forward closure), `applyGraphEdit` per-kind semantics (node-added → nothing; connection change → target+downstream stale, source + its live value kept; params-changed → node+downstream; node-removed → state forgotten), the editor's classified `onGraphEdited` stream (incl. `notifyNodeParamsChanged`, delete = connection-removed then node-removed, clear), **the user-side click guard** (opening the panel 3× per value-bearing node — Constants, Select Table/Column, Output, Int Input, OpenFile — emits zero `params-changed`; a real textarea edit reports once, same value again nothing; validated by disabling the guard → 3 spurious edits), `pendingNodes` (never-run/stale + downstream; fresh flow → empty), `expandToLiveBoundary` (captured boundary vs full ancestry), `AutorunScheduler` (debounce coalescing + dirty union, disabled/non-invalidating edits never run, `kick` on enable, busy → retry with kept set, skipped → wait for next edit, toggle-off cancels, re-entrant `hold`/`release` suspends firing and fires the backlog) |
| `creation-script-emit-tests.ts` | Flow: creation script emit | round-trips (producer assignment, bare-call mutators in order, bare variable refs, join-by-name, friendly-name ref, no leaked optionals, full chem), `emit→import→emit` idempotency, an Output node anchors the producer like a SetVar, warn-and-skip for JS-only nodes (needs a live backend) |

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

**Output ⇄ SetVar unification**: an Output node and a `SetVar` func node compile to the same
thing (`isSetVarNode` in [scheme.ts](src/rete/scheme.ts)). The emitted JS gives an Output node a
run-context registration (`grok.functions.call('SetVar', {variableName: <paramName>, …})` plus the
dataframe-runtime-name variant) and gives a SetVar node a script output (`//output: <type> <name>`
with the type inferred from the connected source socket — `setVarAsOutput` in
[script-emitter.ts](src/compiler/script-emitter.ts) — plus the `<name> = <value>;` assignment;
skipped when the variable name isn't a literal valid JS identifier). The creation-script emitter
anchors producers through Output nodes exactly as through SetVars (`computeAnchors`), and the
post-run auto-preview (`autoSelectFirstOutputNode`) treats SetVar terminals as outputs.

### Utility Nodes ([rete/nodes/utility-nodes.ts](src/rete/nodes/utility-nodes.ts))

| Node | Generated Code |
|------|----------------|
| Select Column | `let v = df.col('name')` |
| Select Columns | `let v = [df.col('a'), df.col('b')]` |
| Select Table | `let v = grok.shell.tableByName('name') ?? grok.shell.getVar('name') ?? …` (tries the exact, no-spaces, and lower-camel name variants), then **throws** `Select Table: no open table or variable named "…"` when all resolve to null — fail fast at the node instead of a cryptic downstream error |
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

`==`, `!=`, `>`, `>=`, `<`, `<=`, `Contains`, `Starts With`, `Ends With`, `Is Null`. **Hidden from the
toolbox for now** (still registered so saved flows load).

### Viewer Nodes ([rete/nodes/viewer-node.ts](src/rete/nodes/viewer-node.ts))

Manual nodes that build a Datagrok viewer from a wired table — **replacing the default viewer
*functions*** (which need a `TableView` lifecycle). A `ViewerNode` is `dgNodeType:'utility'` (so the
compiler's utility path resolves its `table` input) carrying `properties.viewerType` (a `DG.VIEWER`
value), `properties.viewerLook` (the accumulated options, look keys minus `#type`), and
`properties.viewerOptionSpecs` (the curated fields the panel renders). Output socket type is `viewer`.

- **Registration** ([node-factory.ts](src/rete/node-factory.ts)): `CORE_VIEWER_SPECS` register
  synchronously in `registerBuiltinNodes` (so saved viewer flows always deserialize); non-core package
  viewers (`DG.Func.find({meta:{role:'viewer'}})`, by `friendlyName ?? name`) register in
  `registerAllFunctions` via `registerDiscoveredViewers`. Both populate `VIEWER_NODE_TYPES` (the
  Viewers pane). Type name = `Viewers/<label>`.
- **Emit** ([script-emitter.ts](src/compiler/script-emitter.ts) `emitViewerStep`, dispatched in the
  utility branch on `properties.viewerType`): `let v = await <table>.plot.fromType('<Type>', {});` then
  `v.setOptions(<cleanViewerLook>)` (empty values dropped). Async `fromType` is **awaited**; the
  instrumented path summarizes `v` as a `'viewer'`. No table wired → emits nothing.
- **Live editing**: the preview's **“Edit settings”** button (added by `buildPreview` when an
  `onEditViewer` callback is threaded through `OutputPreviewPanel`) calls
  `FuncFlowView.editViewer(nodeId, viewer)` → `grok.shell.o = viewer` + a `DG.debounce(viewer
  .onPropertyValueChanged, 300)` subscription that writes `getOptions().look` (minus `#type`) back into
  `node.properties.viewerLook`. The property panel's `addViewerNodePane` edits the curated subset
  directly. (Empirically: `fromType`/`scatter` are async; passing a look to `fromType` throws —
  `setOptions` is the only reliable applier.)

### Function Nodes ([rete/nodes/func-node.ts](src/rete/nodes/func-node.ts))

`FuncNode` is built dynamically per `DG.Func`. Header uses `func.friendlyName` (split by `|`, last segment). Color comes from `ROLE_COLORS` map. Generates `await grok.functions.call('Pkg:funcName', {...})`. Primitive inputs seed `inputValues` with the **declared default** — `getParamDefault(prop)` = `defaultValue ?? initialValue`, read defensively, with one pair of wrapping quotes stripped (`unquoteDefault` — annotation defaults arrive double-encoded like `"'inner'"`), string-encoded bools/numbers coerced to their declared type (`'false'` must not compile to `true`) — else the type's zero value. In the panel, every DG input is initialized from the stored value via the guarded `stringValue` setter (`PropertyPanel.initInputValue`) because `ui.input.forProperty` / the `value` init option don't reliably load the editor.

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

`shouldExcludeFunc` drops a function when **any** holds: **`meta.includeInFlow: false`** declared on the function (checked first — the author opt-out; meta surfaces as `func.options.includeInFlow`, honoured as boolean `false` or string `'false'`; Flow's own `openCreationScriptFlowDialog` uses it); no inputs *and* no outputs; package ∈ `EXCLUDED_PACKAGES`; `nqName` ∈ the curated denylist [excluded-funcs.ts](src/rete/excluded-funcs.ts); name is/starts-with `test`; **role** ∈ `EXCLUDED_ROLES` (exact) **or a tag** ∈ `EXCLUDED_TAGS` (case-insensitive — sketchers/renderers/`Internal`/`Viewers`/… usually declare their kind as a *tag*, so the tag check is the biggest declutter lever); a `funccall` **or** `semantic_value` input; a `view`/`viewer` **or** `tablerowfiltercall`/`colfiltercall` output; or primitive-only signature. **`EXCLUDED_TAGS` deliberately omits `panel`/`widget`/`widgets`/`tooltip`** — widget-producing functions are usable in Flow (Widgets pane + preview); the `semantic_value` rule is what keeps right-click/context widgets out. The `EXCLUDED_FUNC_NQNAMES` denylist (incl. core plumbing: cache drops, project/publish, raw DB-query builders, UI-container builders) is empirically derived and **meant to be edited by hand** — see [docs/func-catalog-snapshot.md](docs/func-catalog-snapshot.md). This cut the catalog ~568 → ~283.

`registerBuiltinNodes()` populates the `FACTORIES` map with all built-in types. `registerAllFunctions()` discovers DG functions via `DG.Func.find({})` and registers a per-func factory under name `DG Functions/<role>/<funcName>` (or `DG Functions/<role>/<pkg>:<funcName>` on collision). It ends with `registerVariableFuncs()`: **SetVar / GetVar are force-registered** (via `ensureFuncNodeType`) even though the primitive-only rule excludes them from the catalog scan — every imported creation script terminates in SetVar nodes and Flow treats SetVar as an output, so a saved `.ffjson` containing them must deserialize without a prior import having registered them as a side effect.

`createNode(typeName)` looks up the factory and stamps `dgTypeName` on the new instance — this is what the serializer persists.

## Type System (`types/type-map.ts`)

- `DG_TYPE_MAP`: DG type string → `{slotType, color}`. The slot color is what the React Socket component fills the dot with.
- `FUNC_NAME_COLORS`: per-function title-bar color, keyed by simple function name (case-insensitive). `getNodeColors(role, funcName, category?)` checks this **first** (e.g. `SetVar` → red `#EF5350`, `GetVar` → light red). Add an entry to pin any function.
- `ROLE_COLORS`: DG role → title-bar color (white body always).
- `CATEGORY_COLORS`: task-category → title-bar color (incl. **Cheminformatics** pink / **Bioinformatics** deep-purple), the **fallback after role** so the role-less majority (JoinTables, AddNewColumn, chem properties, …) is colored by what it does instead of all gray. `FuncNode` computes the category as `domainCategory(pkg, inputTypes) ?? categorizeBySignature(...)` (same routing as the browser's `categorizeFunc`) and passes it to `getNodeColors`. Precedence: func-name → role → category → default gray.
- `CHEMINFORMATICS_PACKAGES` / `BIOINFORMATICS_PACKAGES` + `domainSection(pkg)` / `domainCategory(pkg, inputTypes)` / `isDomainOperation(inputTypes)`: the by-package domain routing (gated to operations that consume data), shared by node coloring and the toolbox grouping.
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
2. **Compile** — every node becomes a `CompiledStep` with `inputs: Map<key, expr>`, `outputs: Map<key, expr>` (the value expression a downstream connection reads, **not** always a bare variable — see below), `properties`, `inputValues`. Variable names: camelCase of node label + first real output; collisions deduplicated by suffix; a **leading digit is prefixed with `_`** (`uniqueVarName`) since `toCamelCase` keeps digits and `let 2xFlow…` is a syntax error. **Multi-output funcs**: `grok.functions.call` returns the value directly for a single-output func but an **object keyed by the declared output names** for a multi-output one — so `outputs`/`outputVarMap` store `<callVar>` for a lone output and `<callVar>.<outputName>` (`propertyAccessor`, bracket-quoted if not an identifier) for each of several. The old `<var>_<key>` scheme referenced variables that were never declared (2nd+ outputs came out `undefined`). The instrumented per-node output summary is keyed by the **slot key** (`"result1": __ff_summarize(<callVar>.result1)`), matching `labelOutgoingConnections`/`__ff_stash`. Func input resolution runs in two passes: ordinary inputs (connections, primitive `inputValues`) first, then unconnected **`column`/`column_list`** `inputValues` — these inline to `table.col('name')` / `[table.col(…), …]`, where the table comes from the node's `properties['columnTables']` association (so column args need no Select Column node; `tableExprForColumnParam`/`columnSelectionExpr`).
3. **Emit** — steps become JS lines; dataframe inputs first; `//input:` / `//output:` headers from properties + qualifiers.

This pipeline targets **JavaScript**. For the alternate **creation-script** target (grok-language cascade, via the platform's own serializer), see **Creation-Script Emit** above — it walks the same topological order but emits `prepare().toString()` lines instead of JS.

### Validation ([compiler/validator.ts](src/compiler/validator.ts))

- Empty graph → warning
- Cycles → error
- Column input without Table input → error
- Disconnected non-input nodes → warning
- Output node with no incoming connection → warning
- Empty / invalid JS-identifier `paramName` → error
- Duplicate variable name across input/output nodes **and SetVar `variableName`s** (one shared
  namespace — a SetVar doubles as an output) → error

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

**Dataframe input snapshots (in-place isolation)**: instrumented emission compiles with `CompileOptions.cloneDataframeInputs` — every **connected dataframe input of a func step** is snapshot-cloned before the call (`<snapVar> = __ff_clone(<srcExpr>)`, declared hoisted / assigned inside `try` since a `_ffLive` source can throw). The call args, the pass-through outputs, the `__ff_stash` entries, and inlined `table.col(...)` args (the compiler rewrites the input expression **before** column pass 2) all use the snapshot. Why: many platform functions transform tables **in place**; without the clone they'd mutate the very instance stashed as the *upstream* node's live value — upstream previews would show downstream columns, a shell table picked via Select Table would get modified, and an autorun slice re-run would apply the transform twice (non-idempotent). The mutated snapshot flows on through the passthrough, so downstream sees the transformed table; the upstream capture stays "the state at that node". Clean scripts don't clone (they run once from scratch; platform in-place idiom preserved). Viewers/utilities/SetVar don't clone (non-mutating; SetVar's `value` slot is `dynamic`, and its registered output should be the live final instance).

**In-place / threaded table capture**: when a func node has a dataframe input but **no real *dataframe* output**, the wrapper emits a synthetic entry `'<inputName> (modified)': __ff_summarize(<inputExpr>, 'dataframe')` so the post-execution table is previewable — one per **every** connected dataframe input (a two-table mutator with no output previews both). This covers both a pure in-place mutator (no outputs at all) AND a node whose real output isn't a table but still threads one through its passthrough — e.g. AddNewColumn returns a *column*, yet a viewer / the column picker wired to its `table →` passthrough needs that modified table (it's keyed off `!hasDataframeOutput`, not zero-outputs, so `cloneForNode` finds it).

**SetVar preview**: `SetVar` declares no output, but the instrumented wrapper captures its incoming `value` as a synthetic output keyed by the variable name (`'<varName>': __ff_summarize(<valueExpr>)`), so clicking a SetVar node opens the bottom output panel and renders the stored value by type (table → grid, column → sample, …) — same as any output-bearing node.

**Live-value registry (single-node re-run)**: every instrumented step also stashes its live outputs into a tab-global registry — `__ff_stash('<nodeId>', {<outputKey>: <value>, ...})` (real outputs → their variable, dataframe passthroughs → the threaded post-execution value), keyed by output socket key. **“Rerun this node only”** ([`ExecutionController.rerunNode`](src/execution/execution-controller.ts)) runs a *one-node* slice with `EmitOptions.liveExternalInputs`: `compileGraph(flow, liveBoundary)` resolves any connection whose source is outside the slice to a `_ffLive(nodeId, outputKey)` registry read (defined in the preamble) instead of an in-script variable — so nothing upstream re-executes. It's `preserveState` (other nodes untouched) and opens the node's preview on completion. Gated by `canRerunNode(nodeId)`: a func/viewer/utility node whose required inputs are all satisfied AND every connected input has a captured value (`hasLiveValue`). The registry is cleared on any non-preserve run (fresh full/slice run) and on structural graph change, so stale values never drive a re-run. Menu wiring: `FlowEditorCallbacks.onRerunNode` / `canRerunNode` → `showNodeContextMenu`.

## Execution Visualization

KNIME-inspired live feedback. The script runs in the same browser tab and communicates back via custom events.

### Architecture

- **ExecutionController** ([execution/execution-controller.ts](src/execution/execution-controller.ts)) — orchestrates runs. Subscribes to `funcflow.exec.<runId>`, updates `ExecutionState`, drives `ExecutionVisualizer`, pushes outputs into `OutputPreviewPanel`. Exposes callbacks `onBreakpointHit`, `onRunEnd`, `onNodeStateChanged`.
  - **Run readiness gate**: no run (autorun or manual) executes an **unready** node — one with a missing requirement per `nodeMissingRequirements` (see below) — nor anything downstream of it. `invalidNodes()` returns `{roots, cone}` where `cone = ⋃ sliceDownFrom(root)`; `runnableNodes()` = all nodes − `cone`. The full run (`executeFull`, behind `runInstrumented`/`debugInstrumented`) validates hard errors, then runs `runSet = allNodes − cone` (via `onlyNodeIds`), `warning`-ing which roots it skipped (or `error`-ing when nothing is ready). `runAutorun` prunes `cone` from its slice (or the would-be full set) and skips silently when the result is empty. Pruning the *whole downstream cone* keeps the emitted set closed under ancestors (no surviving step references a dropped one). Note a table-less viewer already emits nothing (`emitViewerStep` returns `[]`), but a table-less **Add Table View** / **Select Column** would emit `addTableView(undefined)` / `undefined.col(...)` and fault — the gate is what stops those.
- **ExecutionState** ([execution/execution-state.ts](src/execution/execution-state.ts)) — per-node status tracking (`idle` / `running` / `completed` / `errored` / `stale`) keyed by string node IDs.
- **ExecutionVisualizer** ([execution/execution-visualizer.ts](src/execution/execution-visualizer.ts)) — sets `node.dgStatus` and calls `flow.updateNode(id)` to re-render. The React Node component reads `dgStatus` and writes it to a `data-status` attribute. CSS does the rest (status circle color, pulse animation, body tint).
- **ValueInspector** ([execution/value-inspector.ts](src/execution/value-inspector.ts)) — runtime-value section in the side panel. DataFrame summaries embed a full `DG.Viewer.grid` preview. **"Add to workspace" (`addWorkspaceButton`, `data-testid="ff-add-to-workspace"`) is overlaid top-right on every rich preview** (dataframe, column, viewer), not only the context-panel meta row: dataframe/column → `grok.shell.addTableView(df.clone())`; viewer → `addViewerToWorkspace` opens a table view over a clone of the viewer's `dataFrame`, recreates the viewer on it (`plot.fromType(type, {})` + `setOptions(getOptions().look)` — a **new** viewer) and docks it `DG.DOCK_TYPE.RIGHT` of the view. On a viewer preview the button sits at `right:32px` (left of the gear) when the gear is present, else `right:6px`. A **column** output previews as a **one-column DataFrame grid**: `__ff_summarize` captures `clone: DG.DataFrame.fromColumns([col.clone()])` and `buildPreview` renders it as a grid (falling back to a small text sample if no clone was captured). **In-place column outputs preview the whole table, scrolled to the produced column**: `__ff_col_summary(col, inputTable, type)` (used by the func wrapper only when the step has **exactly one** connected dataframe input) tests, **by Dart-handle identity** (`.dart ===`, since each wrapper access is a fresh JS object), whether the output column is one of `inputTable`'s columns; if so it adds `tableClone` (a clone of the whole table) + `scrollToColumn` (the name), and `buildPreview` renders the table grid and `grid.scrollToCell(name, 0)` (deferred a frame — the grid isn't attached yet). When a node outputs a column, `buildValuePreviews` **suppresses** the threaded `"<input> (modified)"` passthrough dataframe preview — that table is still captured in the state (the column picker / inspect read it), just not shown, since the column preview (table-scrolled or one-column grid) is the meaningful result. **`buildValuePreviews` renders every renderable output, not just the first**: it collects each `buildPreview` block and, when there is more than one (a multi-output func's two dataframes, or a no-output mutator's several modified input tables), lays them out side by side in a `ui.splitH` with draggable dividers; a lone value fills the panel directly. A **widget/viewer** output keeps the live object (`{type:'widget'|'viewer', value}`, captured by reference during the in-tab run by `__ff_summarize`) and `buildPreview` mounts its `.root` directly; the property-panel meta names the kind (`widget`/`viewer`) instead of `[object Object]`.
- **OutputPreviewPanel** ([execution/output-preview.ts](src/execution/output-preview.ts)) — the bottom output panel, a **real pane of the view**: `FuncFlowView.initUI` mounts `panel.root` as the second item of a `ui.splitV([canvasBox, panel.root], …, true)`, so it is resizable via the splitter divider and can never linger over other views (no `grok.shell.dockManager` involvement — see Rules). The canvas pane is `flex: 1 1 0`; the panel pane is `flex: 0 0 auto`, making its explicit `height` the single source of truth (min/max clamps keep `splitV`'s own resize handling from drifting a minimized strip). Three states (`panelState`): **hidden** (start, and after `clear()` on graph change / new run), **expanded**, **minimized** (slim header strip; **only the caret** at the right edge of the header toggles — a fully clickable header would swallow near-miss clicks aimed at the splitter divider right above it). It is **not closable**: the first renderable output (clicking a completed node, or a run's focus node) expands it; once the user minimizes it the choice is **remembered** — later content updates in place and never pops the panel back up; only an explicit caret click restores it. `onStateChanged` lets the view show the divider only when expanded **and** minimize the overview minimap (`flow.setMinimapCollapsed(true)`) the first time the panel opens (hidden → visible edge only; the two share the bottom corner). Embedded hosts (the creation-script dialog) construct the view with `{outputPanel: false}` → the panel is disabled and never shows; `enableOutputPanel()` re-enables it (Open In Editor). **Re-render is identity-gated**: `showForNode` remembers the last `(nodeId, NodeExecState)` pair and skips rebuilding when both match and the panel is visible — `ExecutionState.setNodeStatus` always builds a fresh state object, so state reference identity IS value identity (re-clicking a node doesn't re-mount its grids; a re-run's new state does rebuild). The cache clears on `clear()`.

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

### Invalidation (classified, downstream-only)

The editor emits a **classified `GraphEdit`** (`onGraphEdited` callback in [flow-editor.ts](src/rete/flow-editor.ts): `node-added` / `node-removed` / `connection-added` / `connection-removed` / `params-changed` / `cleared`) for every result-affecting change; the coarse `onGraphChanged` stays as the "refresh UI" hook (status bar, hints, minimap) and also fires for cosmetic changes (annotations). `ExecutionController.applyGraphEdit(edit)` invalidates **only the downstream cone** (`sliceDownFrom` in [graph-compiler.ts](src/compiler/graph-compiler.ts) — forward mirror of `sliceUpTo`, follows data/pass-through/order edges alike) and returns the affected node-id set:

- `node-added` → nothing (not wired yet); `node-removed` → forget its state/visuals/live values (its connections' removal events already invalidated downstream).
- `connection-added/removed` → the **target** and downstream go stale; the source keeps its completed result (and its captured live value — still eligible for slice-boundary reuse).
- `params-changed` → the node and downstream. Reported by the **property panel** from every semantic editor helper — but only through a per-editor **`changeReporter`** guard (DG inputs fire `onValueChanged` on initialization, so an unguarded report would make every node click an "edit"; see the Property Panel section) → `PropertyPanel.paramsChanged()` → `FlowEditor.notifyNodeParamsChanged(nodeId)`. **Title/Description are cosmetic** and don't report. The function-editor dialog writeback reports from the view.
- Invalidation touches: `ExecutionState.markStale`, `ExecutionVisualizer.markStale` (node + incoming edges), outgoing-wire labels of invalidated nodes, their `__ffFlowLive` entries, and the output preview **only if it shows an invalidated node** (the node id is remembered so autorun can re-open it fresh).

There is no graph-version counter anymore — invalidation is entirely event-driven.

### Autorun ([execution/autorun.ts](src/execution/autorun.ts))

Ribbon **bolt icon** (`data-testid` ribbon/autorun; `.ff-autorun-toggle` = off: 0.8 opacity, default weight (outline glyph); `+ .ff-autorun-on` = on: blue + **font-weight 600**, which renders the FA bolt filled — the weight must NOT apply to the off state; dynamic tooltip) toggles `AutorunScheduler` (view-owned, off by default). **Turning it on immediately schedules `ExecutionController.pendingNodes()`** (every node without a completed result — never-run/stale/errored — plus downstream) via `AutorunScheduler.kick(dirty)`, so a fresh flow runs at once instead of waiting for the first edit. The scheduler accumulates the invalidated ids from `applyGraphEdit` and, `AUTORUN_DEBOUNCE_MS` (1 s) after the last edit, calls `ExecutionController.runAutorun(dirty, settings)`:

- `expandToLiveBoundary(flow, dirty, hasLive)` grows the dirty set upstream past nodes whose outputs are **not** captured; if the boundary is fully captured, only the slice runs (`onlyNodeIds` + `liveExternalInputs` + `preserveState` — same machinery as single-node rerun), else full run.
- Outcomes: `'started'` (dirty consumed) / `'busy'` (run in progress → retry next interval, set kept) / `'skipped'` (validation errors, or the run would prompt for script inputs — an input node inside the run set → wait for the next edit, set kept).
- Silent by design: no toasts, no dialog, no `autoSelectFirstOutputNode`; the only UI side effect is re-opening the output preview (content only, no selection change) if invalidation had closed it.
- `FuncFlowView.detach()` resets the scheduler so a pending debounce can't fire into a closed view.
- **`hold()` / `release()`** (re-entrant): suspend firing while a modal interaction is in progress — the view holds around the whole function-editor round-trip. Why this MUST happen: the editor dialog intercepts the global **`d4-before-run-action`** event, which fires for **every client funccall** (verified: `grok.functions.call(...)` and script `fc.call(...)` both fire it), and the interception matches **by func** — an autorun re-running the same function mid-dialog gets its call canceled and resolves the dialog round-trip early with the wrong funccall (stale values; the user's OK writes nothing — the "works the second time" race). Belt-and-braces beyond the hold: `FuncEditorLauncher.open` awaits `waitForRunIdle()` before opening, and passes `ignoreEvent: () => exec.state.isRunning` to `createFuncCallEditor` so events fired by an executing Flow run are never mistaken for the dialog's run action (the `onClose` fallback still resolves with the dialog `fc`, which holds the edited values).

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

### Bundled demos: `files/*.ffjson` vs `scripts/*.flow`

Two forms of the same two demos ship in the package:
- **`files/Workflow Demo.ffjson` / `files/Sequence demo.ffjson`** — raw ffjson docs, loaded by the Start panel's template cards (`FLOW_TEMPLATES` in [funcflow-view.ts](src/funcflow-view.ts), read via `_package.files.readAsText`).
- **`scripts/Workflow Demo.flow` / `scripts/Sequence demo.flow`** — the same graphs as **first-class flow scripts** (the platform auto-registers everything in `scripts/`; the `flow` scriptHandler declares `extensions: flow`). Each file is a `.flow` body: the annotation header (`//name` / `//language: flow` / `//tags: flow` / execution-ordered `//output:` lines) followed by the ffjson JSON — exactly what `flowScriptText` writes. They are **hand-committed** (generated from the ffjson), so they can drift; the **Flow: bundled flow scripts** test category regenerates the canonical body from each `.ffjson` via `flowScriptText` and asserts the header/output-ordering the committed files carry. If a future change alters that output (e.g. header ordering), the test fails — regenerate the `.flow` files to match.

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
| `guide-model.ts` | Types (`Guide`/`GuideStep`/`GuideContext`/`GuideHost`; a step may set `highlights` for multiple pulse targets and `skipIf` to skip a satisfied prerequisite) + DOM-light helpers: `poll`, `waitForClick`, `untilClick`, `untilNodeType`/`untilFuncNode`, `untilMore{Nodes,Connections}`, `untilMoreCollapsed`(delta), `untilSectionExpanded`, `untilValueContains`/`untilValueMatches`(space-insensitive)/`untilValueNonEmpty`, `untilNodeSelected`/`untilNodeSelectedOfFunc`, `untilNodeRightOf`, `untilExists`; **Files-tree** resolvers `byFileTreeConn`/`byFileTreeConnTri`/`byFileTreeFile` + gates `untilFileTreeConnExpanded`/`untilScrolledIntoView`/`untilFuncNodeWithInput`/`untilColumnCountAtLeast`(comma-separated pick count); predicates `hasFuncNode`/`hasNodeType`/`nodeCount`; target resolvers `byTid`/`bySel`/`byNodeFunc`/`byNodeFuncNth`(n-th identical node)/`byNodeType`/`byBrowserFunc`/`byParam`(+`paramFieldSelector`)/`socketOf`; **dialog-aware** `openDialogEl`(top `.d4-dialog`)/`preferDialog`(re-anchor the card to an open dialog so it isn't hidden behind it — the runner also drops the card below z-index 3000 while a dialog is open); `copyToClipboard`/`prefillSearch`; and the pure `computePlacement`. |
| `guide-runner.ts` | `GuideRunner` — runs one guide at a time: highlight → own positioned card → `Promise.race(action, exit)` → re-anchor timer + cleanup; single ✕ exits; Skip; completion card. |
| `guide-content.ts` | The 5 tutorials (flagship `load-data-add-column` first — the Start-panel tour; `interface-tour` walks every toolbox/ribbon/canvas/context-panel control, auto-adding a sample node as a prerequisite for the canvas/panel steps) + how-to questions (`TUTORIALS`, `QUESTIONS`). How-tos include `how-visualize` (demog → Scatter/Bar/Pie viewers, editing columns), `how-join` (two demog tables → Join Tables → output, filling keys via the column picker), and post-launch feature coverage pinned by the content test: `how-autorun` (bolt toggle, incremental reruns), `how-out-of-date` (downstream-only invalidation), `how-func-editor` (the “Open editor” dialog round-trip, gated on `.d4-flow-function-funccall-editor`), `how-reuse-flow` (saved flows → Workflows section), `how-rerun-node` (single-node rerun from captured values). The interface tour includes Autorun and the output panel; Save/Open steps describe the **platform** save (`.ffjson` import/export lives in the Flow menu). |
| `guide-launcher.ts` | `createHelpButton` + `openGuideMenu` (DG.Menu with Tutorials / How do I…? groups). |

The view ([funcflow-view.ts](src/funcflow-view.ts)) implements `GuideHost` (`getFlow`,
`showFunctionBrowser`, `anchorEl = helpButton`), mounts the button, and wires the Start-panel tour to
`TUTORIALS[0]`. Conditions, placement, and a live playthrough are tested in
[tests/guide-tests.ts](src/tests/guide-tests.ts). **Adding a guide**: append a `Guide` to `TUTORIALS`
or `QUESTIONS`; reuse the `until*` helpers and target a concrete element via `byTid`/`byNodeFunc`/
`byBrowserFunc`/`byParam`. Concrete targets (params, demo file) should be verified empirically (e.g.
`grok test --host localhost`) before being hard-coded.

## Auto-summaries (U12) — [`src/summary/`](src/summary)

Plain-language captions so a flow documents itself. `summarizeNode(node)` returns a short caption via,
in order: a built-in type summary (`BUILTIN_SUMMARIES`, keyed by registered type name) → a **curated**
function summary (`CURATED_FUNC_SUMMARIES`, keyed by bare lower-cased func name, reading the node's own
`inputValues`) → the function's **`description`** (same text as the context panel's Function pane) → a
deliberately-set **friendly name verbatim** (used whenever it differs from the raw name or contains a
space — never humanized, so acronyms like "InChI" survive) → a humanized identifier. The curated set
([`summary-defs.ts`](src/summary/summary-defs.ts), ~80 funcs) was chosen empirically from the live
catalog. `summarizeFlow(nodes, connections)` groups the graph into **disjoint pipelines** (union-find
over all connections), orders each left-to-right and pipelines top-to-bottom, and renders a numbered
list when there's more than one. Rendered on the node ([node-component.tsx](src/rete/node-component.tsx),
`.ff-node-summary`, only when no manual `description`) and via ribbon **Advanced → Describe this flow…**.
Pure/offline — tested in [tests/summary-tests.ts](src/tests/summary-tests.ts). **Add coverage**: drop a
key into `CURATED_FUNC_SUMMARIES`.

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
| Start panel | `ff-start-overlay`, `ff-start-panel`, `ff-start-bg`, `ff-start-template-<file>`, `ff-start-blank`, `ff-start-open`, `ff-start-first-flow` (flagship build tutorial), `ff-start-ui-tour` (hint link → interface tour) |
| Function browser | `ff-browser`, `ff-browser-search`, `ff-browser-groupby`, `ff-browser-tree`, `ff-browser-section-<title>`, `ff-browser-item-<typeName>` |
| Files / Queries panes | `ff-browser-files`, `ff-browser-queries`, `ff-browser-query-conn-<conn>` (+ `data-query-conn`); tree rows `ff-files-conn-<name>` / `ff-files-folder-<name>` / `ff-files-file-<name>` (+ `data-conn`/`data-folder`/`data-file`/`data-file-path`) |
| Canvas node | `ff-node` (+ `data-node-id`, `data-node-type`), `ff-node-title`/`-title-text`/`-status`/`-caret`/`-statusline`/`-hint`/`-description`/`-body`, `ff-exec-in`/`ff-exec-out`, `ff-socket-input-<key>`/`ff-socket-output-<key>`, `ff-socket-<type>` |
| Connections | `ff-connection`, `ff-edge-count`, `ff-waypoint-<i>` |
| Property panel | `ff-property-panel`, `ff-property-content`, `ff-property-title-row`, `ff-property-type-badge`, `ff-prop-input-<label/param>`, `ff-prop-pick-columns-<param>` (column picker icon on a column/column_list input) |
| Editor extras | `ff-minimap`(+`-header`/`-toggle`), `ff-guide-overlay`, `ff-suggest-popup`/`-search`/`-list`, `ff-suggest-item-<typeName>`, `ff-hover-docs`, `ff-port-preview`, `ff-output-panel`, `ff-value-inspector`/`ff-value-previews`/`ff-preview-block-<name>`, `ff-add-to-workspace`, `ff-viewer-edit` |

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
|  Output panel (splitter pane; hidden → minimized ↔ expanded) |
+--------------------------------------------------------------+
|  Status bar: Nodes / Links / Validation                      |
+--------------------------------------------------------------+
```

Property panel goes into Datagrok's native context panel via `grok.shell.o = propertyPanel.root`. The function browser is the view's toolbox; toggling `grok.shell.windows.showToolbox` shows/hides it (ribbon icon `list-ul`).

The bottom **Output panel** is a pane of the view's vertical splitter (not a dock) and is *lazy*: hidden until the first time the user clicks a node that has captured runtime values from a prior run ([`ExecutionController.showOutputsForNode`](src/execution/execution-controller.ts)); subsequent clicks update it in place, respecting a remembered minimized state. Emptied and hidden when the graph changes or a new run starts. See **OutputPreviewPanel** under Execution for the full state contract.

### Property Panel ([panel/property-panel.ts](src/panel/property-panel.ts))

- Title row at top (editable label) + node-type badge.
- Accordion with type-specific panes:
  - Func nodes: **Function** (description, full name, role) + **Input Parameters** (per-input editor via `node.inputValues`). **Primitives** (`string`/`int`/`double`/`num`/`qnum`/`datetime`/`bool`, in `PRIMITIVE_INPUT_TYPES`) are edited with **native Datagrok inputs** built from the parameter — `createPropertyInput` → `ui.input.forProperty(param, null, {value, tooltipText, onValueChanged})` — so the control honours the declared type, numeric range, **choices** (a choice-bearing string auto-renders as a combo — no `propertyChoices`/`stringChoiceOptions` wiring; those helpers are retained only for tests), and nullability. The row wrapper (`funcflow-dg-row`), not the DG input root, carries the `data-param`/test-id. `column` → a column-name **`ui.input.string`** (its own native caption), `column_list` → the same field for a comma-separated list, via `createColumnRow` → `createColumnFieldRow`. The picker icon and — for multi-table funcs (≥2 dataframe inputs, writing `properties['columnTables']`) — a table-chooser `<select>` (`funcflow-col-table-select`) are appended **inside** the input as trailing controls via `InputBase.addOptions(el)`; "connected only" label otherwise. (Columns use `ui.input.string`, **not** `forProperty` — a column value is just a name string; `forProperty` would build a live-table-bound column picker we can't back here.) The **picker icon** (`ff-prop-pick-columns-<param>`) — when `onPickColumns` is wired — opens a `ui.input.column(s)` dialog seeded by the upstream table (see Column Picker below). The picker is built by one shared `createColumnFieldRow({nodeId, label, isList, getValue, setValue, tableParam|tableSelect})` whose storage is fully delegated, so it serves **every** node with a column field — func inputs (`inputValues`), **viewer** look options (X/Y/Color/Size → `viewerLook`, table = the viewer's dataframe input), and the **Select Column / Select Columns** utilities (`columnName`/`columnNames` props). `dataframeInputKeys(node)` finds the table input(s) to pick from. Note the test-id slug is lowercased (`X`→`ff-prop-pick-columns-x`); the row's `data-param` preserves case. A `string_list` input renders a comma-separated **`ui.input.string`** (`createStringInput` — `forProperty` on a `string_list` builds the same Text input, so this is just explicit); the compiler runs it through `stringListToArrayLiteral` (trim, drop blanks → JS array) and **omits** an empty value so the function keeps its default. A plain **`list`** input — which is how **`list<string>` params surface at runtime** (`propertyType: 'list'`, `propertySubType: 'string'`; only core funcs like Aggregate report `string_list`) — goes through the **`forProperty` pathway** (`createPropertyInput`): DG builds its native **List** input whose `.value` is a JS **array**. `FuncNode` seeds `string_list` as `''` (`isStringListType`) and `list` as `[]`, so both are editable, not connection-only. The graph compiler emits a `list` value as a JSON array literal (empty → omitted); `creation-script-emitter.resolveInput` passes the array (empty → skip); the importer inlines an imported array-of-primitives `list` arg back into `inputValues` (round-trips) but lets arrays of calls/objects (e.g. column resolvers) fall through to generic resolution.
- **Column Picker** ([panel/column-picker.ts](src/panel/column-picker.ts)) — opens a real column / columns dialog so users pick from a list instead of typing names. Resolves the column's table input (`getTableParam()` — fixed for one dataframe input, the per-row select for multi-table funcs) → `flow.getInputSource()` → the upstream node. If unconnected: a hint to connect a table. If the upstream already ran: `ExecutionController.cloneForNode()` returns its captured `.clone`. If connected but not yet run: confirm, then `ExecutionController.produceTableForNode()` runs a **headless slice** (`sliceUpTo` + `onComplete`) and resolves with the produced `DataFrame`. The slice run uses `preserveState` (additive) — it does **not** `reset()` the run state / visuals, so it only recomputes its own slice and leaves other nodes' captured results intact; `cloneForNode` only reuses a `completed` (non-stale) result, so a second pick against the same table reuses its cached `DataFrame` instead of re-running, while a graph edit forces a recompute. Before the dialog opens, the resolved table goes through **`detectSemanticTypes`** (shared helper in [func-editor-launcher.ts](src/panel/func-editor-launcher.ts), also applied to every editor-seeded table) so semtype-filtered column inputs (Molecule, …) are populated — best-effort, a failure never blocks the dialog. Wired in `funcflow-view` as `propertyPanel.onPickColumns`.
- **Custom function editors** ([panel/func-editor-launcher.ts](src/panel/func-editor-launcher.ts) + [utils/func-editor-utils.ts](src/utils/func-editor-utils.ts)) — a function that ships its **own editor dialog** (`shouldUseFunctionEditor`: an `editor:` meta, surfacing as `func.options.editor`, or the explicit `EXPLICITLY_SUPPORTED_EDITABLE_FUNCTIONS` allowlist — e.g. `core:AddNewColumn`) gets a light-blue **"Open editor" button in the Input Parameters pane header** (`ff-prop-func-editor`, rendered by `decorateEditorHeader` only when the view wired `onEditFuncParams`). Clicking it runs `FuncEditorLauncher.open(node)`, which follows the **column picker's table ladder for every dataframe input at once**: all table inputs must be connected (balloon otherwise) → each resolves via `cloneForNode` (captured result) → anything missing gets one confirm dialog, then `produceTableForNode` slice runs. A `FuncCall` is then `func.prepare()`d from: the resolved tables, live captured values for other connected inputs (`ExecutionController.liveValue`, the `__ffFlowLive` registry), and the panel `inputValues` (column names resolved to real columns of their seeded table via `tableParamForColumn` — the same `columnTables`→suffix→first ladder as the compiler; `string_list` split to an array). `createFuncCallEditor(fc)` opens the platform dialog (`fc.edit()`, dialog polled via `pollDialogCreation`; the run action is intercepted by `d4-before-run-action` + `fc.status = 'Canceled'` so the values are saved without executing; the `d4-flow-function-funccall-editor` class locks the table selectors). On resolve, `applyEditorResult` writes values back into `node.inputValues`: **connected inputs are never overridden**, dataframe inputs skipped, only panel-editable slots touched, and columns come back as **name strings** (`editorValueToPanelValue`: column→name, column/string lists→comma-separated names, plain `list` stays an array, dataframes rejected). The view re-renders the panel on completion, holds the autorun scheduler for the whole round-trip, and reports the writeback via `notifyNodeParamsChanged`. **Never let a Flow run execute while the dialog is open**: the dialog intercepts the global `d4-before-run-action` (fired by EVERY client funccall) matching by func — see the Autorun section for the full race; `open()` waits for run-idle and passes `ignoreEvent` to `createFuncCallEditor`. Tested in [tests/func-editor-tests.ts](src/tests/func-editor-tests.ts).
  - Input nodes: **Input Configuration** (paramName, description, defaultValue, nullable, caption, type/semType filters, choices, min/max, showSlider).
  - Output nodes: **Output Configuration** (paramName + outputType combo for ValueOutput).
  - Utility nodes: **Configuration** for non-underscore properties (bool/number/text auto-detected).
  - **Connections** pane (collapsed by default): Inputs / Pass-through / Outputs grouped, with connection status from `flow.isInputConnected()` / `flow.getConnections()`.
- Editor helpers: primitives use native DG inputs (`createPropertyInput`); the bespoke helpers (auto-resizing textarea, number, toggle, combo) remain for the non-func panes (Title/Description, Input/Output config, Utility/Constants — which edit ad-hoc `properties`, not `DG.Property` objects) and all support optional tooltips.
- **Every semantic editor reports edits through a per-editor `changeReporter(initial)`** — NEVER call `paramsChanged()` directly from a handler. Creating/initializing a DG input (`ui.input.forProperty`, `initInputValue` setting `stringValue`) can fire `onValueChanged` immediately (not guaranteed for every input type, so "skip the first event" is wrong too); unguarded, merely clicking a node (which rebuilds the panel) counts as an edit — invalidating results and, with autorun on, rerunning the flow per click. The reporter keeps the last-seen value and reports only a REAL change (`sameValue`: scalars by string form — `5`≡`'5'`, `null`≡`''` — arrays/objects by JSON). Any new input wired into the panel must follow this pattern. Title/Description are cosmetic and never report. User-side regression test: `invalidation-tests.ts` opens the panel 3× per value-bearing node (incl. OpenFile with a fullPath) and asserts zero `params-changed`.
- After each property change, calls `flow.updateNode(node.id)` to re-render visible state (label, etc.).

### Function Browser ([panel/function-browser.ts](src/panel/function-browser.ts))

- **Files pane** (first, `ff-browser-files`, expanded by default): the KNIME-style file browser
  ([utils/files-browser-tree.ts](src/utils/files-browser-tree.ts), `getFilesBrowser`) listing every
  file connection + its folders/files. Built **once** and reused across renders (so expanded folders /
  scroll survive a search keystroke). **Double-click or drag a file → an `OpenFile` node** pre-set to
  that file (`onFileDoubleClick` callback → `addOpenFileNode`; drag uses the existing `DG.FileInfo`
  droppable). Every row carries a name-based `data-testid` (`ff-files-conn-<name>` /
  `ff-files-folder-<name>` / `ff-files-file-<name>`) plus raw `data-conn`/`data-folder`/`data-file`/
  `data-file-path` attributes.
- **Queries pane** (`ff-browser-queries`): every `DG.DataQuery` from the registry, **excluded from the
  category list** and grouped into one sub-accordion **per connection** (`queryConnectionName` =
  `connection.friendlyName ?? connection.name`; `ff-browser-query-conn-<name>` + `data-query-conn`),
  **regardless of the group-by mode**.
- **Viewers pane** (`ff-browser-viewers`): the manual **viewer nodes** (`VIEWER_NODE_TYPES`, core charts
  first, then discovered package viewers) — see [Viewer Nodes](#viewer-nodes-retenodesviewer-nodets).
- **Widgets pane** (`ff-browser-widgets`): functions whose output is a `widget` (`funcOutputsWidget`),
  also kept out of the categories.
- Search input (with a **clear ✕**, `ff-browser-search-clear`) + Group-by selector. **Default mode is `category` ("what it does")**; other modes are `role` / `tags` / `package`.
- **`categorizeFunc(func, role, packageName?)`** buckets by **domain then signature** (validated against the live catalog — see [docs/func-catalog-snapshot.md](docs/func-catalog-snapshot.md)). `FUNC_CATEGORIES` is also the **tree order** (Data Sources first):
  - **Domain wins (for operations only)**: a function from a `CHEMINFORMATICS_PACKAGES` / `BIOINFORMATICS_PACKAGES` package (in [type-map.ts](src/types/type-map.ts)) groups under **Cheminformatics** / **Bioinformatics** — via `domainCategory(pkg, inputTypes)`, which routes it there **only if it takes a dataframe/column input** (`isDomainOperation`); a chem/bio *source* (produces a table from scalars — a DB fetch, generator, or query) falls back to its signature category (Data Sources). It gets the matching title-bar color (pink / deep-purple). `DG.DataQuery` instances are filtered from the tree entirely into the **Queries pane**. Everything else falls through to the signature categories:
  - **Data Sources** — outputs a table, **0 table inputs** (OpenFile, DB queries). A join is *not* a data source.
  - **Combine Tables** — **≥2 table inputs** (JoinTables, LinkTables).
  - **Transform Tables** — 1 table in → table out, or table in → no output (Aggregate, Unpivot, FilterRows).
  - **Column Operations** — outputs column/column_list (AddNewColumn, descriptors).
  - **Compute Values** — outputs a scalar. **Visualize** — viewer/view/widget/graphics (or role `viewer`). **Other** — the rest.
  - **Workflows** — saved flows (`DG.Script` with language `flow`, `isWorkflowFunc` in
    [node-factory.ts](src/rete/node-factory.ts)); checked before everything else (their signature
    would misfile them under Data Sources) and forced in **every** group-by mode, not just
    category — a flow's role/tags/package say nothing useful.
- Built-in sections (**Inputs / Outputs / Constants / Utilities / Debug**) are **all collapsed by default** — building-blocks reached for deliberately, not the first scan. (The **Comparisons** group is currently hidden from the toolbox — its nodes stay registered so saved flows still load.) The DG function categories below them are what leads.
- DG functions follow, ordered by `FUNC_CATEGORIES` in category mode (alphabetical in other modes); each category section is collapsed until clicked. **Exception:** in category mode the **Cheminformatics** and **Bioinformatics** sections (`DOMAIN_CATEGORIES`) are rendered first, right after the Queries pane (before Viewers/Widgets/built-ins), then the remaining categories — so the domain science leads the toolbox.
- **Parameter descriptions**: `FuncNode` captures each slot's description (`getParamDescription` → `inputDescriptions`/`outputDescriptions`) and the source package (`dgPackageName`) from the live `DG.Func`. The node component shows them as socket-row `title` tooltips; the property panel's input rows use `buildFuncInputTooltip` (which reads the same description) and the **Function** pane shows a **Package** row (`ff-prop-func-package`).
- **Parameter captions (display only)**: `getParamDisplayName(prop)` resolves an input's **display label** — its declared `caption` (`options.caption`, else the `friendlyName`/Dart caption getter) when non-empty, else `prop.name`. Used for the **node input slot label** (`FuncNode` passes it as the `ClassicPreset.Input` label; the slot **key** stays `inp.name`) and the **context-panel** field labels: the "(connected only)"/"connected" info rows, and the `ui.input.string` caption for `column`/`column_list` (`createColumnFieldRow` `caption ?? label`) and `string_list` (`createStringInput`'s `caption` arg). `ui.input.forProperty` (primitives + `list`) already reads the caption itself. **Identity is never affected** — `inputValues`, connections, `data-param`/test-ids, and compiled args all stay keyed by `prop.name`.
- Item tooltips: built-ins use `<description>. Double-click to add`; DG functions show `func.description (packageName)`.
- **Persistence**: the group-by mode and which sections are expanded are saved to `localStorage` (`funcflow.browser.v1`) and restored on open — keys are `b:<title>` for built-ins, `<mode>:<group>` for function groups. Search forces all open without touching the saved state (`isExpanded`/`setExpanded`).

### Catalog exclusions ([node-factory.ts](src/rete/node-factory.ts))

`shouldExcludeFunc(func, role, tags, pkgName)` drops the raw catalog (firehose ~2,357) to ~285 usable nodes. Excluded when **any**: package ∈ `EXCLUDED_PACKAGES` (`Dbtests, ApiTests, UiTests, DevTools, Tutorials, ApiSamples, UsageAnalysis`); **`func.nqName` ∈ [`EXCLUDED_FUNC_NQNAMES`](src/rete/excluded-funcs.ts)** (the curated, hand-editable denylist — helpers, internal twins, demo/test, plumbing); name is/starts-with `test`; a **role OR tag** ∈ `EXCLUDED_ROLE_TAGS` (the combined lower-cased set: `EXCLUDED_ROLES` + `widget, viewers, internal, @editors, design, filehandler` — checked against **both** the role field and every tag, since panels/widgets/sketchers usually declare their kind as a *tag*); it takes a **`funccall`** or **`semantic_value`** input; it **outputs** a `view`/`viewer` or a `tablerowfiltercall`/`colfiltercall` (filter-DSL builder); or it is **primitive-only**. Locked in by [tests/function-browser-tests.ts](src/tests/function-browser-tests.ts). The nqName denylist + tag/structural rules were derived empirically (live catalog dump + per-package source assessment); see [docs/func-catalog-snapshot.md](docs/func-catalog-snapshot.md).

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
- **"Needs input" hint** (U6/U5b): when a node is idle/stale and a *requirement* is unmet, an amber hint line ("Requires: table, molecules") replaces the status line and the status dot turns amber (`data-attention="true"`). Requirements are two kinds: **required inputs** — a non-optional `dataframe`/`column`/`column_list` socket (`FuncNode.requiredInputs`; viewers list `['table']`) neither connected nor filled — and **required properties** (`FlowNode.requiredProps`: panel values that must be set, e.g. Select Column's `columnName`, Select Columns' `columnNames`, Select Table's `tableName`; Add Table View adds `requiredInputs:['table']`). Computed at render time by `nodeMissingRequirements(node, isConnected)` = `missingRequiredInputs` ⊕ `missingRequiredProps` ([scheme.ts](src/rete/scheme.ts)) — the **same** predicate the run-readiness gate uses (a hinted node is exactly one the run skips). `FuncFlowView.refreshNodeHints()` re-renders all nodes on `onGraphChanged` so it tracks wiring live. **Adding a required input/prop to a new node type automatically flags it in the hint and excludes it from runs.**
- **Row counts on wires** (P0.2): after a run, each data connection is labelled at its midpoint with the count flowing through it (`_count` stuffed into the payload, rendered as `<text class="ff-edge-count">` by `FlowConnectionComponent`; set by `ExecutionController.labelOutgoingConnections`, cleared on edit/re-run via `FlowEditor.clearConnectionLabels`).
- **Inspect anywhere** (P1.1, the slice-compile keystone): right-click an output port → **"Run up to here & preview"** runs only the slice up to that node and opens its data. `sliceUpTo(flow, targetId)` ([graph-compiler.ts](src/compiler/graph-compiler.ts)) = the target + all transitive predecessors; `EmitOptions.onlyNodeIds` filters the emitted steps to that set (the set is closed under ancestors, so no surviving step references a dropped one); `ExecutionController.previewNodeData` runs it (skips global validation) and focuses the node on `run-complete` via `pendingPreviewNodeId`.
- **Context menu (right-click on a node)**: Run up to here & preview (when `onPreviewNode` is wired — same slice-compile as the output-port menu, just more discoverable) · **Rerun this node only** (when `onRerunNode` is wired AND `canRerunNode(id)` — re-executes just this node from captured upstream values, see the live-value registry above) · Collapse/Expand · Duplicate · Delete. Right-click on a connection: Delete connection. Built with `DG.Menu.popup().item(...).separator().show({causedBy: ev})` so the platform handles positioning, dismissal, and styling. Trigger comes from `area-plugin`'s `contextmenu` signal — `data.context` is `'root'` / a `FlowNode` / a `FlowConnection` and we branch on it.
- **Delete key (or Backspace)**: removes every selected node and all connections touching it. The handler is registered on `window.keydown` and skips events whose target is an `<input>`/`<textarea>`/`<select>` so typing in the property panel never deletes nodes.
- **Drag from toolbox**: `ui.makeDroppable` on the canvas container — accepts `DG.FileInfo` (creates an `OpenFile` node with `inputValues['fullPath']`) or `DG.Func` (looks up the registered factory by func name). The node is placed **at the drop pointer** (`args.dropEvent` → `addNodeByTypeAtDrop` → `addNodeByTypeAt(typeName, clientX, clientY)`), matching the native-drag behavior; only double-click / programmatic adds fall back to center. Native HTML5 drag from the FunctionBrowser carries the typeName via the `FF_DRAG_MIME` data type; the canvas drop handler reads it and calls `addNodeAt` with the drop point.
- **Snap-to-grid + alignment guides**: `nodetranslate` is intercepted in `FlowEditor.wireEvents` — `computeSnap` returns either an alignment-snapped position (when an edge/center is within `alignThreshold` of another node's edge/center) or a 20px grid-snapped position. While dragging, dashed guide lines on a screen-space overlay (`.ff-guide-overlay`) mark the alignment axes; hidden on `nodedragged`.
- **Pan-to-new-node**: `addNodeAtCenter` calls `panToNode(id)` after one rAF, translating the AreaPlugin so the freshly-added node is at viewport center.
- **Overview minimap** (bottom-right): a screen-space SVG overlay built by `FlowEditor.installMinimap` (sibling of the guide overlay in `this.container`, *not* in the transformed canvas). Draws each node as a small rect (title color) plus the current viewport rectangle; click/drag pans (centers the viewport on the clicked graph point). Redraws are rAF-coalesced via `scheduleMinimapRedraw`, hooked into the area pipe (`translated`/`zoomed`/`render`/`nodetranslated`) and graph-change events. **Clicking anywhere on the header** minimizes it to the title bar (`data-collapsed`; the chevron is just an affordance). It's always present — there's no hide-entirely option. `FlowEditor.setMinimapCollapsed(bool)` sets the state programmatically; `FuncFlowView.setMinimapCollapsed(bool)` queues it until the (async) editor exists — used by `openCreationScriptFlowDialog`, which opens minimized in the dialog and expands on **Open In Editor**. The fit transform is stashed on the svg `data-*` so click→canvas mapping reuses it. rect-select and double-click-to-fit skip `.ff-minimap`.
- **Connect-mode highlighting** (`beginConnectHints`/`endConnectHints`): on `connectionpick`, the container gets `.ff-connecting` (dims every `.ff-node`), the drag origin gets `.ff-node-source`, and each opposite-side socket that's type-compatible (via `TypedSocket.isCompatibleWith`) gets `.ff-socket-compat` with its node `.ff-node-compat` (green glow). Symmetric: output picks light compatible inputs, input picks (existing-connection tail) light compatible outputs. Exec-order ports are excluded. Cleared on `connectiondrop` and on a safety `pointerup` (covers Esc-cancel).
- **Drop-on-node shortcut**: an output drag that misses every socket but lands on a node body connects to that node's **sole** compatible, unwired input (`soleCompatibleInput` → `addConnectionByKeys`) — no pixel-hunting for the dot, and it works on collapsed nodes whose pins aren't drawn. Zero/several candidates → no-op (aim at a pin); only empty-canvas drops fall through to the suggestion menu. The node under the drop is found via `document.elementsFromPoint`. **Reverse direction too**: an *input* drag (a fresh drag out of an input socket — tracked as `dragInSource` on `connectionpick` side `input`) dropped on a producer node's body connects from that node's one obvious output via `soleCompatibleOutput`: the **sole compatible real output wins**; only when no real output is compatible does the **sole compatible pass-through** qualify (e.g. a table drag onto AddNewColumn → its `table__pt`, since its real output is a column). Zero or several candidates in the winning group → abort; already-wired outputs stay eligible (outputs feed many consumers). Empty-canvas input drops open the **reverse suggestion menu** (`openReverseSuggestionMenu` → `findNodeTypesProducingOutput`): every node type with a compatible **output or pass-through** (probed via `getOutputTypesForType`, split real vs `__pt`), ranked with the same context heuristics but with the terminals swapped — tier 0 the **matching Input node** (Table Input for a table drag), tier 1 the science in play, tier 2 **Data Sources** funcs + `COMMON_NEXT_FUNCS`, then built-ins, then the rest; within a tier **real-output matches precede pass-through-only threaders** (`realOutput`), then exact-over-wildcard, used-func float, alphabetical. Picking one creates the node at the drop point and wires its first compatible output (real over pass-through) into the dragged input.
- **Drag-out suggestion menu**: pointerdown on `.ff-socket-row-output .ff-socket` records the source slot; if pointerup lands on empty canvas (not a node, not a socket) AND no connection was created in between, [`openSuggestionMenu`](src/rete/flow-editor.ts) opens a searchable popup of every node type whose first input is type-compatible with the dragged source. **Ranked by canvas context** (U9; `findNodeTypesAcceptingInput(sourceType, SuggestionContext)` — the context carries the drag-source's package, all canvas packages, and all canvas func names, built in `openSuggestionMenu`): tier 0 Value Output; tier 1 **the science in play** — Cheminformatics/Bioinformatics funcs (by `funcCategory`) when the drag source's package is of that domain (`domainSection`), falling back to any domain already on the canvas for domain-less sources (OpenFile, utilities); tier 2 common next-step funcs (`COMMON_NEXT_FUNCS` — Join, AddNewColumn, Aggregate, Filter, …); tier 3 other built-ins; tier 4 remaining DG funcs. Within a tier: **exact type match** (an input slot equals the dragged type, vs a `dynamic`/`object` wildcard — the `exact` flag) → **already-used-on-canvas** funcs → alphabetical. **Labels match the toolbox**: a DG func shows its **friendly name** with the **"what it does" category** in parentheses (`funcCategory(info)` — the cached domain→signature routing, so `Add New Column  (Column Operations)`, never the raw `funcName (Uncategorized)` role segment from the typeName); built-ins show their plain name. The popup search also matches the `typeName`, so typing the raw func name still finds it. Selecting an item creates the node at the drop point and auto-connects the source to its first compatible input. Compatible-types lookup is cached per node type — see `findNodeTypesAcceptingInput` / `_sampleInputTypesCache` in `node-factory.ts`.
- **Start panel** (U1): an overlay over the empty canvas ([`buildStartPanel`](src/funcflow-view.ts)) — template cards (bundled `files/*.ffjson` via `_package.files.readAsText`), a Blank-canvas card, a **Create your first flow** button (runs the flagship `load-data-add-column` tutorial), an **Open a flow…** button, and a discovery hint whose inline link launches the **interface tour**. (Importing a table's creation script is still available from the ribbon's Flow menu.) `updateStartPanelVisibility()` shows it whenever `getNodeCount() === 0` (hooked into `onGraphChanged`, `initEditor`, `newFlow`); the ribbon **Flow › Templates…** re-opens it. Never face a blank page.

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
