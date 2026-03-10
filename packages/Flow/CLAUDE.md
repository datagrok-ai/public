# CLAUDE.md - Flow Package

**IMPORTANT: Always update this file after adding new features, nodes, or changing architecture.**

## Overview

Flow (FuncFlow) is an interactive visual function chain designer for Datagrok. It uses LiteGraph.js to let users compose Datagrok functions into executable JavaScript scripts via a node-based graph editor.

## Architecture

```
src/
├── package.ts                    # Entry: @funcflowApp, @fileViewer decorators
├── funcflow-view.ts              # Main view: 3-panel layout, ribbon, status bar
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
│   ├── validator.ts              # Pre-compilation validation
│   ├── topological-sort.ts       # Kahn's algorithm
│   └── graph-utils.ts            # Safe LGraph node accessor
├── panel/
│   ├── function-browser.ts       # Left sidebar: searchable function catalog
│   └── property-panel.ts         # Right sidebar: node properties editor
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
| String Input | string | nullable, choices, caption |
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
| Value Output | Configurable type selector |

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

Input annotation qualifiers are generated from node widget values:
- `{type: numerical}` from Column type filter
- `{semType: Molecule}` from Column semType filter
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

## Key Dependencies

- `litegraph.js` ^0.7.18 - Graph canvas library
- `datagrok-api` - Platform API (external, not bundled)

## Development Guidelines

1. **Adding new input nodes**: Add class in `input-nodes.ts`, register in `registerInputNodes()`, add to `inputNodes` array in `function-browser.ts`, handle in `buildInputLine()` in `script-emitter.ts` if special qualifiers needed
2. **Adding new utility nodes**: Add class in `utility-nodes.ts`, register in `registerUtilityNodes()`, add to `utilityNodes` array in `function-browser.ts`, add case in `emitUtilityStep()` in `script-emitter.ts`
3. **Adding new types**: Add to `DG_TYPE_MAP` in `type-map.ts`, add compatibility rules to `COMPATIBLE_TYPES`
4. **After any change**: Update this CLAUDE.md file
