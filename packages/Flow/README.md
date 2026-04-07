# Flow

Flow is a visual function chain designer for the [Datagrok](https://datagrok.ai) platform. It provides a node-based graph editor where users can compose Datagrok functions, queries, and scripts into executable JavaScript pipelines — without writing code.

<!-- TODO: screenshot of the full Flow editor with a sample graph (function browser on the left, canvas in the center, property panel on the right) -->

## Getting Started

1. Open the **Flow** application from the Datagrok sidebar or via **Functions > Apps > Flow**
2. Browse available functions in the left panel and double-click to add nodes
3. Connect nodes by dragging between input/output slots
4. Configure node properties in the context panel (right side)
5. Click **Run** in the ribbon to execute with live visualization

<!-- TODO: gif showing the basic workflow: adding a couple of nodes, connecting them, and running -->

## Features

### Visual Node Graph

Build data pipelines by connecting nodes on an interactive canvas. The graph automatically validates connections using Datagrok's type system — only compatible types can be linked. Nodes are color-coded by slot type for visual clarity.

<!-- TODO: screenshot showing nodes connected with colored links (orange for dataframe, blue for column, etc.) -->

### Function Browser

The left sidebar provides a searchable, groupable catalog of all available nodes:

- **Inputs** — Script input parameters (Table, Column, String, Number, Boolean, DateTime, File, Map, and more)
- **Outputs** — Table Output and configurable Value Output
- **Constants** — Literal values (String, Int, Double, Boolean, List)
- **Comparisons** — Equality, ordering, string matching, null checks
- **Utilities** — Select Column, Add Table View, Log, Info, Warning, type conversions (ToString, FromJSON, ToJSON)
- **Debug** — Breakpoint node for pausing execution
- **Datagrok Functions** — All registered platform functions, queries, and scripts, grouped by role, tags, or package

### Drag-and-Drop

Drag functions, queries, or scripts directly from the Datagrok browse tree onto the canvas to create nodes. Drag files to create pre-configured OpenFile nodes.

### Script Generation

Flow compiles the visual graph into a standard Datagrok JavaScript script with proper `//input:` and `//output:` annotations. The generated script can be:

- Executed directly with live visualization (**Run** button)
- Opened in the Datagrok script editor for further editing (**Run Script**)
- Saved and shared as a `.ffjson` file

The compiler pipeline: **Validation** (cycle detection, required inputs, type checks) → **Topological Sort** (Kahn's algorithm) → **Compilation** (resolve variables and expressions) → **Emission** (generate JavaScript source).

### Pass-Through Outputs

Function nodes automatically get pass-through output slots (marked with `->`) that mirror each input. These enforce execution order for mutating functions — connect a pass-through output to the next node's input to guarantee sequential execution.

<!-- TODO: screenshot showing pass-through outputs on a function node, with an arrow connecting table pass-through to a downstream node -->

### Live Execution Visualization

When running a graph, Flow provides KNIME-style real-time feedback:

| Visual State | Meaning |
|---|---|
| Amber (pulsing) | Node is currently executing |
| Green | Node completed successfully |
| Red (with "!") | Node encountered an error |
| Gray | Results are stale (graph was modified after execution) |

Click on a completed or errored node to inspect its runtime values (DataFrame dimensions, column samples, error messages with stack traces) in the context panel.

<!-- TODO: gif showing a graph being executed with nodes lighting up green one by one, then clicking a completed node to see its output values in the context panel -->

### Debug Mode

Add **Breakpoint** nodes to pause execution at specific points. Run in **Debug** mode to enable breakpoints — execution halts at each breakpoint until you click **Continue** in the ribbon.

<!-- TODO: gif showing debug mode: execution pauses at a breakpoint (amber), user clicks Continue, execution resumes -->

### Property Panel

Select any node to configure it in the native Datagrok context panel. Input nodes expose qualifiers (type filters, semantic types, nullable, choices, min/max, sliders). Function nodes show editable default values for unconnected inputs with rich tooltips describing each parameter.

<!-- TODO: screenshot of the context panel showing a function node's properties with tooltips -->

### Saving and Loading

Flows are saved as `.ffjson` files that capture the complete graph state. Open a `.ffjson` file from anywhere in the platform to restore the flow in the editor.

## Node Types

### Input Nodes

Define script input parameters. Each becomes an `//input:` annotation line in the generated script.

| Node | Type | Configurable Qualifiers |
|---|---|---|
| Table Input | `dataframe` | - |
| Column Input | `column` | type filter, semType filter |
| Column List Input | `column_list` | type filter, semType filter |
| String Input | `string` | nullable, choices, caption, semType |
| Number Input | `double` | nullable, min, max, showSlider, caption |
| Int Input | `int` | nullable, min, max, showSlider, caption |
| Boolean Input | `bool` | nullable, caption |
| DateTime Input | `datetime` | nullable, caption |
| File Input | `file` | nullable, caption |
| Map Input | `map` | nullable, caption |
| Dynamic Input | `dynamic` | nullable, caption |
| String List Input | `string_list` | - |
| Blob Input | `blob` | nullable, caption |

### Output Nodes

| Node | Description |
|---|---|
| Table Output | Marks a dataframe as script output |
| Value Output | Configurable output type (string, int, double, bool, dataframe, column, object, dynamic, and more) |

### Utility Nodes

| Node | Generated Code |
|---|---|
| Select Column | `df.col('name')` |
| Select Columns | `[df.col('a'), df.col('b')]` |
| Add Table View | `grok.shell.addTableView(df)` |
| Log / Info / Warning | Console and balloon notifications |
| ToString / FromJSON / ToJSON | Type conversion helpers |
| Constants | String, Int, Double, Boolean, List literal values |

### Comparison Nodes

`==`, `!=`, `>`, `>=`, `<`, `<=`, Contains, Starts With, Ends With, Is Null

### Function Nodes

Dynamically generated for every registered Datagrok function. Each node's inputs and outputs match the function's signature with proper type coloring.

## Example

A typical cheminformatics workflow:

1. Add a **Table Input** node (or drag a file from the browse tree)
2. Add a **Select Column** node to pick the molecule column
3. Add **Chemical Properties** (from the Chem package) to compute descriptors
4. Add an **Add New Column** node for custom expressions
5. Add a **Table Output** node to return results
6. Connect everything and click **Run**

<!-- TODO: screenshot of the complete cheminformatics example graph described above -->

## Development

### Build

```bash
cd packages/Flow
npm install
npm run build    # grok api && grok check --soft && webpack
```

### Publish

```bash
grok publish           # debug mode
grok publish --release # release mode
```

### Key Dependencies

- [litegraph.js](https://github.com/jagenjo/litegraph.js) `^0.7.18` — Node graph canvas library
- [datagrok-api](https://datagrok.ai/help/develop) — Platform API (external, provided at runtime)

See also: [Datagrok Packages](https://datagrok.ai/help/develop/develop#packages)
