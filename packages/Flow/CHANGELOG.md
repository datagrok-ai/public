# Flow changelog

## v.next

### Scientist-centered UX (Sprint 1)

* **Function browser — task-oriented & decluttered.** Default grouping is now **"what it does"**
  (`categorizeFunc`): **Data Sources** (table out, no table in), **Combine Tables** (≥2 tables —
  joins/links, *no longer* mis-filed as data sources), **Transform Tables**, **Column Operations**,
  **Compute Values**, **Visualize**, **Other** — with **Data Sources first**. All built-in sections
  (Inputs/Outputs/Constants/Comparisons/Utilities/Debug) now start **collapsed**.
* **Catalog exclusions.** ~⅔ of the raw `DG.Func` firehose is hidden: dev/test/internal packages
  (`Dbtests, ApiTests, UiTests, DevTools, Tutorials, ApiSamples, UsageAnalysis`), `test*` functions,
  UI-fragment roles (`editor/cellEditor/panel/widgets/tooltip`), and `funccall` command/dialog wrappers.
* **Start panel (U1).** An empty canvas now shows a welcome overlay — template cards (bundled demos),
  Blank canvas, Open / Import buttons, and a discovery hint — instead of a blank page.
* **Plain-language node status (U5).** Nodes show a short line under the title: *Running… / Done · 1,204 × 8
  / Error / Out of date* (row×col from the captured output summary).
* **Status dot un-overloaded (U7).** Collapse/expand moved to a dedicated caret; the run-status dot is
  now display-only.
* **De-jargoned ribbon (U8).** Menu regrouped to **Flow / Run / Edit / Arrange / Advanced**; script &
  creation-script tools live under **Advanced**; friendlier labels and tooltips throughout.
* **Smarter suggestions (U9).** The drag-out suggestion menu floats common next-step functions
  (join/add-column/aggregate/filter/…) to the top.

* FuncFlow view: Migrated drop handler to the new `doDrop(args)` signature in `ui.makeDroppable`
* Import flows from table-creation scripts: new `Flow:flowFromCreationScript` package function and a
  `File > Import Creation Script...` dialog (prefills from open tables that carry a creation script).
  Variable reads advance along pass-through outputs so graph execution order reproduces the script's
  line order; the first assigned variable is wired to an output node; layered auto-layout.
  Column arguments parse to `ResolveColumn(value, parentTable)`; since the platform `ResolveColumn`
  misbehaves at runtime, they are substituted with the built-in **Select Column** utility
  (`table.col('name')`) — the column name becomes the node's `columnName` and the `table` input is
  wired to the enclosing call's table (or the explicit `parentTable` when present).
  `ResolveColumnList` maps to **Select Columns**. The graph is built by a pure, synchronous
  `buildCreationScriptGraph` (DOM-free) and applied to the editor separately. Imported nodes start
  collapsed and are laid out on a compact wrapping grid (max 4 nodes per row) so the flow fits a view.
* Script order now respects the canvas: the topological sort drains **disjoint subgraphs one at a
  time, top path first** (ranked by topmost node), and picks ready nodes top-to-bottom within a
  component — so a lower path that implicitly consumes an upper path's result runs after it finishes.
* Creation-script import: **removed Output nodes** — each variable's single terminal is now a real
  **`SetVar(variableName, value)`** call (labeled `set: <name>`), registering the value at run time
  under its original name.
* Creation-script import: auto-layout now places **one band per disjoint path, ordered by
  dependency** — a path producing a table sits above the path that reads it via `Select Table` (even
  when defined later in the script), so disjoint paths no longer interleave and the visual top-to-
  bottom order matches execution.
* Select Table emits a tolerant resolver: `tableByName(name) ?? getVar(name)` across the exact,
  no-spaces, and lower-camel name variants.
* Per-function node colors: `FUNC_NAME_COLORS` (type-map.ts) pins a title-bar color by function name,
  checked before role coloring — `SetVar` is now red, `GetVar` light red. Add an entry to pin any function.
* SetVar nodes are now previewable: an instrumented run captures the node's incoming `value`
  (summarized by type), so clicking a SetVar opens the docked output panel showing the stored value
  (table → grid, column → sample, …) even though SetVar declares no output.
* New **Select Table** utility node — resolves an open table by name via
  `grok.shell.tableByName(name)`. The creation-script importer substitutes it for `ResolveTable`
  calls (broken platform-side, like `ResolveColumn`), titled `table: <name>`.
* Creation-script import: replaced the wrapping-grid auto-layout with a connection-aware layered
  layout — build layers become columns (every edge points right), nodes inside a column order by
  predecessor barycenter and greedily align to it, so chains read as straight lanes and branches fan
  out without overlap; column pitch derives from the widest estimated node title.
* Creation-script import: `column_list` arguments (e.g. `JoinTables` keys/values), which parse to
  arrays of `ResolveColumn` calls, now map to a single **Select Columns** utility; numbered params
  pair with the matching table (`keys2` → `table2`). *Every* script variable now gets its own output
  node (previously only the first). Constant nodes are titled after their value (`const: <value>`),
  including when edited in the property panel / inline widget; utility script emission now
  dispatches on the registered node type instead of the user-editable label.
* Fixed: connections attached to a node that starts out collapsed (creation-script import, loading
  an `.ffjson` with collapsed nodes) were invisible until the node was expanded and collapsed again —
  collapsed nodes only render socket DOM for connected sockets, so the editor now re-renders
  collapsed endpoints whenever one of their connections is created or removed.
* Tests: added a Flow test suite (`grok test`) covering type compatibility, the node factory,
  topological sort, the script emitter, validator, serializer round-trip, and creation-script import
  (including the chem-properties example end to end).

## 0.0.1 (2026-03-08)