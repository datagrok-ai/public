# Flow changelog

## v.next

### Edit string-list inputs inline

* `string_list` (and its `list<string>` spelling, which DG folds to the same thing) inputs are now
  **editable in the context panel** as a comma-separated text field — like column lists. The compiler
  trims each entry, drops blanks, and passes a real JS array to the function (`"a, b ,, c"` →
  `["a", "b"]`... → `["a", "b", "c"]`); an empty field is omitted so the function keeps its own
  default. Round-trips through creation-script import/emit, and a wired `String List Input` / `List`
  node still works for the connected case. Plain `list` (which may hold non-strings) is left as-is.

### Easier wiring — compatible-target highlighting + drop-on-node

* **Drag from a pin to see where it can go.** While you drag a connection from a socket, the canvas
  dims and only the sockets (and their nodes) that can legally accept it — compatible type, opposite
  side — light up green. Works from an output pin (lights compatible inputs) and from the tail of an
  existing connection (lights compatible outputs); the drag origin stays bright.
* **Drop on the node, not just the tiny pin.** When you drop an output drag anywhere on a node that
  has a single compatible, unwired input, Flow connects to it — no more pixel-hunting for the dot.
  Nodes with zero or several candidates still require aiming at a specific pin; empty-canvas drops
  still open the suggestion menu. Bonus: this works even on collapsed nodes (whose pins aren't drawn).

### Pick columns from a list — no typing names from memory

* Every `column` / `column_list` input on a function node now has a **picker icon** in the context
  panel (next to the field). Clicking it opens a real **column / columns dialog** seeded by the
  *actual* upstream table, so you choose from a list instead of recalling column names. For
  multi-table funcs (Join Tables: `keys1`/`keys2`/`values1`/`values2`) each column resolves against
  its own table input (`keys1`→`table1`, `keys2`→`table2`), and the picker uses the right one.
* The picker now extends to **every** node with a column field, not just DG functions: **viewer
  nodes** (X/Y/Color/Size axis columns, picking from the wired table) and the **Select Column /
  Select Columns** utilities. One shared `createColumnFieldRow` drives all of them — any node with a
  column-valued field plus a dataframe input gets the picker for free.
* **Fixed:** the picker (and "inspect anywhere") failed for a table reached through a node's
  passthrough output when that node's *real* output isn't a table — e.g. a Scatter Plot wired to
  AddNewColumn's "table →" passthrough (AddNewColumn returns a *column*). The instrumented run now
  captures the threaded, post-execution table (`<input> (modified)`) whenever the node has no real
  dataframe output, so the picker can read it instead of erroring with "no table produced".
* Three cases, handled automatically: the table input isn't connected → a hint to connect one; it's
  connected and already run → pick immediately from its captured output; connected but not yet
  computed → offer to **run the flow up to that point**, then pick from the produced table
  (`ExecutionController.produceTableForNode` runs a headless slice and returns the live `DataFrame`).
* The slice run is **additive** (`preserveState`): it only (re)computes its own slice and leaves every
  other node's captured result intact — so picking `keys1` for table1 and then `keys2` for table2 no
  longer wipes table1, and a later pick against the same table reuses the cached result instead of
  re-running it. Reuse is limited to *fresh* (completed, non-stale) results — a graph edit still forces
  a recompute.
* New **how-to: “How do I add visualization nodes?”** — loads demog and wires it into Scatter Plot,
  Bar Chart, and Pie Chart, changing columns along the way.
* New **how-to: “How do I join two tables?”** — adds two demog tables and an output, wires up Join
  Tables, and uses the column picker to fill `keys1`/`keys2` (both `USUBJID`) and `values1`/`values2`
  (all columns) without typing. The two source tables get distinct, ordered socket highlights
  (`byNodeFuncNth`), and the column-pick steps **re-anchor the hint card to the open dialog**
  (`preferDialog`/`openDialogEl`) — previously the dialog was hidden behind the card (z-index 3000 vs
  5000). The guide card now also drops below any open platform dialog so it can never block it.

### Node colors by what they do

* Function nodes with no DG role (the gray majority — Join Tables, Add New Column, chemical
  properties, …) now take a **title-bar color from their task category** (Data Sources / Combine /
  Transform / Column Operations / Compute Values / Visualize / Other) via `CATEGORY_COLORS`, so the
  canvas reads at a glance. Precedence is unchanged where it matters: a pinned function color
  (`FUNC_NAME_COLORS`) > an explicit role color (`ROLE_COLORS`) > the category color > gray. The
  browser's `categorizeFunc` and the coloring now share one `categorizeBySignature`.

### Viewers & Widgets — first-class visualization in the toolbox

* **Viewers pane** with **manually-built viewer nodes** that don't need a TableView lifecycle. Wire a
  table into a viewer node, run, and the live `DG.Viewer` renders in the preview panel. Core charts
  (Scatter Plot, Histogram, Line/Bar/Pie/Box, Heat Map, Grid, Trellis, Network) come from the
  `DataFrame.plot` namespace; non-core package viewers (Radar, Sunburst, Word cloud, Scaffold Tree, …)
  are discovered via `DG.Func.find({meta:{role:'viewer'}})` and built generically. Emits
  `let v = await table.plot.fromType('<Type>', {}); v.setOptions(<look>);`.
* **Edit any setting, persistently.** Each viewer node exposes a few high-value options (X/Y/Color/Size
  columns, title) in the context panel; for everything else, click **“Edit settings”** on the preview —
  Flow does `grok.shell.o = viewer` (Datagrok renders the full settings editor) and **captures every
  change back onto the node** (debounced `onPropertyValueChanged` → `getOptions().look` minus `#type`),
  so a re-run reproduces the exact look.
* **Widgets pane** collects every function that produces a `widget` (info panels, search widgets, …),
  grouped out of the categories.
* **The default viewer *functions* are gone** from the catalog — they required a TableView and are
  replaced by the viewer nodes above.
* **Clear-search ✕** at the right edge of the toolbox search box.
* **Big catalog cleanup:** functions whose I/O is **only scalars** (string/number/bool/dynamic — ~250
  math/string helpers) are hidden, along with viewer/view-producing functions. The catalog drops from
  ~786 to a far more navigable set focused on data-flow steps.

### Visualize results — widgets & viewers in the preview panel

* **A node whose output is a `widget` or `viewer` now renders live** in the bottom preview panel.
  The instrumented run captures the actual object by reference (same in-tab mechanism as the
  DataFrame clones), and clicking the node mounts its `.root` in the docked panel.
* **Context panel** no longer shows `[object Object]` for a widget/viewer output under Execution —
  it now reads plainly `widget` / `viewer` (the live preview is in the docked panel).
* **The docked preview behaves like it belongs to the view:** switching to another view closes it
  (via `grok.events.onCurrentViewChanged` + a `grok.shell.v === this` check), and if you manually
  close it then click another previewable node (results still fresh), it **reopens** — the panel
  re-checks whether its dock is still alive before reusing it.
* **“Run up to here & preview” is now on the node's right-click menu**, not just the output port's —
  the more discoverable, intuitive place for it.
* **Catalog cleanup:** functions whose output is a whole **view** are filtered out (they can't be
  composed or previewed), and the **Comparisons** built-in group is hidden from the toolbox for now
  (still registered, so existing flows keep loading).

### KNIME-style Files & Queries in the toolbox

* **Files pane (open by default), first in the toolbox.** A real file browser (the shared
  `getFilesBrowser` tree) listing every file connection and its folders/files. **Drag a file onto the
  canvas** — or **double-click it** — to drop an `OpenFile` node already pointing at it (no path
  typing). Its expanded folders and scroll position survive a search keystroke (the tree is built once
  and reused).
* **Queries pane, grouped by connection.** Every `DG.DataQuery` is collected here into one
  sub-accordion **per data connection** (`connection.friendlyName ?? connection.name`) — Chembl,
  Chembl Sql, Unichem, Biologics, … — **regardless of the “Group by” mode**. Queries are **removed
  from the Data Sources category** so the catalog isn't drowned by hundreds of DB queries. The
  per-connection sub-accordions are visually nested (indent + left accent rail + lighter headers +
  card shadow) so the hierarchy reads clearly as a level below the top sections.
* **Name-based test-ids on every tree row** — `ff-files-conn-<name>` / `ff-files-folder-<name>` /
  `ff-files-file-<name>` (+ raw `data-conn` / `data-folder` / `data-file` / `data-file-path`), and
  `ff-browser-files` / `ff-browser-queries` / `ff-browser-query-conn-<name>` on the panes — so a
  specific connection, folder, or file (e.g. the Demo connection's demog.csv) is addressable from
  guides and UI tests.
* **Guides now teach the real data-loading gesture.** The flagship tutorial and the
  “How do I bring data in?” / “preview a node's data?” how-tos no longer paste a path string — they
  walk the user through the Files browser: **open Files → expand the Demo connection → scroll to
  demog.csv → double-click / drag it**, each step properly gated (`untilFileTreeConnExpanded`,
  `untilScrolledIntoView`, `untilFuncNodeWithInput`). Prerequisites are detected and skipped when
  already satisfied. While scrolling, the highlight tints the **whole Files pane** (not the
  still-off-screen target row) and snaps to the file only once it's actually in view.

### New tutorial — “Tour the interface”

* A guided walkthrough of **every UI control**: each toolbox pane (search, group-by, Files, Queries,
  built-in blocks, function categories), each ribbon group (run/debug/stop, view-script, save/open,
  undo/redo, layout, zoom, toggle-toolbox, help), the canvas + a node's anatomy (title, caret, status
  dot, sockets), the overview and status bar, and the context panel (title row, type badge,
  parameters). Mostly read-and-Next steps; selecting a node is interactive. A sample Table Input node
  is **auto-added as a prerequisite** so the canvas/panel sections always have something concrete to
  point at.

### Auto-summaries — the flow documents itself (U12)

* **Every node shows a plain-language caption** when you haven't written your own description —
  "Loads file demog.csv", "Adds column “AgeMonths” = ${AGE} * 12", "Computes MW, HBA, logP",
  "Calculates logP", "Joins two tables (left)". Generated heuristically (`summarizeNode`) from, in
  order: a built-in type summary, a **curated** function summary, or a humanized fallback from the
  function's friendly name — so unknown functions still read as words, not `camelCase`.
* **Curated definitions** in [`src/summary/summary-defs.ts`](src/summary/summary-defs.ts): ~80
  high-traffic functions (core transforms/joins/columns, Chem properties/space/clustering, Bio
  alignment/space/regions) + Flow's built-in nodes, chosen empirically from the live 800-function
  catalog. Templates read the node's own inputs (file path, formula, enabled property flags, join
  type, …) so the caption reflects what *this* node does.
* **"Describe this flow…"** (ribbon → Advanced) summarizes the whole canvas, **grouping disjoint
  subgraphs into separate numbered pipelines** (`summarizeFlow` over connected components). Each
  pipeline is numbered in the **exact order the script executes** — it shares the compiler's
  `topologicalSortNodes` (extracted as a pure core of the emitter's own sort), so the description and
  the generated/run script are guaranteed to agree line-for-line — and every step lists **what feeds
  it** — e.g. *"Joins two tables (inner) (← table1 from step 11, table2 from step 9)"* — so the
  description reads as a real recipe, not a flat unordered list. The dialog renders each
  pipeline as a card with a numbered, full-width (untruncated) step list. SetVar/GetVar plumbing nodes
  get informative captions ("Stores result as “SPGI”").


### Guide system — concrete, well-placed, hands-on

* **Tutorials are now concrete and end-to-end.** The flagship tutorial **“Load data and add a
  column”** has the user *search* for each function (case/space-insensitive — "open file" or
  "openfile"), opens a real dataset (`System:DemoFiles/demog.csv`, pasted from the clipboard into the
  Open File node's **File path** field), drags **Add New Column** clear so it stops overlapping
  (gated until it's actually moved right), **wires specific pins** (Open File's `result` → Add New
  Column's `table`, both pins highlighted), names the column (`My New Column`, gated on that exact
  text), pastes the formula `${AGE} * 12`, adds a **Table Output** and connects Add New Column's
  `table →` pass-through to it, then runs. Every step highlights a *specific* element (browser item,
  canvas node, socket pin, `data-param` input row, ribbon icon), never a vague "the canvas". All
  targets/params/socket keys/the demo file were verified empirically against a live server.
* **Steps can highlight several elements.** A step's optional `highlights` returns multiple targets,
  each getting its own pulsing dot — used to light up *both* pins the user must connect.
* **The other tutorials are concrete too.** *Find the right function* now adds **Join Tables** and the
  final step highlights its `result` output and is gated on the user actually picking a function from
  the drag-out suggestion popup (was an ungated "tip" highlighting the whole canvas). *Organize your
  canvas* now adds two specific nodes (Table Input + Table Output), has the user **drag them apart so
  they don't overlap**, then collapse, tidy, **undo, redo**, navigate the overview, and zoom to fit.
* **How-to answers now have prerequisites and action-based gates.** Each "How do I…?" first ensures
  its preconditions are met — a step with `skipIf` silently adds the needed node(s) only when they're
  missing (e.g. *How do I collapse a node?* adds a Table Input first; *How do I connect two nodes?*
  adds a Table Input + Table Output and highlights the two specific pins; *How do I preview?* adds an
  Open File and pastes the demo path so there's real data to preview). Gates now detect the *action*,
  not pre-existing state: collapse waits for *another* node to fold (`untilMoreCollapsed`) instead of
  resolving instantly if any collapsed node already existed. Highlights point at the exact element
  (the specific node, pin, caret, or input row).
* **The "All done" note auto-dismisses after 5s** (and on its Done button) so it never lingers.
* **No lingering highlight.** A step whose condition was already satisfied could leave its orange
  highlight stuck (a queued frame re-applied it after cleanup). Fixed at the source (guarded/cancelled
  re-anchor frame) and hardened with `GuideRunner.clearAllHighlights()` on every step boundary and on
  finish/exit.
* **Function search matches the real name, case- and space-insensitively.** The toolbox shows the
  friendly name ("Open File"), but searching "OpenFile"/"openfile"/"open file" — or "tableoutput" for
  the built-in "Table Output" — now all match (`funcMatchesSearch`/`nameMatchesQuery`). Previously
  "OpenFile" returned nothing because the list only matched the friendly name.
* **Smart popup placement.** The instruction card is now our own element with a pure, unit-tested
  placement function (`computePlacement`): it honors the preferred side when it fits, else flips to
  the opposite, else picks whichever side has room, and **always clamps fully on-screen** — so a hint
  on a bottom element opens *above* it, one on a far-left element opens to its *right*, and nothing
  spills off the viewport. The card re-anchors on a timer (nodes re-render, the context panel shifts
  layout, the user drags), and the old `scrollIntoView` that yanked the whole UID upward is gone.
* **Single close affordance.** The platform `ui.hints.addHint` injected its own ✕ that overlapped our
  Exit link; the card now owns one ✕ (top-right) that exits the guide. Targets get a clear pulsing
  outline (`.ff-guide-target`) instead of the stray corner blob.
* **Launch the tour from the Start panel.** A primary **“Take a 2-minute tour”** button on the empty-
  canvas Start panel runs the flagship tutorial.
* **More how-to answers**, including context-panel ones (“How do I add a calculated column?”, “How do
  I edit a node's settings?”, “How do I bring data in?”), each highlighting a concrete element.

### Empty-canvas overview

* The **overview minimap is hidden on an empty canvas** (nothing to overview) and reappears the moment
  the first node lands (`FlowEditor.refreshMinimap`, driven from `onGraphChanged`).

### Interactive tutorials & how-to help

* **Guide system.** A non-invasive floating help button (bottom-left) and a ribbon icon open a menu
  of **4 multi-step tutorials** (Build your first flow · Find the right function · Organize your
  canvas · See & reuse the generated script) and **14 how-to answers** (add/connect/run, preview a
  node, change settings, delete, layout, categories, collapse, undo, save, view script, open,
  navigate, bring data in). Each step **highlights the relevant element** (via `ui.hints`) and
  **waits for the real action** — a click, a typed value, a node added/connected/collapsed — before
  advancing, with Skip/Exit always available. Built on the new `data-testid`/`data-*` attributes.
  Lives in [`src/guide/`](src/guide).

### Testability

* **`data-testid` on every UI surface.** The canvas (nodes, sockets, status dot, caret, hints,
  connections, edge counts), ribbon, start panel, function browser (sections + items), property
  panel (per-input rows), minimap, suggestion menu, and output/value previews all carry stable
  `data-testid` attributes for deterministic UI testing. Dynamic ids are generated from what the
  element *is* (function name, param name, socket key) via [`utils/test-ids.ts`](src/utils/test-ids.ts)
  (`tid`/`setTid`). See the **Test IDs** table in CLAUDE.md.

### Scientist-centered UX (Sprint 2)

* **Inspect anywhere (slice-compile).** Right-click any output port → **"Run up to here & preview"**:
  Flow compiles and runs just the slice needed to produce that node (the node + its upstream
  ancestors, `sliceUpTo`), then opens its data — no full run and no Output node required.
  `EmitOptions.onlyNodeIds` restricts emission to the slice; the controller focuses the node on
  completion.
* **"Needs input" hints.** Each node now shows a plain-language amber hint listing the structural
  inputs (a table, a column) still to be connected or filled — continuous, pre-run, on the node
  itself (`missingRequiredInputs` + `FuncNode.requiredInputs`). The status dot turns amber too.
* **Row counts on wires.** After a run, each data connection is labelled at its midpoint with the
  row/value count flowing through it (Make/n8n-style), cleared on edit / re-run.
* **Function browser remembers its state.** The group-by mode and which sections you expanded are
  persisted to `localStorage` and restored next time (`funcflow.browser.v1`).

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