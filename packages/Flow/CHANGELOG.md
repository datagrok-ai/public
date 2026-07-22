# Flow changelog

## v.next

* Toolbox: restructured into clear zones — a global **Search** on top that filters everything (the function catalog AND the Queries/Workflows/Favorites tabs, whose headers show live match-count badges while 0-match tabs dim and the Suggestions strip collapses to a dimmed header so it isn't mistaken for results), then the collections tab strip, then a labeled **Functions** catalog whose "Group by" became a compact "by: …" popup button on the zone header, with the Suggestions strip at the bottom; fixed the tab pane's vertical scrolling (the platform tab host's fixed 400×300 min-size and flex-shrink:0 prevented both sizing and scrolling inside the narrow toolbox), made the pane the single scroll container (the files tree no longer nests its own scrollbar), and made the tab strip hold up on narrow/zoomed screens (labels ellipsize before a count chip can clip; a container query compacts headers and chips on tight widths)
* Toolbox: the collections tab strip — **Files**, **Queries**, **Workflows** (saved flows), and the new **Favorites** — persists the selected tab, hugs its content height (capped at a third of the panel, then scrolls), and the Files tab leads with a "drag a file from your computer" hint
* Toolbox: introduced **Favorites** — hover any node row and click its ★ to pin the node into the Favorites tab (stored in localStorage); starred rows show a gold star everywhere, and favorites can be added/dragged/unstarred straight from the tab
* Toolbox: redesigned the Suggestions pane — a clearly distinct assistant surface (accent top edge, tinted background, lightbulb header) instead of blending into the catalog; suggestions render as thin one-line rows (action + inline ellipsized reason) so more fit, the pane takes 25% of the toolbox height instead of 30% and shrinks to its hint when empty, and the hint speaks plain language ("Click a step on the canvas…")
* Toolbox: searching "chart", "plot", "graph", or "viewer" now surfaces the whole Viewers section, so charts are findable by the words scientists actually type
* Toolbox: clearer zone separation — the collections tabs sit on a grey tray with a soft recess under the white Functions catalog; every section header (and each query-connection row, via a uniform database glyph) got a muted Font Awesome icon, so the Queries tab no longer reads as more function categories; the Suggestions lightbulb warmed to a dimmed orange

* Local file upload: drop a file from your computer onto the canvas (or several at once) and it becomes a replayable, shareable **Uploaded File** node — bytes stay in memory until the flow is saved, then upload to the server's GUID-addressed file store and the node references the stable id; saving/sharing the flow grants read on the files to the same people — share handling runs globally via the `flowShareSync` autostart, so sharing a flow from Browse with no Flow view open works too; 100 MB per-file cap; CSV/TSV/TXT/d42 parse natively, other formats go through the platform's registered file importers (the platform's "drop to open" overlay is suppressed while a Flow view is active via the new `grok.events.onFileDragEnter`)

* Function catalog: inverted the inclusion logic from deny-based (structural filters + role/tag rules + a curated denylist) to an explicit **allowlist** (`INCLUDED_FUNC_NQNAMES`, frozen from the exact catalog the old filters produced on a live server) — any other function now needs `meta.includeInFlow: true` to appear in the toolbox; saved flows and queries are always included by kind (queries from dev/test packages like Dbtests excepted), and `meta.includeInFlow: false` still opts out unconditionally

* AI: the flow view exposes its operations to the AI assistant as registered package functions (`flowViewFunction` tag, returned from `FuncFlowView.getFunctions()`) — list/inspect nodes and connections, search the curated node catalog with filters (text query, accepted/produced DG type), add nodes with prefilled inputs, connect ports with type checking, set node input values, select a node, and run the flow with validation and per-node failure reporting; each function takes the generic view and acts on `view.jsView`
* AI: guides — `listFlowGuides` / `startFlowGuide` let the assistant find a matching interactive tutorial or "how do I…" walkthrough and launch it (after asking the user) instead of explaining in text
* AI: the view sets `aiDescription` — a briefing shown to the assistant with every prompt naming the Flow view functions to use

* Node groups: select several nodes and Group them (Ctrl+G or the context menu) into a titled, collapsible frame — dragging the frame (title bar or background) moves all members together; the caret (or a double-click on the card) minimizes the group into a single collapsed-node-like card that shows the title, description, and an aggregate run-status dot (running / done / error / out of date across the members), hides the member nodes and their internal wires, and re-anchors boundary wires to dots on the card's edges; maximizing restores everything — wherever the card traveled, the contents follow; groups are purely visual (the graph, compilation, and execution are untouched), persist in the `.ffjson`, and offer Ungroup (Ctrl+Shift+G), Tidy layout (re-arrange just the members), Remove from group on a member node, rename/description editing via double-click, and Delete group and nodes
* Annotations now carry their contents: dragging an annotation moves every node whose center sits inside it — plus annotations fully contained in it and the waypoints of wires between carried nodes; drag a node out of the frame and the next annotation drag leaves it behind
* Start panel: added a "Recent flows" section — the 10 most recently updated flows visible to the user (own + shared), newest first, rendered as one-line entity rows like in the Browse panel; clicking a row opens the flow; the list refreshes whenever the start panel reappears
* Copy/paste and multi-node duplicate: Ctrl+C copies the selected nodes (with the connections among them) into the editor clipboard, Ctrl+V pastes them offset and selected (repeat pastes fan out); the context menu's Duplicate now duplicates the whole multi-selection with its internal connections — and right-clicking a selected node no longer clears the multi-selection; duplicated outputs/SetVars get a unique variable name instead of an instant duplicate-name error
* Escape now clears the node selection
* Preview panel: added a pin — pinning freezes the preview on the shown node so clicking other nodes (to adjust their settings) no longer replaces it, while fresh results of the pinned node re-render live (e.g. pin a chart, tweak an upstream formula with autorun on, watch the chart update); the pin survives invalidation and re-runs, drops when the node is deleted or the graph is replaced, and unpinning jumps the preview to the currently selected node; while the pinned node recomputes the panel no longer hides and re-docks — the stale content stays under a "Recalculating…" indicator and the fresh result swaps in place
* Introduced integrated output views (Spotfire-style pages): the status bar leads with a tab strip — Canvas plus one tab per table output (Table Output, dataframe-typed Value Output, and SetVar terminals carrying a table) — and each tab hosts a full detached TableView of that output (ribbon, toolbox, viewer gallery via PowerPack, dockable viewers); tabs are named after the table, show a grey/green/amber state dot (empty / ready / out-of-date), display a centered "run the flow" message with a Run button until a value lands, refresh in place on re-runs (viewers survive by column matching), and shrink with ellipsized labels when crowded; switching to a table tab swaps the view's ribbon panels and toolbox to the TableView's (led by Flow's Save pill, since saving the flow is what persists the tab layout; the TableView's own core-hidden Save is dropped instead of copying as an empty panel) and Canvas restores Flow's own
* Output-view tab layouts (viewers, docking, grid settings) save with the flow (a new optional `outputViews` section in the `.ffjson`, keyed by the output's paramName so layouts survive node-id remapping) and re-apply when the tab is next opened with a value; layout-only edits don't light the Save button (v1 limitation — the layout is still captured on every save)
* Save now always opens the combined save dialog: script name/description/space plus a run-aware Dashboard section — if the flow hasn't produced tables yet it offers a "Run the flow" button (the section refreshes when the run ends); with computed outputs, saving can chain straight into the platform's **standard Save-project dialog** seeded with the output tables and their tab layouts (data sync toggles, creation-script dependencies, share link, upload — all from core via the new `DG.Project.showSaveDialog`); each table's creation script calls the saved flow (`param = Flow(...)`), with the output accessor (`.param`) appended only for multi-output flows, so a data-synced dashboard replays the flow once on open and every table binds its own output
* The published dashboard project is bound to the flow (persisted in the saved script): saving again updates the SAME project instead of creating a new one per save — a "publish as new" link in the save dialog breaks the binding; and a tab's saved layout now ships with the dashboard even when the output view was never opened that session
* Fixed the output-view tab ribbon disappearing on the second switch to the same tab (`setRibbonPanels` moves elements — the TableView's ribbon content is now captured once as stable inner-element references and re-set on every swap)
* Updated the in-app tutorials and how-tos for output views and dashboard publishing: the flagship tutorial and interface tour now cover the result tabs, the save-related steps describe the combined save dialog, and a new "Publish your results as a dashboard" tutorial plus four how-to answers (open a result as a full table, layout persistence, publish a dashboard, update/unbind a published dashboard) were added
* Rewrote the README for non-technical users — a feature walkthrough (building blocks, suggestions, running, result tabs, glass-box script, saving, dashboard publishing) illustrated with 16 real screenshots captured from a live server (`img/`, excluded from publishes via `.npmignore`)
* Wire data-count labels now read "N × K" (rows × columns) instead of "N rows", and appear on every table-carrying wire — including pass-through (`table →`) ports and utility nodes like Select Table, which previously never got a count
* Re-clicking or grabbing the already-selected node no longer re-fires selection events, so the context panel and the suggestions pane stop rebuilding on every click; the suggestions pane also skips re-rendering when a recompute returns the identical suggestion set

* Introduced the Outputs strip: a thin vertical "Outputs" column outside the canvas viewport (graph content can never pan or fit behind it) that owns every flow output — Table/Value Output nodes render as small fixed-size chips inside the strip (one-letter type; the tooltip names the output, its type, and what feeds it) instead of free-floating cards; chips are screen-space, so zoom and pan never move or resize them, with wires plugging into the canvas edge at the chip's row; any output socket dragged onto the strip publishes that value as a flow output (auto-named, auto-typed), and the strip lights up as a drop target during output drags; the data model, serialization, and compilation are unchanged, so existing flows load as-is with their outputs docked
* Context panel: Title, metadata chips, and Description now share one aligned header block; the node and function descriptions are combined into the single editable Description input (function text is the fallback)
* Context panel: Replaced the Function pane with compact chips (full name, package, roles, tags); the parameters pane is titled with the function name and expanded by default
* Order ports: now hidden unless wired, hovered, or mid-order-drag (with pointer-events disabled while invisible), and carry plain-language tooltips explaining run-order wiring
* Order drags: dimming + green highlight of every other node's opposite order square while dragging, and dropping on a node body auto-connects its order port (no aiming at the small square)
* Humanized caption-less parameter names everywhere (node slot labels, context-panel rows, connection endpoints) mirroring core's propertyNameToFriendly — maxNumOfSomething reads "Max Num Of Something", matching ui.input.forProperty; declared captions still win; all-caps acronyms (MW, HBA) are preserved instead of core's Mw/Hba folding
* Node title bars are now pastel (identity hue mixed 60% with white via `pastelize`) so nodes blend with the canvas; vivid hues stay on small surfaces (minimap, sockets, connections)
* Node identity colors now come from the platform's categorical palette (`DG.Color.categoricalPalette`) instead of bespoke Material hues, so nodes wear the colors users know from categorical coloring across Datagrok
* Live-by-default nodes: Open File, Add New Column, and all viewers now run automatically — on change, on a ready drop (file dragged onto the canvas), and on flow load — even with the autorun toggle off; only the listed nodes run, each gated on satisfied inputs, never the rest of the canvas (editable lists `AUTORUN_FUNC_NAMES` / `AUTORUN_NODE_TYPES` in execution/autorun.ts)
* Socket dots replaced with Column Manager-style type-letter chips (`t` table, `s` string, `i` int, `d` double, …) — column-data types use core's exact `Color.typeColors` letter+color pairs, unknown types get a gray first-letter chip; hover/connect animations and pass-through dashing preserved
* Suggestions pane: the bottom of the toolbox (minimizable via caret) now ranks the ~10 most likely next steps from the full canvas context — detected semantic types in a node's data suggest domain operations (Molecule columns → Chem descriptors/properties, Macromolecule → Bio To Atomic Level/MSA), a molecule clicked in the output preview suggests similarity/substructure searches prefilled with that value, selecting two tables suggests Join Tables wired to both, and a single table offers the common next steps, matching viewers, and outputs; suggestions add like toolbox nodes — double-click, or drag onto the canvas (wiring and prefills apply either way)
* Selection now matches the platform's selectRows conventions: Shift+drag draws the area-select marquee (was Ctrl+drag) that adds the covered nodes — Ctrl at mouse-up removes them; Shift+click adds, Ctrl+click toggles, Ctrl+Shift+click removes a node; Ctrl+A / Ctrl+Shift+A select / deselect all; the marquee wears the platform look (black stroke, light-blue fill)
* Suggestions: with two tables selected, the combiners (Join Tables & co) now rank above per-table suggestions like cheminformatics matches
* Required inputs now cover every non-optional parameter without a declared default (not just tables/columns): the node shows "Needs input" until the value is connected or filled, and no run — live-by-default, autorun, full run, or rerun — executes it before that; fixes a bare Open File live-running with no `fullPath`
* Context panel: Connections pane shows only wired slots — each row naming the node and slot at the far end (`table ← Open File · result`), order edges as "runs after/before" facts — flags missing required inputs/values in amber, and auto-expands when something is missing
* Introduced per-function input overrides (`HIDDEN_FUNC_INPUTS` / `CUSTOM_FUNC_INPUT_EDITORS` in utils/func-input-overrides.ts): Add New Column's `subscribeOnChanges`/`errorBehavior` no longer clutter the node or the context panel (still compiled and round-tripped in creation scripts), and Open File's `fullPath` is now edited with a proper file input instead of a plain text field
* Introduced function wrappers (`FUNC_WRAPPERS`): a function with an awkward signature can expose reshaped, Flow-friendly node inputs that fold back into the real arguments at compile time — Append Tables (previously unwirable: it takes a `dataframe_list`) now exposes two plain table inputs, joins the two-table suggestions, and auto-wires

### Toolbox sections are now platform accordions

* The function browser's collapsible sections (Files, Queries, Viewers, Widgets, the built-ins, and
  the function categories) are now a standard **`ui.accordion`** instead of hand-rolled collapsible
  divs — platform look & feel, lazy pane content (a collapsed section never builds its items), and
  **self-persisted expand/collapse state** (`localStorage['Accordion:funcflow.toolbox']`, keyed by
  pane name; the per-connection query groups are a nested accordion under
  `Accordion:funcflow.toolbox.queries`). Category counts render as the pane's count badge instead of
  being baked into the title.
* While a search is active every matching section is still forced open — on an **unkeyed** accordion,
  so the forced states never overwrite the user's saved ones.
* Pane headers keep `data-section` and the `ff-browser-section-<title>` test ids; the guide now
  detects expansion via the accordion's `expanded` header class.
* **Every section (except Files) now shows its item count** — Queries, Viewers, Widgets, and the
  built-ins got count badges too, matching the categories.
* **Section order reorganized**: Files, Queries, Cheminformatics, Bioinformatics, the task
  categories (Data Sources, Workflows, Combine Tables, Transform Tables, Column Operations,
  Compute Values), Viewers, Widgets, Inputs, Outputs, Constants, Utilities, Other — and Debug last.

### Save button reflects whether there's anything to save

* The ribbon **Save** button greys out **and becomes non-clickable** (`pointer-events: none`, steel
  colors) when there is nothing to save — an **empty canvas** ("Nothing to save yet — the canvas is
  empty") or **no change since the last save** ("No changes to save since the last save"). Its
  wrapper carries the "why" tooltip so it still shows while the button is disabled. It re-enables the
  moment the graph differs from the saved baseline.
* Change detection compares a serialization of the graph + settings with the volatile
  `created` / `modified` / `author` fields stripped — `serializeFlow` stamps fresh timestamps on
  every call, so leaving them in made every comparison report "changed" (Save never greyed, and the
  same flow could be saved back-to-back). Undoing back to the saved state now disables Save again.
* The state refreshes on **both** structural changes (`onGraphChanged`) **and parameter edits**
  (`onGraphEdited` — `notifyNodeParamsChanged` reports only there, so editing an input value now
  updates Save). Baselines are recorded on save, on opening a server-backed flow, and for a fresh
  empty flow.

### Auto-pin the preview view so the toolbox appears

* A flow app / flow script opens as an unpinned **preview** view, and the platform hides the toolbox
  until a view is pinned. Flow now **pins itself on the first interaction** (a click on the canvas,
  an edit, a keypress) when it is the current view — so the toolbox shows up without the user having
  to pin manually. One-shot: the listener removes itself once it fires (and on `detach`).

### Preview shows every renderable output, side by side

* A node that surfaces **more than one renderable value** now previews them all, laid out side by side
  with draggable vertical dividers (`ui.splitH`) instead of showing only the first. Covers a
  **multi-output func** (e.g. two dataframes) and an in-place **mutator with no output but several
  modified input tables**. A lone value fills the panel exactly as before.
* The instrumented run now captures **every** connected dataframe input as a `"<input> (modified)"`
  summary when the node has no dataframe output (previously only the first) — so a two-table mutator
  previews both tables it transformed.

### Creation-script import fixed for the new SetVar signature

* Importing a creation script (and everything downstream: per-table split, order-edge inference,
  auto-layout) worked again. The platform's `SetVarFunc` grew optional `outputName` / `outputIndex`
  parameters (multi-output binding), so a parsed assignment (`X = f(...)`) now surfaces up to four
  inputs. The importer's `asAssignment` gated on **exactly two** inputs, so it silently rejected
  every assignment — no variables were registered and the whole import collapsed to nothing. It now
  validates the `variableName` / `value` inputs directly and ignores the extra optional params.

### Multi-output funcs compile correctly

* A func node with **more than one output** now compiles to valid, runnable code. `grok.functions.call`
  returns the value directly for a single-output func but an **object keyed by the declared output
  names** for a multi-output one — so each output is now read as `<callVar>.<outputName>`. Previously
  the emitter referenced per-output variables (`<var>_result1`, `<var>_result2`) that were never
  declared, so any downstream node wired to the 2nd+ output got `undefined` (or a ReferenceError).
* **Variable names are always valid identifiers.** `toCamelCase` keeps digits, so a func/node named
  e.g. "2 Inputs Flow" produced `let 2InputsFlow… = …` — an immediate syntax error. `uniqueVarName`
  now prefixes a leading digit with `_`.
* The instrumented per-node output summary is now keyed by the **output slot key** (a valid object
  key, and what `labelOutgoingConnections` / the live-stash look up) instead of the value expression
  — which for a multi-output func is `<var>.<name>` and can't be an object-literal key.

### Save verifies the bound entity actually exists

* **Save** no longer silently writes to a bound script id that isn't really on the server. A flow can
  carry a `boundScript` with an id yet have no server entity — a template opened as a new flow, or a
  flow whose entity was deleted elsewhere — and Save used to update that phantom id under whatever
  name the template carried. It now confirms the entity resolves (`grok.dapi.scripts.find(id)`, via
  `scriptExistsOnServer` — a throw *or* a nullish result counts as "not saved") before a silent
  update; otherwise it opens the **Save As** dialog to ask for a name.

### Minimap minimizes when the output preview opens (restored)

* Opening the bottom **output preview** again auto-minimizes the overview **minimap** to its header
  strip, so the two don't crowd the same corner. This one-shot behavior (fired only on the panel's
  hidden → visible edge; reopening the minimap by hand while the preview stays open sticks) was
  dropped in the dock → splitter rework and is now wired back in `FuncFlowView` via
  `OutputPreviewPanel.onStateChanged` → `flow.setMinimapCollapsed(true)` (the live minimap only — the
  remembered initial state is untouched).

### Column picker: menu instead of a second dialog

* The column / column-list picker (the list icon next to a `column` / `column_list` input) now
  drops a **`DG.Menu` column selector right under the icon** instead of opening a modal dialog.
  When the upstream table hasn't run yet, the run-confirmation dialog still appears; on **OK** the
  selector menu pops immediately — one interaction, no back-to-back dialogs. Single inputs use
  `singleColumnSelector` (click a column → set → close); lists use `multiColumnSelector`, writing the
  checked set back live on each toggle.
* **Type / semantic filtering.** When the input is a DG func param carrying a `semType` and/or
  `columnTypeFilter` (`numerical` | `categorical` | `int` | `double` | `string`), the menu only
  offers matching columns (`buildColumnMatchFilter` → `Column.matches` / `semType`). Both attributes
  together require both; an unrecognized `columnTypeFilter` is ignored. Non-func nodes (Select
  Column utilities) are unfiltered.

### Precise run-state invalidation

* Graph edits are now **classified** (`GraphEdit`: node added/removed, connection added/removed,
  parameters changed, cleared) and invalidate **only the affected downstream cone** — previously
  *any* change (even adding a disconnected node) flipped every node to "Out of date".
  * Adding a node invalidates nothing (it isn't wired to anything yet).
  * Adding or removing a connection (data, pass-through, or order edge) invalidates its **target**
    and everything downstream; the source keeps its completed result.
  * Editing a node's parameters or inputs — property panel, function-editor dialog, column picker —
    invalidates that node and everything downstream. This previously invalidated **nothing** (the
    panel wrote values silently); title and description edits stay cosmetic and invalidate nothing.
  * Removing a node drops its state; its connections' removal events handle the downstream cone.
* Invalidation is precise across all run artifacts: node status/visuals, incoming-wire styling,
  outgoing wire labels, captured live values (`__ffFlowLive` — so single-node rerun and column-picker
  reuse stay available for still-valid nodes), and the output preview (closed only when it shows an
  invalidated node).

### Autorun

* New **ribbon toggle (bolt icon)** — faded outline when off (0.8 opacity, default font-weight),
  colored **and filled** when on (font-weight 600 renders the FA bolt as its solid variant), with a
  state-aware tooltip. When on, the flow **reruns automatically** (debounced, 1 s after the last
  edit) on any result-affecting change.
* **In-place transforms are isolated per node** (instrumented runs): every dataframe crossing into a
  function step is snapshot-cloned (`__ff_clone`) before the call; the pass-through, the live-value
  stash, and inlined `table.col(...)` args all use the snapshot. Previously an in-place function
  (descriptor calcs, AddNewColumn, …) mutated the very instance the upstream node had captured — so
  a node's preview showed columns added by *downstream* nodes, an open shell table picked via Select
  Table got modified, and an autorun slice re-run applied the transform **twice** (the boundary
  value it read had already been mutated by the previous run). Now every node's captured value is
  the state *at that node*, and re-runs are idempotent. Clean (exported) scripts are unchanged —
  they run once from scratch and keep the platform's in-place idiom.
* **Switching autorun on runs immediately**: everything without a fresh result (a never-run flow
  entirely; a half-run flow just the missing part + downstream) is scheduled at once — no need to
  make an edit first (`ExecutionController.pendingNodes` → `AutorunScheduler.kick`).
* **Clicking a node is not an edit**: building the context panel initializes its inputs, and DG
  inputs can fire `onValueChanged` during initialization — previously each selection reported a
  parameter edit (invalidating results and, with autorun on, rerunning the flow on every click).
  Every property-panel change handler is now guarded by a per-editor **value-change reporter**
  (`PropertyPanel.changeReporter`): an edit is reported only when the value actually differs from
  the last seen one (scalars compared by string form, arrays/objects by JSON). Covered by a
  user-side test (repeatedly opening the panel for value-bearing nodes, incl. OpenFile, emits zero
  `params-changed`; a real edit reports exactly once, re-entering the same value not at all).
* Reruns are **incremental**: only the invalidated slice runs, with its boundary inputs read from
  values captured by the previous run (the single-node-rerun machinery); falls back to a full run
  when nothing was captured yet. Consecutive edits coalesce into one run; a run in progress
  postpones the next one.
* Autorun is deliberately silent: no validation toasts (an invalid mid-edit graph just waits for
  more edits), no run dialog (flows whose run would prompt for script inputs are skipped), no
  selection stealing — but if the invalidation closed the output preview, the autorun re-opens it
  with fresh values once its node completes.
* **Fixed: an autorun firing while a function-editor dialog was open hijacked the dialog.** The
  dialog intercepts the global `d4-before-run-action` event to save values without executing — but
  that event fires for **every** client funccall, and the interception matched by func. An autorun
  slice re-running the same function (e.g. the AddNewColumn you were editing) mid-dialog was
  mistaken for the dialog's run action: the rerun's call got canceled and the round-trip resolved
  early with the wrong funccall — the user's OK then wrote nothing back (the "works the second
  time" race). Three-layer fix: the autorun scheduler is **held** while an editor dialog is open
  (edits accumulate; the rerun fires right after the writeback), the launcher waits for any
  in-flight run to drain before opening, and the interception **ignores events while a Flow run is
  executing** (`createFuncCallEditor`'s `ignoreEvent`). Reproduced deterministically in a test: a
  concurrent same-function call while the dialog is open must execute uncanceled and must not
  resolve the round-trip.

### Only ready nodes run

* **Runs (autorun and the manual Run/Debug) never execute an unready node** — one missing a
  required input or property — nor anything downstream of it. Previously a plot with no table, a
  Select Column / Add Table View with no table, or a function missing a required column would be
  emitted and fault at run time (or silently do nothing); now the whole downstream cone of every
  unready node is pruned from the run and the ready subgraph runs on its own.
  * **Readiness = required inputs + required properties.** Node requirements were widened beyond the
    func/viewer required *inputs* (`requiredInputs`) with a new **`requiredProps`** (panel
    properties that must be filled). Declared on the custom nodes that needed it: **Select Column**
    (`table` + `columnName`), **Select Columns** (`table` + `columnNames`), **Select Table**
    (`tableName`), **Add Table View** (`table`). `nodeMissingRequirements` unions both and drives
    the run gate **and** the on-node **"Needs input"** hint (which now flags an empty Select Table
    name, an unset column, …, not just unwired sockets).
  * **Autorun** silently drops the unready cone; if nothing is left ready, it waits for the next
    edit. **Manual Run/Debug** runs the ready subgraph and shows a one-line warning naming the
    nodes it skipped; if *nothing* is ready it says so instead of running an empty script.
  * `ExecutionController.runnableNodes()` exposes the ready set (pure — used by the gate and tests).

### Output preview: add to workspace, in-place column view

* **"Add to workspace" is now on every preview**, not just the context-panel dataframe row — a small
  `+` overlaid in the preview's top-right corner (same corner treatment as the viewer gear).
  * **DataFrame / column** → opens a fresh table view over a clone
    (`grok.shell.addTableView(t.clone())`).
  * **Viewer** → every viewer is DataFrame-backed, so it opens a table view over a clone of the
    viewer's dataframe, recreates the viewer on that clone (a **new** viewer, same look via
    `getOptions().look` → `setOptions`), and docks it to the **right** of the table view. Enabled for
    any viewer-output node.
* **Column outputs of in-place transforms now preview the whole table, scrolled to the produced
  column.** Most column-producing functions (e.g. Add New Column) mutate their input table and return
  the added column. Showing that column alone hid the context it belongs to. The instrumented run now
  detects — **by instance** (Dart-handle identity, since each wrapper access makes a fresh JS object)
  — whether the output column is one of the node's *single* input table's columns; if so it captures a
  clone of the whole table plus the column name (`__ff_col_summary` → `tableClone` / `scrollToColumn`)
  and the preview renders the table scrolled to that column. Zero or several input tables, or a
  genuinely standalone column, keep the one-column grid. Add-to-workspace works in both modes.

### Demo flows shipped as `.flow` scripts

* The two bundled demos (`Workflow Demo`, `Sequence demo`) are now also committed as first-class
  **flow scripts** under `scripts/` (`.flow` — the language the `flow` scriptHandler registers), so
  they register as runnable/reusable functions on publish, not just Start-panel templates. Each
  carries the proper annotation header (`//name` / `//language: flow` / `//tags: flow` /
  execution-ordered `//output:` lines) followed by the ffjson body — the exact `flowScriptText`
  format. A new **Flow: bundled flow scripts** test regenerates the canonical body from each
  `files/*.ffjson` and asserts the committed headers match, so the hand-committed files can't drift.

### Tutorials & how-tos updated for everything above

* The **interface tour** gained stops for the Autorun bolt and the bottom output panel; the status
  step explains the downstream-only "Out of date" marking; the Save & Open step now describes the
  **platform** save (with `.ffjson` import/export in the Flow menu).
* Stale how-tos fixed: `how-save` / `how-open` no longer claim Save downloads a `.ffjson` file.
* Five new how-tos: **rerun automatically** (autorun toggle, hands-on), **why "Out of date"**,
  **edit parameters in the function's own dialog** (full demog → AddNewColumn → Open editor ladder,
  gated on the real dialog), **reuse a saved flow** (Workflows section), and **rerun just one node**.
* The tutorials mention the output panel, autorun, and Workflows where relevant ("Load data…"
  finale, "Find the right function" group-by step, "See & reuse the generated script" save/finish
  steps). The guide content test now pins this coverage (the five question ids + the two new tour
  stops must exist).

### SetVar ⇄ Output unification

* **SetVar nodes and Output nodes now compile to the same thing.** A SetVar node also declares a
  script output: `//output: <type> <name>` (type inferred from the connected source socket, like
  Value Output's on-connect auto-typing) plus the `<name> = <value>;` assignment. An Output node
  also registers its value in the run context via `SetVar` (and, for a dataframe, under its
  runtime name), so name-based consumers (Select Table, downstream scripts) resolve it either way.
* Same in **creation scripts**: an Output node anchors the producer exactly like a SetVar —
  `T = OpenFile(...)` with no intermediate variable and no redundant `T = T` line.
* After a run, the first-output auto-preview now lands on SetVar terminals too (an imported
  creation-script flow, whose only terminals are SetVars, opens its first stored value just like
  a flow with Output nodes).
* The validator's duplicate-name check now covers the shared namespace: two SetVars, two
  outputs, or a SetVar and an output using the same variable name is an error (they would
  silently overwrite each other).

### Workflows section

* Saved flows (a `DG.Script` with language `flow`) are themselves usable as functions inside
  Flow. They no longer masquerade as Data Sources: they get their own **Workflows** section in
  the function browser — in **every** grouping mode (category / role / tags / package) — and the
  drag-out suggestion menu labels them `(Workflows)`.

### Bug fixes (this round)

* **Select Table fails fast**: when no open table or context variable matches the configured
  name, the emitted lookup now throws `Select Table: no open table or variable named "…"`
  instead of passing `null` downstream and failing far from the node that caused it (in
  instrumented runs the throw surfaces as that node's error).
* **SetVar / GetVar are always registered**: both fall to the primitive-only catalog exclusion
  and were previously registered only as a side effect of importing a creation script — so a
  saved `.ffjson` containing SetVar/GetVar nodes silently dropped them when opened in a fresh
  session. `registerAllFunctions` now force-registers them.

### Output panel rework

* The run-output preview is now a **real part of the Flow view**: a bottom pane of a vertical
  splitter (`ui.splitV`) above the status bar — resizable via the divider — instead of a
  `grok.shell.dockManager` dock with custom show/hide/lingering event handling.
* It is **not closable**: it opens on the first renderable output (clicking a completed node, or
  a run's focus node) and **minimizes to a slim header strip** via the caret at the right edge of
  the header (only the caret is clickable, so near-miss clicks aimed at the splitter divider never
  collapse the panel). The choice is remembered — while minimized, new content updates in place
  and never pops the panel back up; an explicit caret click restores it. Cleared and hidden when
  the graph changes or a new run starts (values are stale).
* The panel exists only in the real editor view — embedded hosts (the creation-script dialog)
  create the view with `{outputPanel: false}`; **Open In Editor** re-enables it.
* Removed the dock-era hacks: the `onCurrentViewChanged` close subscription (an in-view pane
  can't linger over other views) and the `onDocked` → minimap auto-minimize workaround (the
  minimap lives inside the canvas pane, which now shrinks with the splitter).
* Multiple-views hardening: `ExecutionController` no longer wipes the page-global live-value
  registry (`__ffFlowLive`) wholesale — it deletes only its own flow's node ids, so a run in one
  Flow view no longer disables single-node rerun in another. Added a standing rule to
  [CLAUDE.md](CLAUDE.md): no page-global mutable state — several Flow views can be alive at once.

### Bug fixes

* Fixed instrumented runs failing with `log is not defined` when the flow contains a **Log**
  (or Info / Warning) node: side-effect-only utilities declare no variable, but the compiler
  gave them a phantom `value` output anyway, so both the instrumented output summary
  (`__ff_summarize(log)`) and the live-value stash (`__ff_stash(…, {"value": log})`) referenced
  an undeclared identifier. `compileUtilityNode` now declares no outputs for output-less
  utilities, and the instrumented wrapper summarizes only steps that actually declare their
  variable.
* Fixed giant scrollbars around the canvas after the splitter rework: as a direct `.ui-box`
  child, the canvas container is forced to `overflow: auto !important` by core css (both the
  `div.ui-box > div.ui-div` rule and the huge `:not()`-chain fallback rule for other children).
  A higher-specificity package rule (`div.ui-box > div.ui-div.funcflow-canvas-container`)
  restores `overflow: hidden`; covered by a computed-style regression test.
* Fixed connections detaching from collapsed nodes (edges frozen at stale positions after
  moving nodes) and collapse carets going dead. Root cause: the editor↔React bridge was a
  page-level global (`window.__ff_editor`) that any newer `FlowEditor` construction rebound
  to itself (file previews, Browse entity previews, the creation-script dialog) and that any
  `destroy()` deleted (detached compile editors, e.g. running a saved flow as a function).
  The bridge now lives on each node (`FlowNode.editorBridge`, stamped by the owning editor
  on `nodecreate`), so any number of editors coexist on a page
  ([node-component.tsx](src/rete/node-component.tsx), [flow-editor.ts](src/rete/flow-editor.ts)).

### Editor

* Every way of opening a flow — a saved flow entity (editor, preview, Open from platform),
  an `.ffjson` / `.flow` file, the Import .ffjson dialog, a template, or a creation script —
  now fits the graph to the screen before it is shown. When the flow loads while the view is
  not yet attached (canvas not laid out), the fit is deferred to the canvas's first real
  layout instead of zooming against a 0×0 viewport (`FuncFlowView.fitToScreen`).

### First-class Flow entities

* Introduced a hierarchical **space picker** ([space-picker.ts](src/ui/space-picker.ts)): browse
  root spaces, drill into subspaces of any depth (lazy-loaded), and create subspaces in place.
* Reworked the **Save Flow** dialog: tooltipped name/description, an optional **Bind to space…**
  step using the space picker (default stays a plain script in the user's namespace), an advisory
  duplicate-name warning (space-scoped when a space is chosen), and Save disabled until a name
  is given.
* Ribbon: the saving section now leads with a text **Save** button; **Save** opens the Save Flow
  dialog for never-saved flows (including template-created ones that previously saved silently)
  and quick-saves bound entities. In creation-script mode Save routes to the per-table creation
  scripts dialog and the platform save/open options are hidden.
* Fixed double-click on a flow (Browse tree node or gallery card) adding the editor view twice.
* Opened/saved flows now set the `/script/<id>` URL — the address updates like it does for
  scripts, and the link routes back into the visual editor.

* Introduced the **Flow entity**: flows now save to the platform as Script entities with
  `language: flow` — shareable, searchable, space-organizable, and runnable like any function.
  The body is an annotation header (name/params derived from Input/Output nodes) followed by the
  lossless `.ffjson` document ([flow-script-format.ts](src/serialization/flow-script-format.ts));
  the header parses even on clients without the Flow package.
* Added the `flow` **script-language handler** (`flowScriptHandler`, `meta.role: scriptHandler`):
  running a flow script compiles the graph to its JS twin on a detached editor and executes it
  in-browser, copying inputs/outputs through the FuncCall
  ([flow-entity-handler.ts](src/entity/flow-entity-handler.ts)).
* Added `flowScriptEditor` / `flowScriptPreview` / `flowScriptWidget` — the entity's visual editor
  (consumed by core through the new `scriptHandler.editorFunction` seam), Browse preview, and
  context-panel pane; all thin wrappers over one `FlowEntityHandler` class.
* Save/Open now target the platform: **Save** updates the bound entity (or asks for a name and an
  optional Space), **Open from platform…** lists flow entities; `.ffjson` import/export moved to
  dedicated menu items.
* Added a `fileViewer` for the `.flow` extension.
* Fixed the editor-load race: `deserializeFlow` now registers the function catalog
  deterministically (`ensureFunctionsRegistered`), and `loadFromJson`/`loadFromDoc` await editor
  construction instead of retrying on timers.

### Custom function editors on the context panel

* Functions that ship their **own editor dialog** (an `editor:` meta, or the explicit allowlist in
  [func-editor-utils.ts](src/utils/func-editor-utils.ts) — e.g. **AddNewColumn**) now show a light-blue
  **"Open editor"** button in the **Input Parameters** pane header that opens that editor for the node
  ([func-editor-launcher.ts](src/panel/func-editor-launcher.ts)). The launcher follows the column
  picker's ladder: every **table input must be connected** (balloon otherwise); an already-computed
  upstream table is reused from the captured results; otherwise it offers to **run the flow up to
  that point** first. The editor's `FuncCall` is seeded from the connections (live captured values)
  plus the panel-edited values (column names resolved to real columns of the seeded table), and on
  close the edited values are written back into the node — connected inputs are never overridden,
  and columns come back as name strings (the panel edits names, not live columns).

### Drag-out suggestion menu reads like the toolbox

* The "Add node…" popup (drag a node's output to empty canvas) now shows each function's
  **friendly name** with its **"what it does" category** in parentheses — `Add New Column
  (Column Operations)`, `Chemical Properties  (Cheminformatics)` — instead of the raw function name
  with the role-based `(Uncategorized)` segment. Searching by the raw name still works.
* Suggestions are now **ranked by canvas context**, not just alphabetically: Value Output stays
  first; then **the science you're doing** — Cheminformatics / Bioinformatics functions when the
  drag comes out of a chem/bio node (or, for a domain-less source like OpenFile, when the canvas
  already holds chem/bio nodes); then the common table next-steps (Join, Add New Column,
  Aggregate, …); then the rest. Within a tier, an **exact type match** beats a `dynamic` catch-all
  (a column drag offers real column consumers first) and **functions already used on the canvas**
  float up (pipelines repeat their ops).

### Output preview stops jumping on repeat clicks

* Clicking the same node again no longer rebuilds the bottom **Node Output** panel — same node,
  same captured results → the existing preview DOM is kept (grids don't re-mount, scroll doesn't
  reset). The panel re-renders only when the node or its captured state actually changed (another
  node clicked, or a re-run produced new results).

### Minimap minimizes when the output preview opens

* The overview minimap and the bottom-docked **Node Output** panel share the same corner — the
  minimap now auto-minimizes to its title bar when the preview panel is first docked, so the two
  never overlap. Expand it back anytime with a click on its header.

### Semantic types detected before pickers and editors open

* The tables resolved for the **column picker** and the **function-editor dialog** now go through
  `detectSemanticTypes()` before the dialog opens, so semtype-filtered column inputs (Molecule, …)
  are populated instead of empty — a captured upstream clone may not have been through detection
  yet. Best-effort: a detection failure never blocks the dialog.

### "What produces this?" — suggestions for input drags too

* Dragging an **input** socket to empty canvas now opens the same "Add node…" popup as an output
  drag — listing every node type with a **compatible output or pass-through**, real producers
  ranked above pass-through threaders. The matching **Input node** leads (a table drag offers
  *Table Input* first — the "make this a script parameter" producer), then the science in play
  (chem/bio), then **Data Sources** + common table funcs. Picking one creates the node at the drop
  point and wires its output (real over pass-through) into the dragged input.

### Drop an input drag anywhere on a producer node

* The drop-anywhere-on-node connection shortcut now works in **both directions**. Dragging out of
  an **input** socket and dropping on another node's body connects from that node's one obvious
  output: its **real output** when compatible, else its **sole compatible pass-through** (e.g. a
  table input dropped on *Add New Column* wires to its `table →` pass-through, since its real
  output is a column). Zero or several equally-good candidates → the connection is simply aborted
  (aim at a pin) — no guessing, and no suggestion menu for the upstream direction.

### Domain sections sorted by relevance

* The **Cheminformatics** and **Bioinformatics** toolbox sections are now ordered so the flagship
  package leads — **Chem** functions come first under Cheminformatics, **Bio** under Bioinformatics —
  and within each, the "most-used" operations (descriptors, properties, fingerprints, similarity, …)
  float to the top, then the rest alphabetically. (`orderDomainSection` in
  [function-browser.ts](src/panel/function-browser.ts).)

### Toolbox trim + friendlier chem inputs

* **46 more functions removed** from the toolbox (added to the `excluded-funcs.ts` denylist) — heavy
  analyses / dialogs / panels / scripts that need their own app UI or a `FuncCall`, not a pipeline
  node (Chem substructure/R-group/MMP/ChemProp/scaffold, Bio activity-cliffs/align/HELM, the
  Chemspace/PubchemApi/ChemblApi search panels, Docking/Proteomics/Dendrogram/ClinicalCase scripts,
  `core:ColumnGridWidget`, `core:FilterToColumn`, …). The live catalog drops ~283 → **228**.
* A batch of **Chem** functions gained proper input **captions**, **descriptions**, and **`Molecule`
  semantic types** on molecule-string inputs (`getMorganFingerprints`, `getSimilarities`,
  `getDiversities`, `findSimilar`, `similarityMatrixTopMenu`, `clusterMCSTopMenu`,
  `addChemPropertiesColumns`, `addChemRisksColumns`) so their nodes read cleanly in the toolbox and
  context panel. (Source-only metadata; takes effect on the Chem package's next build.)

### Function defaults load into the node

* When a node is added, each primitive input now seeds its **declared default** —
  `prop.defaultValue ?? prop.initialValue` — instead of a bare zero value. Annotation defaults that
  arrive double-encoded (`"'inner'"`) are unquoted, and string-encoded booleans/numbers are coerced
  to their declared type so the compiler emits correct literals. Since `ui.input.forProperty`
  doesn't initialize its editor itself, every DG input in the context panel (primitives, `list`,
  `string_list`, column fields) is now explicitly initialized from the stored value via the
  `stringValue` setter (guarded — an unparseable value just leaves the editor blank).

### Node auto-summary uses the real description

* The plain-language line under a node's title (shown when you haven't written your own annotation)
  now prefers the function's **`description`** — the same text the context panel's Function pane
  shows — and, failing that, uses a deliberately-set **friendly name verbatim** instead of running
  it through the camelCase humanizer (which mangled acronyms: "InChI" showed as "In ch i"). Curated
  parameterized summaries (e.g. "Loads file demog.csv") still take precedence; the humanizer now
  only handles functions with no display metadata at all.

### `meta.includeInFlow: false` opt-out

* Any function can now opt itself out of the Flow toolbox by declaring **`meta.includeInFlow:
  false`** (in a decorator's `meta`, a `//meta.includeInFlow: false` header, or script/query
  annotations). The meta surfaces as `func.options.includeInFlow`; both boolean `false` and the
  string `'false'` are honoured, checked first in `shouldExcludeFunc`. Flow's own
  `openCreationScriptFlowDialog` (an internal dialog opener that leaked into the toolbox) is the
  first user.

### Parameter captions on slots and panel fields

* Where a function's input property declares a **`caption`** (distinct from its name — e.g. `minPts` →
  "Minimum points", `dB` → "Diameter"), that caption is now shown as the **input slot label on the
  node** and as the **field label in the context panel** ("connected only" rows and the
  `ui.input.string` column / column-list / string_list fields; `ui.input.forProperty` already did
  this). Display only — the property **name is unchanged** as the slot key / `inputValues` key /
  compiled argument. A null or empty caption falls back to the name (`getParamDisplayName`).

### `list` / `list<string>` parameters are editable

* Function inputs of type **`list`** — which is how **`list<string>`** parameters surface at runtime
  (`propertyType: 'list'`, `propertySubType: 'string'`) — were stuck as "connected only". They are
  now seeded editable and rendered through the same **`ui.input.forProperty`** pathway as the
  primitives: DG builds its native **List** input, whose value is a real JS array. The compiler
  emits the array as a JSON literal (empty → omitted so the function default applies), the creation
  script emit/import round-trips it, and only `string_list` keeps its comma-separated text field.

### Files drop where you drop them

* Dragging a **file** from the file browser (or any `DG.FileInfo` / `DG.Func` from the Datagrok
  browse tree) now drops the node **at the pointer**, like dragging a function from the toolbox —
  instead of always landing in the canvas center. Double-click still adds to the center.

### Column-output preview shows the column as a grid

* A node whose real output is a **column** now previews that column as a **one-column DataFrame
  grid** — the instrumented run captures a `DG.DataFrame.fromColumns([col.clone()])` and the docked
  preview renders it like any table (with semantic types detected). Previously it showed a small
  text sample **and** the threaded "table (modified)" passthrough grid; that passthrough table is now
  **suppressed** in the preview when a column is output (it's still captured for the column picker /
  inspect). Falls back to the text sample if no DataFrame was captured.

### Native Datagrok inputs in the context panel

* Primitive function parameters (`string`, `int`, `double`/`num`, `qnum`, `datetime`, `bool`) are now
  edited with **native Datagrok inputs** built straight from the parameter via
  `ui.input.forProperty` — instead of the bespoke textarea / number / checkbox / select controls. The
  input honours the property's declared type, numeric range, **choices** (a choice-bearing string now
  renders as a combo automatically — no special-casing), and nullability. **Column** and
  **column-list** fields use a plain `ui.input.string` with its **own native caption** (a column
  name is just a string; `forProperty` would build a table-bound column picker we can't back here);
  the column-picker icon — and, for multi-table funcs, the table chooser — are appended **inside**
  the input as trailing controls via `InputBase.addOptions()`. `string_list` is likewise a native
  `ui.input.string` (comma-separated), so every field in the panel shares one consistent look.

### Parameter descriptions & package in the node/panel

* Input and output slots now show their **parameter description** on hover — on the node's sockets
  and on the context-panel input rows — read from the function's `@grok.decorators.param({options:
  {description}})` (or a query/script `[description]`). The context panel's **Function** pane also
  shows the **Package** the function comes from, disambiguating a vague name (e.g. which
  "Descriptors"). `FuncNode` captures per-slot descriptions + the package name from the live
  `DG.Func`.

### Domain sections lead the toolbox

* **Cheminformatics** and **Bioinformatics** sections now float to the **top of the function list,
  right after the Queries pane** — the science a chemist/biologist reaches for first — ahead of the
  Viewers/Widgets panes, the built-in building blocks, and the general task categories.

### Cleaner, domain-aware function toolbox

* The **“what it does”** grouping gains two domain sections: **Cheminformatics** and
  **Bioinformatics**. Functions are routed there by source package (Chem/Chembl/Admetica/… →
  Cheminformatics; Bio/SequenceTranslator/Helm/… → Bioinformatics), which wins over the
  signature-based task category, so a scientist finds all chem/bio steps together. The node
  title-bar color follows the same domain (pink / deep-purple). **Only operations belong there:**
  a chem/bio function is placed in a domain section only if it takes a dataframe/column input; a
  pure *source* (DB fetch, generator, or query that produces a table from scalars) falls back to
  **Data Sources**, and every **`DG.DataQuery`** stays in the **Queries** pane — so the domain
  sections hold things that *do something to your table*, never queries. Package sets +
  `domainSection()`/`domainCategory()` live in `type-map.ts`.
* The catalog was **de-cluttered from ~568 to ~283 nodes** by fixing the exclusion filter and adding
  a curated denylist:
  * The exclusion check now also inspects **tags** (case-insensitively), not just the role field.
    Sketchers, cell renderers, folder viewers, semantic-type detectors and
    `Internal`/`@editors`/`Viewers`-tagged functions almost always declare their kind as a *tag* —
    the old role-only check let them leak in. **Widget-producing functions are deliberately kept**
    (they flow to the Widgets pane and can be previewed); only right-click widgets (`semantic_value`
    input) and specific dashboard-chrome widgets are dropped.
  * Two new signature rules: functions taking a **`semantic_value`** input (right-click/context
    actions) or emitting a **`tablerowfiltercall`/`colfiltercall`** (internal filter-DSL builders)
    are dropped.
  * A hand-reviewable **`nqName` denylist** (`src/rete/excluded-funcs.ts`) removes the case-by-case
    residue rules can't catch — internal engine/helper getters, internal near-duplicates of a kept
    canonical (e.g. `Chem:getInchis` over `Chem:addInchisTopMenu`, `Eda:apply*` over `Eda:train*`),
    demo/test/autocomplete helpers, core plumbing (cache drops, project/publish, raw DB-query
    builders, UI-container builders), and state-plumbing RPCs. Produced empirically (live catalog
    dump + per-package source assessment) and meant to be edited by hand.

### Rerun a single node

* Right-click a node → **“Rerun this node only”** re-executes just that node using the values its
  upstream already produced — no re-running the whole slice. The option appears only when the node's
  required inputs are all connected or filled AND every connected input has a captured value from a
  prior run. Under the hood, an instrumented run stashes each node's live outputs into a tab-global
  registry (`__ff_stash`); the single-node re-run compiles with `liveExternalInputs`, so connections
  from outside the one-node slice resolve to `_ffLive(nodeId, outputKey)` registry reads instead of
  re-running upstream. Works for function, viewer, and utility nodes (a viewer wired to a
  passthrough table re-plots from the captured modified table). The registry is cleared on a fresh
  full run and on any structural graph change, so the option is offered only when the values are valid.

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