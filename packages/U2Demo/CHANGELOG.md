# U2Demo changelog

## v.next

* U2 Designer: the app registers the platform viewer tags (u2-viewer-*) and imports viewers.css; Master–detail (grid) sample
* Inputs: added an `IconInput` row to the Inputs tab; the package imports `grid.css` and `icon-input.css`
* Added `u2Record(text)` — the published function to wire a button to from the designer: it takes a literal argument and leaves a trace that outlives the balloon by appending to the demo form's Run log. Nothing published both took an argument and showed anything, so `wire a button with a literal argument` could not be completed at all: `info` and `balloons` take nothing, `test` runs the test harness, and the `u2Record` the e2e checks used was registered client-side by the check itself
* U2 Designer: the demo context's run log is `demoRunLog` rather than `log`, and every command and function that fires writes it through one place — a key named `log` on a brand-new empty form reads as the designer's own plumbing rather than as the host app's data

* U2 Designer: `demoOrders` is published, so the local-mode fixture knows it as a package function — the e2e checks that build a data-bound form now reach it through `DG.Func.find` like any other function, with no client-side registration standing in for it
* U2 Designer: the demo context is platform-backed (`platformContext`), so a button the designer's function picker wires runs the function it names
* Added `demoOrders` — a small orders table filtered by a `days` parameter, the function the u2 designer's data sources are demonstrated and tested against
* U2 Designer: the demo's `cmd:save` also appends to a bound 'Run log' text area in the form — the balloon autohides in seconds, and two acceptance passes misread the command as silent
* U2 Designer: a new app (Dev) hosting u2's read-only spec inspector over a representative form — pick a node on the canvas or in the structure tree and the context panel shows its properties, events and bindings; the ribbon opens another spec, copies the current one, and toggles Design/Run
* Input convergence: where the pane cannot hold both editor columns the grid scrolls sideways on its own instead of the page — the prose, the notes and the status line stay where they are, and the page says so and what to close to get both columns back
* Input convergence: the intro points at the Files & columns tab for the pickers that need an open table

* Input convergence: the grid now carries the whole input matrix — a property record per editor the platform routes to (radio, textarea, password, color, font, image, molecule, bounded int/float, slider, bigint, qnum, switch, datetime, multi-choice, list, map, file), each rendered by both generators
* Input convergence: a leading `…` per row opens that property's own metadata in the context panel — a u2 PropertyEditor whose fields follow the property's type — and every edit reconstructs both editors of the whole grid, resetting the row's value when its type changes
* Input convergence: the page names the platform's hover chrome (each side reveals its own, so the two cells are compared one at a time) and carries inline notes on the rows that degrade here — the molecule row without Chem, the file row in local mode, the empty bigint row, and the Tags row left out because the platform's own editor for it cannot add a tag, each note where its row sits
* Input convergence: `…` opens the context panel when it is closed, names the property it is editing (and keeps naming it through a rename), and highlights the row being edited; a rename onto an empty or taken name now says so on the field instead of being dropped in silence
* Input convergence: the image preview and the map editor stay inside their column instead of bleeding over the next one and pushing a page scrollbar, so their labels stay beside them
* Object form: SAVE and Delete report the failure in words on the result line, and the tab says up front that persisting needs a server
* The app reopens on the tab it was left on, so a reload for the dual-run toggles comes back to the convergence page
* Input convergence: the page says how the two sides' display formats differ and what local mode cannot do (the app card's error chip, no deep link back)
* Input convergence: the A/B is one grid with a row per property, so each Dart editor and its u2 counterpart share a baseline instead of drifting apart down two columns
* Input convergence: the sync counters now count user edits that crossed the bridge — a fresh load reads 0 · 0 and a write-through echo no longer counts as an edit from the other side
* Input convergence: toggling a dual-run type shows what will happen and offers a reload, and a type that registered at startup is marked (active)
* Input convergence: the molecule row says what it needs instead of disappearing where Chem is absent, and the page prose renders its code spans instead of showing backticks
* Imports every u2 stylesheet the library ships (the DateTime editor was rendering UA chrome for want of css/date.css)
* The Molecules tab asks for the Chem package instead of building bridged widgets without it — five 'Unable to get project asset drawMolecule' errors on app open in local mode

* Added the Input convergence tab: one DG.Property set rendered twice over one target object — the platform's InputBase.forProperty on the left, u2's propertyForm on the right — cross-synced with value-comparison echo suppression
* Added the Files & columns tab: the platform-bound u2 inputs — fileInput, filesInput, rsaInput, the rich columnInput, tablesInput, columnsInput, columnsMapInput and aggregatedColumnsInput — each with a value readout
* Added opt-in dual run: per-type toggles register u2 inputs as platform `role: valueEditor` funcs (int/double/bool/datetime/string, persisted in the `u2.valueEditors` local-storage key), so unmodified Dart forms and func-param dialogs resolve u2 editors with the property's own metadata
* Input convergence: the Dart → u2 sync is a plain callback again, and the value editors use the library's own `assumeWritable` — both workarounds for library bugs the bridge now fixes

* Reports Browser: rich list rows per the list-item-rendering recipe — hover-revealed copy/open-ticket actions with a right-click menu of all actions, compact timestamps with full date-time on hover, and adaptive hiding of the reporter as the pane narrows

* Both apps open from the Browse tree on single click (app functions return the view) and ride the shell chrome via u2's appView: search and refresh in the view ribbon, live "N of M reports" in the per-view status bar

* Introduced the U2 Demo app: a tabbed tour of @datagrok-libraries/u2 (inputs, combobox/async, virtual lists and trees, layout, popups, form and property grid, platform bridge)
* Added the Entities tab: user/group pickers over dapiSource + ObjectHandler rendering, entity chips with handler tooltips and current-object wiring
* Added the Reports Browser app: a master-detail browser over grok.dapi.reports on AsyncPager/dapiPager — virtualized incremental list, smart filter, and a detail pane with reporter/assignee chips
* Added the Spaces & dashboards tab: space and project pickers with a live ObjectHandler renderCard preview (entityCard), chip, and Open action
* Added the Molecules tab: the platform sketcher bridged as a u2 input (moleculeInput), semType-driven editor selection in propertyForm via the new Editors registry, and structure chips/typeahead rows via moleculeRenderer
