# Domain UI changelog

## v.next

* GROK-20298: The unsaved-changes gate now holds on EVERY path: `EntityListWidget.create` resolves null when the first load is refused (so a page keeps the list — and the pending batch — it already has), `DomainAppView.load()` aborts on a refusal and republishes the URL of what is actually on screen (a cancelled browser Back no longer destroys the batch), overlapping loads are generation-guarded, `DomainEntityAppView.back()` and the post-Delete close go through the prompt (a Delete discards the child editors explicitly and says so), `AppView.confirmDiscard` is re-entrant (two racing gestures share one dialog), and a `beforeunload` handler is registered exactly while the page is dirty, so a tab close or a page reload asks too
* GROK-20298: A value on a column the caller cannot write is now marked `kind: 'error'` and refuses the save on NEW rows as well as modified ones — a `DomainGridOptions.defaults` parent FK the user may not write used to be dropped from the insert payload silently, saving a child row with a null FK
* GROK-20298: `DomainFrameEditor`'s per-row caches and running change count are reset when the frame's rows are added or removed by anyone else (they are keyed by row index, which every external removal shifts); `validate()` is refused during a save like the other mutators
* GROK-20298: `EntityListWidget` search fixes: the identity columns are narrowed to declared, non-reference STRING columns as the built-in Domain View does (a `like` against a uuid column errors server-side), the term's own `%`/`_` are escaped, and the box follows the CURRENT query — a filter it cannot AND into hides it and drops the term instead of showing one that is not applied
* GROK-20298: Cap banners where they were missing: grid mode now runs the count like the gallery modes, and an entity page's child grid asks for the list's row cap instead of the server's default page size
* GROK-20298: `DomainAppViewOptions.query` seeds the page (the "my campaigns" pattern) when the URL carries no query of its own, and `confirmDiscard` is no longer among the options — the page IS the gate
* GROK-20298: A successful form save keeps the page's loaded values and version current from the server's own answer, so a cancelled reload behind it can no longer produce a phantom version conflict
* GROK-20298: One predicate decides whether an action changed the row (keyed on the platform's ribbon icons, not on captions that flip like 'Watch'/'Unwatch') — the list and the entity page used to disagree about History and Watch
* GROK-20298: An entity page's child grid is mounted in a sized box, so it fills the pane instead of rendering at the viewer's default 400x300; a slow context menu no longer opens over another item (or a detached list); list and row pages carry distinct view types; the prompt releases its dialog subscription

* GROK-20298: Added the app framework — `EntityListWidget` (the "list of entities" an app opens on: platform-rendered cards, a brief list or an editable `DomainGrid`, search over the table's identity columns, per-item commands from the handler resolved on demand, the built-in Domain View's row cap and count banner, New gated on `canInsert`) and the view hierarchy `AppView` / `DomainAppView` / `DomainEntityAppView` — a browse/CRUD app whose subclass overrides nothing: list page with the `DomainQuery` ⇄ URL round trip (plus the reserved `view=` / `entity=` params), a row page with the reflective form, a tab per child table (editable `DomainGrid` with the parent FK prefilled), history, and the row's own permission-gated actions
* GROK-20298: Added the save/discard prompt (`promptUnsavedChanges` / `confirmDiscardChanges`) — the refresh policy the editing controls deliberately do not have: every rebuild an `AppView` or an `EntityListWidget` performs asks first (save / discard / cancel), a failed save cancels the navigation, and a save in flight refuses instead of prompting; the app's status bar carries the pending-change count with Save and Discard
* GROK-20298: Added `domainHandler(table)` — the registered `DomainObjectHandler` subclass for a table, or the reflective default; `DomainGridOptions.defaults` prefills rows added through the toolbar (how a detail grid sets its parent FK)
* GROK-20298: A grid built outside the DOM no longer keeps the viewer's default 400x300 size when mounted in a page

* GROK-20298: Initial release — `DomainFrameEditor` (batch editing of domain rows with the
  `~state` / `~changes` / `~errors` service-column state, one-transaction save, standard
  version-conflict dialog, rebuild-on-refresh) and `DomainGrid` (platform-decorated editable
  grid with change highlighting, per-cell validation markers, deleted-row exclusion, and
  capability-driven read-only degradation); re-exports of the core domain surface
* GROK-20298: `onDirtyChanged` now emits the true → false transition a `refresh()` produces
* GROK-20298: Undeleting a row added in the same batch restores its `'new'` state instead of
  dropping it from both the save and the discard
* GROK-20298: A version-conflict RELOAD drops the row's pre-reload edit snapshot, so a
  following edit records the value the cell actually held
* GROK-20298: The editor is closed while a save is in flight — writes, `discard()` and
  `refresh()` are refused with a warning (and `DomainGrid` locks its editing) instead of
  being silently lost or shifting the rows the results address
* GROK-20298: `DomainGrid.decorate` keeps a resolved handler that overrides `renderGrid`
  under a type of its own, matching the platform's `isPlatformDefault` collapse rule
* GROK-20298: Changing a column the user cannot write is marked as a blocking cell error
  instead of leaving the editor permanently dirty after a "successful" save
* GROK-20298: Conflict dialogs name the row by the registry's declared name column
* GROK-20298: Per-row parse cache for the `~changes` / `~errors` columns and a running
  change count; `domain-ui-*` stylesheet
