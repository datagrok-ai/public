# Domain UI changelog

## v.next

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
