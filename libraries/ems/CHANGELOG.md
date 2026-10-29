# EMS changelog

## v.next

* Create Domain Schema: Recognizes what text columns hold from a sample of their values (a SMILES column is drafted as a Molecule), so bound tables draw structures with nothing declared
* Lists: A table opened at `/domains/<schema>/<table>` shows its rows as cards with the structure (any column a display rule draws) beside the name and two text lines, or as a table of columns where no column has a picture, with the condition builder beside them
* Tables: A column a display rule draws (a structure) is drawn in its cells, the rows tall enough for it; added the `view` app option to open the list page as a table
* History: Fixed "Show older" turning off with no reason while it loads; it now says "Loading older entries…"
* u2-domain-pick: Fixed a pick bound to a read-only path rendering a spec error instead of a read-only input
* u2 v1: Moved onto u2's v1 API (plain `value`/`peek()`, `build`, owners as options, captions as `text`, `disabled` reasons, entry-point imports). Save and Discard say why they are off ("No changes to save", "Nothing to discard", "Wait for the batch being saved to finish"), as do the authoring dialogs' and the manifest editor's controls that are off, and New says "Loading…" until its table is loaded; `DomainAction` takes `disabled` (the reason) instead of `enabled`, and the `u2-domain-pick` tag's `enabled` prop is `disabled`
* Forms: Fixed text a number, date or color field cannot read ("abc") not counting as an unsaved change: Discard drops it and the gates ask before leaving
* Properties: Added an editor for users' own columns on a plugin's extensible tables, opened from Edit properties… on the list's ribbon (`domains.editProperties`); the plugin's columns stay read-only, and a save that deletes values asks first
* Forms: Fixed a field left invalid going stale in the app: the breadcrumb and the entity page's status name the record as its fields read (an emptied title is "New issue" again, not the last valid one; a saved record's emptied name leaves the name it has, never a blank crumb), and a field over a cell already changed is no second unsaved change ("2 unsaved changes" for one edited title). `holdUnwritten` takes what the fields hold for their row (`UnwrittenEdits`); added `DomainSource.current(row)`
* Import: Fixed writing the key back failing ("Column named 'Sku' already exists") when the source table has a column of that name in another case (`sku`): column names ignore case, so the key goes into a column of its own beside it, and the source's column keeps its values
* Import: Fixed the key written back, and the report's key column, holding each row's caption ("Row 1", its name) instead of its key: a one-column key is its own value, as the column's format spells it (`DG-000060`)
* Import: Fixed closing the dialog while the report reads the keys back building its rows table under the closed dialog (a u2 usage error); the key is still written back
* Apps: Fixed browser Back or Forward inside an app's view (`domains.views`) opening a second view of the table: the app function the platform asks again answers no view, and the open one takes the address; with unsaved changes the gate asks once
* Forms, Import: A required field is starred through the input's own `required`, so it is `aria-required` too (the asterisk was a class set by hand)
* Apps: Fixed "Wait for the batch being saved to finish." showing when a record made with New is saved, then deleted from its own page and saved again: the saved delete no longer reads as a discarded draft, and the app goes back to the list once
* Apps: An app's own parts are `DomainApp` options instead of a subclass — `ribbon(base, app)` and `menu(base, app)` extend the ribbon and the Actions menu the app builds, `setup(app, self)` wires the rest under the app's owner, beside `shortcuts`; removed the `app` class option of `DomainTable.app()`
* Added `domains.configure(address, fn, {id})`: an app's rules apply to every handle over the table
* Bulk edit: Redesigned the dialog as one form with a field picker
* Bulk edit: Refusals and the table's validators show in the dialog
* Import: The Preview button names what lands and says why it is held
* Import: Refused cells are tinted and a skipped row says why
* Child tabs: Grids fit two-line captions with units
* Child tabs: A tab that fails to open shows its error
* Data tables: Fixed a table stopping at its first page
* Filters: An invalid query shows its reason under the box
* Unsaved changes: The dialog reads CANCEL · DISCARD · SAVE
* Trash: Restore uses the undo icon
* Forms: Fixed touching one field of a new record showing every required field's error: a field shows its message once the user edited or left it, or a save validated the form (the form re-reading its row no longer counts)
* Save: Off while the form or the writer refuses the batch, the reason as its tooltip ("Title is required"), and Tab from the last field reaches it; `DomainSource.guard(check, problem)` takes the verdict live
* DomainApp: A read-only page — the trash until a restore is staged — shows no Save or Discard
* Forms: Only a valid edit reaches the row; an invalid one stays in its field and still counts as an unsaved change (`DomainSource.holdUnwritten`): every gate asks before dropping it, Discard drops it, Save stays off saying why — the unsaved-changes prompt's SAVE too. The writer's verdict on the row's value no longer blocks a new value, and a constraint is re-judged on what the fields hold
* Bulk edit: OK says why it is off in its tooltip
* DomainSource: Added `hosted`: the app answers a host filter condition (a structure search) with the ids of its rows, sent as an id list before every read
* DomainApp: Added `filters: {mode: 'builder', placement: 'left'}`: the condition builder in a pane beside the list
* Import: Added `table.importer`: an app derives columns, keeps rows out (already there, repeated in the file, unreadable), maps its own columns and shapes the preview's and report's columns and messages; the report lists every row with its outcome and key, written back to the table on request
* Import: An import is held to the table's validators, as a form is, and a unique refusal reads "Code is already used", not the database's message
* Added `domains.views`: a view per table under one app path, with deep links and tree items; a row opened from anywhere lands in its table's view
* DomainApp: Added `sort`: the list's initial order
* Notifications: Fixed refusal balloons piling up and outliving the fix: a source shows one "Cannot save" balloon, replaced by the next refusal and closed by the edit that settles it, so it no longer sits beside "Saved"
* Forms: Fixed a Save refused by a table validator leaving a new record's other problems hidden: every field's message shows, and the first field to fix takes the focus, as when the form itself refuses
* Filters: Fixed the builder wiping a condition still being written: it sends only its complete conditions and keeps the rest
* History: Fixed json and list columns showing as changed on every update (`[object Object] → [object Object]`); json values read as JSON, and restores show "restored"
* Forms: Showed the column defaults on a new row's form
* Memory backend: Filled the column defaults into a draft, and applied the default on an insert that leaves the column out, as the server does
* Memory backend: Fixed derived singular and plural names reading the friendly name and replacing `_`: like the platform, they come from the table name (`order_details` → `order_detail`), and a declared friendly name is answered as `info.friendlyName`
* Memory backend: Fixed derived table names: a table named in the plural (`categories`) read `categoriess`; singular and plural now follow the platform's rules (`category`/`categories`, `box`/`boxes`, `study`/`studies`)
* Memory backend: Fixed `updateWhere` rewriting the same rows on every run; like the server, it skips rows already holding the values
* History: Loaded a row's trail 50 entries at a time, with "Show older" for the rest; added `{limit, before}` to `DomainTableLike.audit` and `seq` to `AuditEntryLike`
* [public#4320](https://github.com/datagrok-ai/public/issues/4320): Import: Fixed an empty numeric CSV cell being posted as the missing-value marker (-2147483648, 2.68e-34) instead of null
* Styles: `src/dg/styles.ts` imports EMS's own sheets only; u2's stylesheet installs itself once per page
* Tests: Fixed the spec-install test on Windows (path separators, CRLF recipe fences)
* Added `installEms(registry?)` and `emsPlugin`: one install registers `u2-domain-source`, every `u2-domain-*` control, `u2-manifest-editor` and the spec session; `registerEmsComponents` is a deprecated alias; `emsSourcesPlugin` is the platform-free install (the source tag and the spec session), which `registerDomainSources` now installs, returning the uninstall
* Breaking: Importing ems no longer sets a u2 global; the spec session comes from the installed plugin, so a host calls `installEms()`
* Registrations: The construction options of the domain controls and `u2-domain-source` are structural: a bind reads them once and warns, instead of "has no live setter"
* Manifest: Added the kind catalog (`kinds`, `kindsVersion` 1) for the `data-u2` values EMS stamps (domain controls, session buttons, authoring dialogs)
* Manifest panel: Parts (title, form, sections) are keyed and kept while they show what the model holds, so focus and the caret survive an edit; a replaced part hands focus to the element with the same keys
* Bulk edit: The dialog stays open and busy while the write runs; a failed write keeps it open with what was checked and typed
* DomainSource: A disposed source no longer starts a refresh, and the app opens no page once disposed; `afterSave` no longer overwrites a newer refresh's total with a late count
* Docs: Fixed the spec-app recipe (the source sits on the components tray, the spec mounts after `installEms()`); the README and `docs/dg.md` describe `installEms()`; a test renders every package `*.spec.json` and every recipe spec on the memory backend and fails on any placeholder
* Packaging: Replaced `.npmignore` with a `files` allow-list (226 → 98 files), and the build starts from an empty `dist`
* Testing: `pnpm test` fails with a clear message when `u2/dist` is missing or older than u2's sources, and works on Windows
* Split out of `@datagrok-libraries/u2`
* `domainFilterSchema(address)` (was u2's `FilterSchemas.forDomainTable`) is exported from `src/dg`
* Gallery and browser checks: `gallery/index.html`, `npm run e2e:local` (u2's runner with `--checks`), `npm run e2e:stand-domain`
