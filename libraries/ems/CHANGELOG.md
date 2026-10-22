# EMS changelog

## v.next

* Forms: A draft's server-assigned key reads "assigned on save", a failed derivation shows on its input field as it ends, and a create page starts in its topmost form without scrolling
* Forms: `immutable` columns (also a `domains.app` option) are set on a new row only; a readonly reference is a link to its row
* Pickers: An open control beside a picked row (`DomainPick.open`, `DomainPick.opener`)
* Filters: The table's own columns come first, no empty value is offered, and a wrong builder row keeps the applied filter
* Filters: Fixed a host operator removing the core operators from other tables' columns of its semType
* Import: Fixed reference names failing the whole import on the platform: the lookup now sends what the server accepts
* Import: The mapping step is reusable as `ImportMapping`; an unmatched reference is reported once
* Fixed a batch saved with an inline new row showing a warning for each discard attempted while it saved
* Memory backend: Answers the system columns on every row and refuses a request naming one, and refuses `in` with a list, as the server does
* App: Lists open in creation order (`sort`, `DomainApp.defaultSort`), and the ⋯ menu's Sort switches between the table's orders (`sorts`, `DomainApp.sortChoices`)
* App: Delete on a record page shows the pending deletion with Restore, and the saved deletion goes back to the list; the page no longer offers Open for itself, and a missing row no longer shows its id
* App: The left filters pane folds to a strip while it holds no condition, so narrow screens keep the list
* App: The row a user moves to in the list is shown in the platform's context panel (`context` maps it)
* Children: `onNew` lets a child tab's New open the child table's create page with the parent filled in (`DomainApp.create(values)`, `DomainApp.createAt`)
* Captions: A reference to a numbered row reads `DGB-000045` in child grids too, and a numbered key column reads in its display format there
* Search: A bare number (`12`) finds the row whose formatted key it is (`DG-000012`)
* Status: Several tables' changes are named by table — "1 batch and 1 compound saved", "2 unsaved changes: 1 batch and 1 compound"; a refusal over a computed column names the input it is computed from
* Fixed a count from an outdated reload overwriting a newer one, and a structure filter being answered again on every new order or search
* Memory backend: A row without a name column is named by its one-column business key, as on the platform
* Added `table.derived`: columns the browser computes from other columns of a row; Save waits for them, a failure (or no answer in 60 s) refuses it, forms and data tables show them read-only
* Added `table.operators`: filter operators the app answers for a column the server cannot search that way (a structure search) — offered in the builder on their semType, sent to the server as the ids they match
* App: Added `filters: {mode: 'builder', placement: 'left'}` (the condition builder in a pane on the list's left), `itemHeight`, and `inline`: a create page may point a reference at a new row of its table, filled in on the same page and saved in one transaction
* App: An entity page's ⋯ menu holds the row's Delete and actions; Import runs `importRows()`, which a subclass may override; a table with no searchable column has no search box
* Search: A term spelled like the caption column's formatted number (`DG-000123`) finds that row
* Fixed floats read and written with float32 noise: a stored `0.1` read `0.10000000149011612`, a typed `1.23` and a CSV's `0.73` were saved as `1.2300000190734863` and `0.7300000190734863`
* Captions: A row named by a formatted number reads `DG-000123` in lists, cards, data tables, pickers and references, an unresolvable reference reads "(unknown)", and a row keyed by a reference reads by that row's caption; a picker finds a typed id by its number and other text by the target's searchable columns
* Data tables: Numbers show in their display format
* Children: `children.tables` takes the app's table handles, so their validators and renderers apply in the tabs
* Buttons: New is disabled while a save runs
* Status: Counts and messages use the table's declared names
* Lists: A poll that sees the table move, or a save elsewhere, re-reads the loaded rows in place: the scroll position, the current row and the selection stay
* Filters: The builder applies the conditions that are complete; a row still being written no longer clears the builder
* Import: Reference columns are matched to rows by name (`DGB-000045` → that batch); a name that matches nothing is reported on its row
* Import: The wizard runs the table's derivations over the mapped rows, and lists a row one fails on with the server's verdicts
* Memory backend: A column's `format` reaches its property, as on the platform
* Memory backend: Fixed derived singular and plural names reading the friendly name and replacing `_`: like the platform, they come from the table name (`order_details` → `order_detail`), and a declared friendly name is answered as `info.friendlyName`
* Memory backend: Fixed derived table names: a table named in the plural (`categories`) read `categoriess`; singular and plural now follow the platform's rules (`category`/`categories`, `box`/`boxes`, `study`/`studies`)
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
