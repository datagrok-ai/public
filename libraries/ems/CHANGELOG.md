# EMS changelog

## v.next

* Added derived columns (`table.derived`): computed in the browser from other columns, read-only, and Save waits for them
* Added filter operators an app answers itself (`table.operators`), such as a structure search
* Forms: A server-assigned key reads "assigned on save"; `immutable` columns are set on a new row only; a reference opens its row; `sections` group fields under a folding heading
* Filters: The condition builder can sit in a pane left of the list; a row still being filled in is kept and not applied
* Lists: A ⋯ Sort menu with the orders an app gives (`sorts`), opening in the first; a refresh keeps the scroll position, current row and selection
* Record pages: A pending delete shows with Restore; the ⋯ menu holds the row's own actions; Save is disabled while a field shows a problem (its tooltip says which), New while on a create page
* The context panel shows the row a user moves to and the open record (`context`)
* Child tabs: New can open the child's own create page with the parent filled in (`onNew`)
* Keys with a display format (`DG-000123`) read and search in that format everywhere
* Import: References are matched by name, derived columns are computed (and not offered as targets), the table's validators are checked, and the mapping step is reusable (`ImportMapping`)
* Fixed floats read and saved with float32 noise (`0.1` read as `0.10000000149011612`)
* Fixed Import failing on the platform when a file names referenced rows
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
