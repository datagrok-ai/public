# EMS changelog

## v.next

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
