# EMS changelog

## v.next

* Memory backend: Fixed a delete ignoring the child column's `onDelete`: `cascade` and `setnull` act as on the platform instead of refusing, and `restrict` refuses with `restrict`, not `validation`
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
