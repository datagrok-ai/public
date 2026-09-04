# API Tests changelog

## 1.10.3 (WIP)

GROK-20298: `dapi/domains.ts` — `grit.issue.number` is auto-numbered since Grit 2.1.0, so its row property is expected nullable (the server fills it)
GROK-20298: New `dapi/domain-cross-schema-refs.ts` pins the soft cross-schema ref through PlatesFixture's `plate.source_query` (`ref: Core.queries`): an insert with a visible query id is accepted and an unknown id is refused with the per-row `fk` error (409), `source_query.friendlyName` filters and expands through the ref, a categories facet on the ref column carries the query's name, and deleting the query leaves the row readable with the dangling id while expand and filter find nothing.
GROK-20799: New `dapi/entity-properties.ts` pins the `grok.meta` discovery surface — catalog order and labels, the filterable subset, `refType`/`relationKind`, the domain-table branch, the uniform null miss, `coreLocationOf`, no secret in any catalog, the drift probe (every filterable User property filters through the Core query route), the no-oracle refusal of a secret column, the session→user hop agreeing with dinq, and the structured 403 on a Core write.
GROK-20799: `dapi/domain-lifecycle.ts` pins the server-composed `travelableRelations` and `securingTable` of `DomainTableClient.capabilities()` (declared relation travelable for admin and not for an ungranted user; a master-mode junction secured by its delegate target)

GROK-20753: `functions/param-eval.ts` adds the W4 pins — `Functions: ParamValidators` (client-registered sync validators run via `FuncCall.evalParamValidators` against the call's current value: passing results omitted, `false` → the "didn't pass the <friendlyName> check" message with `isError: true`, a string result becomes the message, declaration-order results, unknown-name and async-registered rejections with the sync-required message) and a `ScriptSync` comparison-expression pin over a variables map.

GROK-20753: New `utils/string-utils.ts` pins `DG.StringUtils.levenshteinDistance`/`jaroWinklerDistance` (pure-function values).

GROK-20753: New `widgets/pickers.ts` pins the table-picker dialogs (`ui.pickTableFromFiles`/`ui.pickTableFromQuery` open their platform dialogs and resolve null on cancel) and the `ColumnGrid` popup surface (`Widgets: ColumnGrid`) — construction, `onCurrentRowChanged`/`currentColumn`, the search-box filter, `close()` detaching, and the wrapped-`DG.Column` filter-callback contract on both the `popup` option and the `filter` setter.

GROK-20753: `functions/param-eval.ts` adds the W3 table-param pins (`Functions: TableParams`) — dataframe marshaling into a FuncCall (dart identity kept, same-reference writes collapse to one `onChanged`), the resolver readback shape (a string into a column param reads back as a `FuncCall`), the `ColumnList` surface of a column-array write (`names()`/`toList()`/`length`), the `grok_Property_Get`/`grok_Property_Set` `parentTableParamName` door (implicit annotation link with null `options['table']`, explicit `{table: df}`, set round-trip on a registered func), and `columnTypeFilter` derived from `{type: ...}`.

GROK-20753: `functions/param-eval.ts` adds the W2 pins — the tags∪options union on viewer properties (`.is-legend-property` on `colorColumnName`, writes into the merged copy do not survive a fresh `getProperties()` read), `evalParamChoices` propagate lookup over a client-registered DataFrame provider (`values` key→key, own column excluded from `lookup`), and static list-literal `choices` answering through the evaluator.

Functions: New `functions/param-eval.ts` covers the FuncCall parameter-source evaluators (`evalParamChoices` with `dependsOn` and the null `lookup`, `evalParamSuggestions` receiving the typed text, `evalParamDefault` including rejection on a broken command), `scriptSync` with a variables map (fresh-context isolation) plus its one-argument back-compat, and `Property.options` write-through for FuncParam-backed properties.

Tests: Removed the `grok test` Node pass — the suite runs in the browser only. `package-test.ts` no longer exports `testNode()` or takes an `excludeNodeTests` input, which are the two things `grok test` probes to decide whether to run tests headless, and every `{node: true/false}` annotation is gone. Two workarounds that existed only so the test bundle would evaluate under that pass went with it: `benchmarks.ts` closes its table unconditionally, and `js-viewer.ts` extends `DG.JsViewer` directly instead of a dummy base.

Tests: The standalone Node runner (`package-test-node.ts`) is unaffected and still drives `grok stresstest` and the nightly Stress-Tests job. It selects by `{stressTest: true}`, not by `{node: true}`, so it keeps the `typeof process` self-skips that let browser-bound tests opt out under Node — including four stress-marked ones whose loss would have dropped the sweep off its 100% baseline.

GROK-20774: `Property: Accessors` covers custom `get`/`set` on `DG.Property.js` (accessors win over the default field closures, options metadata still applied, the Dart options push does not clobber them) and a custom `set` on `Widget.addProperty` (the write lands through the setter and `onPropertyChanged` still fires).

AI viewers: Fixed the two legend-visibility tests that had never passed. `legendPresent` asked whether a `.d4-legend` element exists, but hiding a legend that has already been shown collapses it to zero width instead of detaching it — so every "legend is now hidden" assertion after a toggle was doomed, while the same assertion at creation time passed because the element had never been built. It now measures rendered width, and lives in `helpers.ts` so BarChart and BoxPlot share one definition.

Node runner: The stress suite now also loads a sibling package's node tests - `extraTestPackages` merges DBTests' registry into the one the runner filters, so the sweep covers read-only Postgres queries through grok_connect alongside the platform API. Each package resolves its own copy of the test library, so the merge is by registry object rather than by import.

Node runner: Fixed the `datagrok-api/{dg,grok,ui}` aliasing for the ESM graph — `Module._resolveFilename` returned a `dg-runtime:` sentinel that is not a real path, and whichever resolver saw it (tsx's or Node's) read it as a directory import and threw. 25 of the 49 eligible test files (every `dataframe/`, `functions/`, `bitset/`, `stats/`, `property/`, `valuematcher/` file and four `dapi/domain-*` ones) never loaded and their tests never ran. `bindRuntimeGlobals()` now generates a real module per global once `startDatagrok()` has produced them and resolves the specifiers to that, so the alias no longer depends on which resolver wins.

GROK-20752: `ProgressIndicator` covers `onLogUpdated` delivering plain `{level, message, flag, params, time}` objects — a debug-mode call of an ad-hoc JS script, asserting the client `CALL DURATION` start event's types and params.

Build: Fixed the local build — `@datagrok-libraries/domain-ui` now takes `datagrok-api` from `../../js-api` instead of the registry, so its types are the same ones ApiTests compiles against. `build-js-api-tests-local` builds js-api and domain-ui only; the `link-*` scripts are gone.

GROK-20298: `Dapi: domain frame editor` covers the move of the editor into the platform — `DomainFrameEditor.attachTo(frame, schema, table)` adopting a frame its HOST owns (the Dart Domain View's entry point), with `DG.DomainFrameEditor` asserted to be the very class `@datagrok-libraries/domain-ui` exports.

GROK-20298: `Dapi: domain frame editor` covers the ref-column gate lift — a writable `ref` column is editable in `DomainGrid` (the in-place picker's precondition), a cleared ref reaches the update op as an explicit `null`, a required ref clear blocks the save, and a pick back to the original drops the pending change.

GROK-20298: Added `Dapi: domain relations` — the many-to-many surface over the new `tags` relation of `apitests.item` (junction `item_tag`, target `tag`): insert-with-links and the deduped `{id, name}` expand, the `queryDf` chips column plus its `'~tags.id'` companion, relation filter paths (`tags.name`, `tags.id` ANY-of, the empty selection) and the owner-counting `tags.id` categories facet with display names, set-replace `update` semantics (replace, no-op re-send, absent key, `[]` clear, rolled-back missing target), the create-and-link transaction with an element-wise `$ref`, and the batch refusal.

GROK-20298: `Dapi: domain relations` covers the uncheck gesture — `'<relation>.id' != [id, ...]` excludes the owners linked to any of them (a multi-linked owner does NOT survive unchecking one of its links), and an empty exclusion list constrains nothing.

GROK-20298: `JS: domain handlers` covers the render CONTEXT the platform threads into a registered handler — the Dart dispatch (`forEntity` -> `renderMarkupFor`) hands `renderMarkup(x, context)` the `{relation, table, input}` object the many-to-many chips pass.

GROK-20298: `JS: domain handlers` covers the new `ObjectHandler.renderListItem` — the base default (the caption in a div) and the delegating round trip through the platform meta — and `ObjectHandler.renderInput` — the null base default and the `UserMeta` user selector reached through `forEntity`.

GROK-20298: Widened the domain UI suites for the WO-9F fix batch — `Dapi: domain app framework` now proves that `options.query` seeds a page whose URL it publishes, that a CANCELLED prompt leaves the very same list, editor and pending batch in place (the browser-Back path), that two racing gestures share ONE prompt, and that the entity page's Back is gated on everything the page tracks; `Dapi: domain frame editor` adds column security refusing a save whose payload would silently drop the value (a modified cell AND a new row's prefilled FK), an out-of-band row removal not desyncing the per-row state, a running-counter-vs-recount assertion at every state transition, and the collapse branch of `DomainGrid.decorate` (a foreign-typed handler that overrides nothing must fall through to the platform). The decorate test registers its handler and seeds inside the try, the pause seam releases in a `finally`, and the editable-fixture check is an assertion instead of a self-skip. `DG.DomainObjectHandler.rowFrom` garbage rejection is covered in `JS: domain handlers`.

GROK-20298: Added `Dapi: domain app framework` — `@datagrok-libraries/domain-ui`'s `EntityListWidget` and the AppView hierarchy: list rendering and server-side search over the table's identity columns, the capability-gated New button, the three render modes (the grid one exposing its editor), a one-expression `DomainAppView` whose query round-trips through the URL, a deep link restoring the query and the reserved view mode, the row page (header, FK-inverted child tab, a form edit read back from the server), the unsaved-changes prompt driven through its REAL dialog for all three outcomes (cancel keeps the changes, discard drops them without writing, save writes them and then rebuilds), and permission degradation under a restricted session. Also adds the `Domain Items` app — a complete browse/CRUD app over `apitests.item` in one expression, the vehicle for the framework's live verification.

GROK-20298: Widened `Dapi: domain frame editor` for the WO-7F fix batch — the full `onDirtyChanged` sequence across edit/refresh/edit/discard (the transition a save-or-discard prompt binds to), undelete of a row added in the same batch reaching the server as an insert, a version-conflict RELOAD invalidating the pre-reload edit snapshot (the next in-grid edit reverts to the reloaded value), the whole in-flight-save window (writes, `addRow`, `discard()` and `refresh()` refused, grid editing locked, everything released and landed afterwards) driven through a paused-transaction seam, and `DomainGrid.decorate` keeping an overriding handler that claims the table under a type of its own. The one-transaction assertion now demands a non-null `tx_id` on every batch audit row (it could pass vacuously on an all-null set), the restricted-user test fails loudly instead of self-skipping green and revokes its grant in a `finally`, and the suite documents that none of it is pure — every test needs a live server.

GROK-20298: Added `Dapi: domain frame editor` — `@datagrok-libraries/domain-ui`: the `~state`/`~changes`/`~errors` service columns (attached, tagged, and provably absent from `toCsv()`/`toByteArray()` of a DIRTY frame), the state transitions (edit/manual-revert/`revertCell`/add/delete/undelete/`discard`, with an invalid cell staying pending), the transaction-op builder (changed columns only + `expectedVersion`, no op for an add-then-delete row), the live loop (2 edits + insert + delete saved as ONE transaction — asserted through shared `tx_id` in the audit log — with in-frame version bumps and a value-for-value readback), server validation errors landing on the offending cell, all three version-conflict outcomes (reload / overwrite / dismissed), the rebuild-on-refresh contract, cooperative filtering keeping deleted rows hidden across a recompute, server-parity `validateCellValue` messages, and `DomainGrid` decoration + capability gating including read-only degradation under a genuinely restricted user. `withRestrictedUser` is now exported from `domain-lifecycle.ts` so the new suite reuses the one impersonation harness; the package bundles the library and re-exports it as `apitests.domainUi` so library-only surface can be probed from a browser.

GROK-20298: `withRestrictedUser` now runs the signup itself inside the try and resolves the user by login in the finally, so a tokenless or failed signup can no longer leave an unblocked test user behind; the Domain View `refresh()` assertion changes the row on the server first, so it can only pass on a real re-query.

GROK-20298: Widened `Dapi: domain query state` — the URL round trip now carries every list parameter (columns, joins, aggregations, groupBy, offset) and re-checks it through `toParams()`, a lone smart-filter element is asserted to reach the REST spec verbatim, and the malformed-input probe covers the new negative-limit and sub-group shape rejections. The creation-script test forces `grok.shell.settings.dataHistory` on and restores the profile's value, so it no longer depends on the session; row inserts moved inside the try that deletes them.

GROK-20298: Hardened the domain suites — the Domain View test now asserts that the reported `query` carries the invisible `permanentFilter` and the search string, that `refresh()` re-lists, and that assigning `meta` fails with the named registry error; both restricted-user tests run through one `withRestrictedUser` helper (the impersonation harness, cookie mechanics and blocking cleanup in one place) so a setup failure can no longer leave an unblocked test user behind; row inserts moved inside the try that deletes them.

GROK-20298: Added `Dapi: domain query state` — `DG.DomainQuery`: lossless URL round trip with the platform's list-element binding (reserved `view=`/`entity=` ignored, gaps closed), a live Domain View's `query` reconstructed and re-run to the same subset, `run()` matching `queryDf(toSpec())` value-for-value while recording a `DomainQuery` creation script, `fromBuilder` (AND conjuncts split per element) selecting the same rows, and client-side failures (plain `Error`s, not the typed `DomainError` family) for malformed element indices, non-integer limits, unparseable filter elements, and specs that need the function's own grammar parser.

GROK-20298: Added `DG.DomainView` and dialog-opener coverage to `JS: domain handlers` — an embedded Domain View lists exactly the rows its `permanentFilter` allows and reports its `query`, `showFilters()` docks the filter panel, and the platform dialogs are driven for real (create cancelled -> false, edit saved -> true with a server-side version bump, picker cancelled -> null, conflict dialog resolving reload/overwrite/dismiss, delete confirmation cancelled -> false), plus the audit and grants panes and `openRow` navigation through the platform's routing.

GROK-20298: Added JS-side renderGrid fall-through coverage (a plain, non-overriding handler decorates a queryDf grid identically to the platform meta) and a `defaultRowVisibility: "none"` fixture table (`apitests.hidden_item`, schema 1.2.0) with a restricted-session probe proving a table grant does NOT reach its rows — the server semantics `capabilities()`' reaches-rows rule mirrors. Ribbon-action assertions follow the new gating (Open always, Share... only with the row's Share permission, nothing at all for an unsaved row).

GROK-20298: Added `DomainObjectHandler` coverage to `JS: domain handlers` — an override-nothing handler renders DOM-identically to the platform meta (card/markup/tooltip/caption) and `getById` gives JS its first `DG.DomainRow` acquisition path; registering one keeps the Dart meta (owner of the CRUD commands) while winning `forEntity`; reflective properties/detail tabs/editor (inputs = writable columns); `DomainRow.permissions()` and the ribbon actions flip on a real Delete grant round-trip, and an unsaved row holds no permissions; malformed and unknown table addresses fail with clear/typed errors.

GROK-20298: Added `Dapi: domain registry` (rowProperties constraints round-trip incl. grit.issue choices/min/nullable, tableInfo childTables FK inversion, resolveNames business-key identity + null for unresolvable ids, typed unknown-table rejections) and `Dapi: domain capabilities` (admin full capabilities + cache invalidation; canInsert/canEdit flip on a real Edit grant round-trip probed under a throwaway restricted user's session token, cleaned up in finally).

GROK-20298: Added renderGrid coverage to `JS: domain handlers` — the base no-op is sentinel-marked `isPlatformDefault`, a registered JS handler wins the dispatch and decorates a `queryDf` frame (no `~item` assumption), and the per-table Dart meta's `renderGrid` decorates a JS-created grid (system column hidden, ref caption stamped); the suite forces the per-table meta registration so it passes under `grok test`'s unflagged client profile.

GROK-20605: Coverage sweep — `save()` business-key-duplicate regression probe (no fabricated version), typed `opIndex` on transaction rollback, `fetchFields([], fields)` zero-column shape, builder `.skip/.select` combination.

GROK-20604: Added optimistic-concurrency coverage to `Dapi: domain parity` — `save` insert/update round-trip with a typed stale-save conflict, `retryOnVersionConflict` convergence under a deterministic forced conflict, `updateWithRetry` (fresh-row mutate, null-skip, typed not-found).

GROK-20603: Added query-builder coverage to `Dapi: domain parity` — thenable chain equals the spec form, equality-map/3-arg `where`, `.df()/.first()/.count()/.exists()` terminals, apostrophe-containing values through bound conditions (`cond`/`or` helpers), raw-string + condition mixing rejected client-side; the transaction probe drops its positional casts (mapped-tuple result types).

GROK-20602: Added a details-expand datetime pin to `Dapi: domain parity` — child rows' datetimes materialize as dayjs (via `detailDatetimeColumns`) and the instant equals the same row read top-level (guards the naive-wire tz-shift class).

GROK-20601: WO-4b — `Dapi: domain lifecycle` gains schema-level grants round-trip and a real permissions assertion on shareColumn's grant half; the dapi2 domains smoke test now pins the namespace REMOVAL (typed surface reached parity — dapi2 is type-only now, dapi2Init survives).

GROK-20601: Added the `Dapi: domain lifecycle` suite — user-schema round-trip (createSchema → apply dryRun/apply → table client data → whole-schema audit with ddl + row events → delete; typed `DomainVersionConflictError` on stale ifVersion and `destructive-confirmation-required` with the plan attached; self-skips without the CreateDomainSchema privilege), table grants list/grant/revoke round-trip, and column security (shareColumn / idempotent restrictColumn / restoreColumnVisibility).

GROK-20600: Added `deleteWhere` coverage to the `Dapi: domain parity` suite — string + tree filters with a `hasMore` drain loop, empty filter → `DomainValidationError`, restrict → `DomainRestrictError` with zero deletions (Grit fixture, self-skips); `Dapi: domains batch` teardown now uses one `deleteWhere` instead of the query + per-row delete loop.

GROK-20599: Added the `Dapi: domain parity` suite — count/exists (string + tree filters), first, getByKey (composite hit / ambiguity → null), upsert (inserted→updated with version bump; typed validation failures incl. the report-carrying case), fetchFields (requested shape, absent ids, empty ids → empty DataFrame), aggregateDf vs aggregate parity, watch round-trips (table + row), row-watch on an audit:false table (Inventory fixture, self-skips), auditLog newest-first with numeric wire ids.

GROK-20598: Added the `Dapi: domain errors` suite — typed error probes: stale-version update → `DomainVersionConflictError` with both versions, restrict delete (Grit fixture, self-skips without Grit) → `DomainRestrictError`, malformed filter → `DomainFilterError` (400), `errorOnDuplicate` insert → `DomainValidationError.isDuplicate`.

GROK-20598: Adapted the domain suites to the typed js-api domain surface (WO-1): condition-tree filters now compile without `as any`, facet results are per-kind typed (casts on heterogeneous batches), batch report rows / transaction results / audit entries are typed; no behavior changes.

GROK-20591: Added the `Dapi: domain filters` suite — batched facets round trip (categories/histogram/minMax/count under a condition-tree filter), FK-path filtering (dotted tree through `query` + path categories facet), and saved filter presets (save, list, in-place update, entity sharing via `grok.dapi.permissions`, delete) over the package's own `apitests` schema.

GROK-20572: Added the `Dapi: domain visual queries` suite — `TableQueryBuilder` (select/where/sortBy/limit and an aggregation) over the virtual `Domain` connection of the package's own `apitests` domain schema, plus the raw-SQL refusal. Self-skips when the schema or its connection is not deployed.

GROK-20322: Behavior change (WO-B15): a mutation that fails as a whole (rolled back, nothing applied) now REJECTS the `DbTable` promise with the SQL message — callers that polled `res.errorMessage` on a resolved result get a rejection instead. Added `Dapi: connector writes` cases: total-failure update rejects; `allOrNothing: false` per-row errors still resolve with survivors committed.

GROK-20322: Added `Dapi: connector ddl` suite for the `grok.data.db.ddl(...)` fluent DDL builder and `DbTable.uploadAs` (create-table-from-DataFrame): scratch-table lifecycle (create → insert → addColumn → dryRun dropTable → typed `DdlConfirmationRequiredError` with the plan attached → confirmed drop), uploadAs dryRun plan + 5k load with native-type spot checks, and a capability negative (PostgresDart, supportsDdl=false). Round trips self-skip on a stack without the DDL dispatch.

GROK-20358: Reworked the `Dapi: connector writes` object[] cases for the typed-DataFrame boundary: 5-row null preservation, a Date column round-tripping as datetime, int→float widening, an all-null-column loud error (client-side, runs in CI), and a 100k object[] going through the same bulk path. Follow-up: added a dayjs-column → datetime case and a `1e21` large-float-magnitude case (no longer throws).

GROK-20341: Made the `Dapi: connector writes` write-DB target configurable (defaults to the local compose demo `world` DB reachable from the grok_connect container; override via a `dgConnectorWritesDb` global) so the suite runs against a local stack instead of the hardcoded remote DB.

GROK-20341: Added `Dapi: connector writes` suite for the `grok.data.db.table(...)` structured-write surface (insert object[]/DataFrame bulk, upsert by keys, update/delete by where, capability-negative). Write round trips self-skip when the running grok_connect lacks `/mutate`.

GROK-20316: Added a `dapi2.domains` generated-client smoke test (`queryRows` over the wire) to the `Dapi: domains` suite.

GROK-20315: Added `Dapi: domains batch` suite for the phase-2 `grok.dapi.domains` surface (batch upload via CSV/DataFrame/d42/Parquet, upsert counts, partial-success reports, multi-entity transactions with `$ref` + rollback, aggregate, `queryDf` values and column tags) plus the `item_event` detail table in `databases/apitests/schema.json`.

GROK-20307: Added `Dapi: domains` suite for `grok.dapi.domains` (row CRUD, optimistic concurrency, business-key dedup, promote, audit) over the new `databases/apitests/schema.json` fixture schema.
* Tests: added the `Undo` category covering the new `grok.shell.undo/redo/canUndo/canRedo/onUndo/onRedo` and `DG.UndoService` surface — push/undo/redo round trip, multi-level LIFO, event firing, release of records on table close, and one-way records
* Tests: fixed `Undo: multi-level LIFO` / `onUndo fires` — `UndoService.contextCheck` only applies records whose context is the current table, so the tests now take the view before undoing (the Dart-side suite gets the same isolation by nulling `contextCheck`, which JS callers cannot do)
* Tests: split `Dapi: all data sources` into 14 per-source tests (`all data sources: queries` … `environments`) — the combined test summed 14 sequential round trips and drew one misleading 3s+ band in the stress report at high concurrency instead of showing per-endpoint latency
* GROK-20443: Node runner: default `--mode` flipped to `stress` — `grok stresstest` passes no mode flag, so the `functional` default silently tripled the nightly Stress-Tests workload (29 → 105 tests) and fired false regression alerts; functional runs opt in via `--mode=functional` (`npm run start-node` updated)
* Tests: fixed flaky `shell.route returns a View` — `shell.route` returns the current view, which is null when earlier tests close all views; the test now opens its own view and routes to it
* Tests: annotated UI-independent categories (dapi, dataframe, functions, bitset, valuematcher, property, stats, shell) with `node: true` — `grok test` now runs them headless in Node before the browser pass; view/dialog/grid tests keep `node: false`
* Tests: adopted DBTests' browser-bound tests as `DB: Data sync` / `DB: Client cache` / `DB: Benchmarks` (streaming into a table view, IndexedDB client cache, cached benchmarks) so DBTests runs fully headless

* Moved SPGI test datasets (`SPGI-linked1/2`, `SPGI_v2_infinity`) from the Demo files root into the package (`files/datasets`); viewer tests now read `System:AppData/ApiTests/datasets/...`

* Node runner: Expanded coverage from `dapi/` to all UI-independent categories (dataframe, functions, bitset, valuematcher, property, stats, shell); added `--mode functional|stress` (stress keeps the old stressTest-only behavior); per-file load errors no longer abort the run

* Security: rebuilt `apitests-docker-test1` on `python:3.12-alpine` (+ `apk upgrade`, refreshed pip/setuptools/wheel) to clear base-OS CVEs (expat/krb5/openssl/musl) and stale Python tooling.

Add `Dapi: entities.save (polymorphic)` tests covering Project and DataConnection round-trip via `grok.dapi.entities.save`.

Run the Node test runner (`package-test-node`) as ESM via tsx: dynamic `startDatagrok` import, ESM `loadTestFiles`, asset/dayjs shims in `node-test-loader/`, and a shared `test-package` `_package` holder so dapi tests no longer transitively load the browser suite. New `start-node`/`stress-node` scripts.

## 1.10.2 (2026-03-17)

Additonal setOptions tests

## 1.7.17 (2024-08-30)

Add tests for grok.dapi.files.list

## 1.7.9 (WIP)
