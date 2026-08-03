# API Tests changelog

## 1.10.3 (WIP)

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

Add `Dapi: entities.save (polymorphic)` tests covering Project and DataConnection round-trip via `grok.dapi.entities.save`.

Run the Node test runner (`package-test-node`) as ESM via tsx: dynamic `startDatagrok` import, ESM `loadTestFiles`, asset/dayjs shims in `node-test-loader/`, and a shared `test-package` `_package` holder so dapi tests no longer transitively load the browser suite. New `start-node`/`stress-node` scripts.

## 1.10.2 (2026-03-17)

Additonal setOptions tests

## 1.7.17 (2024-08-30)

Add tests for grok.dapi.files.list

## 1.7.9 (WIP)
