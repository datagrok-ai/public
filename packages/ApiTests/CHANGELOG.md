# API Tests changelog

## 1.10.3 (WIP)

GROK-20298: Added `Dapi: domain query state` — `DG.DomainQuery`: lossless URL round trip with the platform's list-element binding (reserved `view=`/`entity=` ignored, gaps closed), a live Domain View's `query` reconstructed and re-run to the same subset, `run()` matching `queryDf(toSpec())` value-for-value while recording a `DomainQuery` creation script, `fromBuilder` (AND conjuncts split per element) selecting the same rows, and typed failures for malformed element indices, non-integer limits, unparseable filter elements, and specs that need the function's own grammar parser.

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
