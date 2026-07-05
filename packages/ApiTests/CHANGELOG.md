# API Tests changelog

## 1.10.3 (WIP)

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
