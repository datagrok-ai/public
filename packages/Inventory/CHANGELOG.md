# Inventory changelog

## v.next

* GROK-20602: BREAKING regen — `grok api` codegen v2: datetime fields are dayjs, `stock_movements.reason` is a literal union, adjustStock rides the typed `inventoryDb.transaction`, db.ts clients are lazy
* GROK-20318: Improved adjustStock — the quantity update and the stock_movements row now apply atomically via `grok.dapi.domains.transaction` (movement carries `moved_on`)
* GROK-20318: Introduced Inventory — reference app for hybrid domain schemas: items + stock_movements manifest with per-department jsonb property schemas, typed clients via `grok api` codegen, CSV/Parquet batch upsert by SKU, optimistic stock adjustments with 409 retry, stock-on-hand-by-location aggregation, column-security demo setup

## 0.0.1 (2026-07-04)
