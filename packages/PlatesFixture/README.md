# PlatesFixture

A minimal engine-correctness fixture for **entity-mapped domain schemas** — not a demo.
The `databases/plates/schema.json` manifest exercises, in one small schema, everything the
domain engine must get right:

* all three security modes — `plate_type` (`table`), `plate` (`row`, lazy promotion),
  `plate_well` (`master`, delegating to its plate)
* constraints — required columns, `unique`, `min`, `choices`, business keys
  (dedup-on-insert), idempotency keys
* referential actions — `cascade` (wells with their plate), `restrict`
  (plate types referenced by plates)
* jsonb property schemas (`bio_plate`, `pricing`) gating dynamic columns
* the per-table `audit` flag (`plate_type` runs audit-off)

`src/tests/domain-fixture-tests.ts` runs these paths in CI through `grok.dapi.domains`.
Cross-user visibility and predicate tests live server-side in
`core/server/datlas/test/services/domain_*.dart`.

See `core/docs/features/db-schemas/ARCHITECTURE.md` (§9.4) in the Datagrok monorepo for the design.
