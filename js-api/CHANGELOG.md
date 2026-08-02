# datagrok-api changelog

## v.next

* GROK-20602: BREAKING: `grok api`-generated domain clients now type datetime columns (including `created_on`/`updated_on`) as dayjs `Dayjs` and pass `datetimeColumns` so JSON reads materialize dayjs objects; untyped `grok.dapi.domains.table('s.t')` clients keep ISO strings. Regenerate `src/generated/db.ts` and fix call sites that treated datetimes as strings.
* GROK-20601: BREAKING: Removed the `domains` namespace from the generated `grok.dapi2` client — the typed `grok.dapi.domains` surface reached parity (reads/writes, batch, transactions, deleteWhere, watch, audit, schema lifecycle, grants, column security). `dapi2` is currently type-only (`dapi2Init` remains for future generated routers); the OpenAPI yaml keeps every route.
* GROK-20601: Added `DomainSchemaClient` (`grok.dapi.domains.schema(name)`: manifest/apply/audit/delete + schema-entity grants) and table-scoped `grants`/`grant`/`revoke`/`shareColumn`/`restrictColumn`/`restoreColumnVisibility`.
* GROK-20600: Added `DomainTableClient.deleteWhere(filter, {limit})` — filtered bulk delete with `{deleted, hasMore}` report.
* GROK-20599: Added domain parity methods: `count`/`exists`/`first`/`getByKey`/`fetchFields`/`aggregateDf`/`upsert`/`auditLog`/`watch`/`unwatch`/`isWatching`.
* GROK-20598: Introduced the typed domain surface (`src/domains.ts`): condition-tree filter types accepted uniformly, typed results, and the `DomainError` class family (`DomainValidationError`, `DomainVersionConflictError`, `DomainRestrictError`, ...) replacing message matching.
