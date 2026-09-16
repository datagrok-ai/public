# The domain backend contract

What every `DomainBackend` must agree on — the platform one (`DgDomainBackend` over
`grok.dapi.domains`, dg) and the in-memory one (`MemoryDomainBackend`, core). A `DomainSource`
and the controls over it are written against this list, not against either backend; the memory
backend is what the headless tests and the gallery run, so a case the platform answers
differently is a bug in one of them, not a difference to code around. Pinned headless in
`tests/memory-domain.test.js`, `tests/domain-source.test.js` and `tests/domain-session.test.js`
(the optional members additionally in `tests/domain-tree.test.js`, `tests/domain-bulk.test.js`,
`tests/domain-import.test.js` and `tests/domain-live.test.js`); the platform side is pinned live
by the U2Demo `U2: domain source`, `U2: domain session`, `U2: domain trash`, `U2: domain bulk` and
`U2: domain import` categories over `apitests.item`, and by the ApiTests `Dapi: domain trash`,
`Dapi: domain bulk` and `Dapi: domain hierarchy` categories against the server itself.

**Optional members** (`audit`, `restore`, `updateWhere`, `ancestors`, `probe`, `batch`) are the
seam's escape hatch, not a licence to diverge: a backend that declares one answers it exactly as
the table below says, and a control that needs one checks for it and refuses by name rather than
degrading silently. The memory backend mirrors all of them except `batch`.

## Table handle

| Case | Both answer |
|---|---|
| `properties` | the system columns first — `id`, `version`, `created_on`, `updated_on`, `author_id` (semType `User`) — then the declared columns, in declaration order; a `ref` column carries its target address as `semType`; `user`/`group` columns `User`/`Group` |
| `info` | `nameColumn` (the `isName` column, else a `name` string column, else null), `businessKey`, `singularName`, `pluralName`; `searchableColumns` (the `searchable: true` columns, else `[nameColumn]`, else `[]`); `constraints` (the `{expr, message?}` entries of the schema by name — SQL `check`s stay server-side); `refFilters` (`{column: expr}` for every ref column declaring a `filter`); `permissions` (the custom names the schema declares, list or map form); `childTables` (`{schema, table, fkColumn, label}` for every ref column pointing at this table, `label` the column's caption); `hierarchy` + `parentColumn` (the `hierarchy: true` declaration and its one self-referencing ref column; absent/false and null otherwise) |
| `access()` | `{can, fields}`: `can.<capability>` for the five plus every custom permission (`can.<name>`, granted by default in memory); `fields[col]` for every column the caller may see — `readonly` for the system columns and autoNumber, `editable` for the rest. Column security only: a row's writability comes from `can.edit`/`can.insert` and the `~can_*` columns |
| `restore(id)` | a soft-deleted row brought back: `~is_deleted` off, `version + 1`, the `'undelete'` audit op; refused with `validation` naming the ref column where the row points at a row still deleted. Optional on the seam (dg fills it from `client.restore`) — a backend that does not declare it cannot answer a `deleted` query either, and `DomainSource` refuses one over it by name |
| `audit(id)` | the row's history, oldest first: `{id, tx_id, op, actor_id, ts, before, after}` per op that touched it — every op of one transaction shares `tx_id`; the memory backend records `actor_id: null` and nothing for seeded rows. Optional on the seam (dg fills it from `client.audit`) |
| `updateWhere(filter, values, {limit})` | `values` written into every LIVE row the filter matches that the caller may edit, as ONE transaction; answers `{updated, hasMore}`. A missing or empty filter is refused (`validation`) — there is no "update the whole table" — and so is an empty `values`; a column the caller may not write is refused BEFORE anything is written; the selection is capped at `limit` (default and maximum 1000, clamped to `[1, 1000]`) and `hasMore` says the filter matched past the cap. Rows are patched one by one inside the transaction, so validation, the version step and the audit line are per row, and a row outside the caller's Edit predicate is silently not selected rather than refused. Optional on the seam (dg fills it from `client.updateWhere`) |
| `ancestors(id)` | the row's parents along `info.parentColumn`, **root first**, **excluding the row itself** — so a root answers `[]`, and so does a row the caller cannot see (no oracle). Only a table declaring `hierarchy` answers; anything else is a `filter` refusal on BOTH backends. The chain stops at the first invisible ancestor, at a cycle, and at depth 64. Optional on the seam (dg fills it from `client.pathTo`) |
| `probe(spec)` | `{count, last}` — how many rows match `filter`/`search`/`deleted` and the newest `updated_on` among them (`null` over an empty match), in ONE request and never a page of rows. What a `live` source polls; a backend that does not declare it is never polled, so `live: true` over one is inert |
| `batch(rows, options)` | whole rows uploaded in one call, which is where the `'upsert'` business-key merge, `allOrNothing` and `errorOnDuplicate` live, and which answers the per-row report `domains.import` renders. Optional and NOT mirrored in memory: the memory backend has no `batch`, so `domains.import` over it falls back to the same rows as one transaction with a business-key upsert — the report is synthesised from the transaction's results |

## Query

| Case | Both answer |
|---|---|
| `filter` | a smart-filter string or the canonical condition tree; the same grammar, evaluated by the same rules (`Filters`); an unknown column or a parse error is a `validation` refusal. A column reference (`end_date >= start_date`) compares row-wise under the six comparators, null on either side false; a `$param` is bound by the caller (`Filters.bind`) before the query — an unbound one is refused |
| `search` | a case-insensitive substring over `info.searchableColumns`, ORed across them and ANDed with `filter`; a table with no searchable column refuses with `validation` |
| `sort` | `'col,!col'`; nulls last ascending, first descending (Postgres) |
| `limit`/`offset` | a page; `limit: 0` answers no rows (an empty or draft source) |
| `columns` | a projection; absent, every column |
| `withAccess` | adds `~can_edit`, `~can_delete`, `~can_share` (and `~can_<name>` per custom permission) to every row: booleans on a row-mode table; off row mode `~can_edit`/`~can_delete` are the table's answer and `~can_share` is **null** in JSON and **absent** from a frame (a bool column holds no null) — `Access.row` treats both as "not carried" |
| `deleted` | which rows a query answers: `'exclude'` (the default — live rows only), `'include'` or `'only'` (the trash). Anything but `'exclude'` projects `~is_deleted` with every row (a bool column on a frame), and the filter, the search and the paging still apply |
| `count(filter, search, deleted)` | the total under the same filter, search and deleted mode |
| `frame(spec)` | the rows as a frame with the writer attached (`DomainFrameLike`): the frame carries `~state` and the access columns as service columns; `append` adds a page into the same frame; `dispose` detaches the writer |

## Transaction

| Case | Both answer |
|---|---|
| Atomicity | all or nothing across every table and schema the batch touches — one refused op refuses the batch, and nothing before it landed |
| Targets | `op.table` is `'<table>'` in the transaction's schema or `'<schema>.<table>'`; the memory backend resolves a bare name against its own schema |
| Ordering | ops run in a stable topological order by request index: a `$ref` use after the insert that declares it (forward references allowed), a child table's delete before its parent's (a ref column from the child to the parent, none back); anything unconstrained keeps the request order. `results[i]` and every error's op index refer to the REQUEST index |
| `insert` | `id`, `version: 1`, `created_on`, `updated_on`, `author_id` (the caller) are stamped; a supplied `id`/`author_id` is ignored by the server (kept by the memory backend for seeding) |
| `update` | `version + 1`, `updated_on`; with `expectedVersion` a mismatch is a `version-conflict` refusal naming both versions |
| `delete` | SOFT: the row keeps its place in storage with `is_deleted` set, `version + 1`, and leaves every query that excludes deleted rows (the writer drops it from the frame all the same). A missing row is `not-found`; a row still referenced by a LIVE child row is refused (`validation` in memory — the FK veto the child-first ordering beats), and a child already in the trash holds nothing back |
| Validation | `required` (empty is null, undefined or `''`), `choices`, `min`/`max` — refused as `validation` naming the column; the memory backend and the platform editor report the same messages per cell before the batch is sent (`Value can't be empty`, `Must be one of: …`, `Must be at least N`, `Must be at most N`) |
| `$ref` | an insert may name itself with `ref`; any op's `'$<ref>'` (also inside lists) is that row's id, earlier or later in the batch; `'$$'` escapes a literal `$`; an unknown reference, a duplicate `ref` or a cycle among references is `bad-ref` (the cycle names its first op) |
| Draft ids | the writer stamps `~new:<uuid>` into a draft's `id` cell (`Rows.draftId()`); its `buildOps` names the insert `ref: <draftId>` and rewrites every value equal to a draft id — its own or another writer's — into `'$<draftId>'`, and a literal leading `$` into `'$$…'`; `id` is never sent. A landed batch answers `onSaved({assigned: {draftId: realId}})` |
| Write lock | every participating writer is CLOSED for the whole `saveAll` (`isSaving`): `setValue` and `markDeleted` are refused, adding a row throws — an edit typed while the transaction is in flight is not in the batch being sent and must not read clean afterwards. The platform closes its editors in the js-api `DomainSession`, the memory backend in `saveAll` |
| `saveAll(edits)` | every writer's pending batch concatenated as ONE transaction (each op's `table` qualified by its writer's table), the results sliced back to each writer's `applyResults`; resolves to whether it landed — the memory backend throws its `DomainBackendError`, the platform (the js-api `DomainSession`) answers false after its own conflict and validation dialogs |
| Non-writable columns | a value for a column the caller may not write (`fields` not `editable`, or the row's `~can_edit` false) is dropped by the writer before the batch is sent; the server refuses one that reaches it (server only — the memory backend does not check) |

## Errors

`DomainBackendError.code`: `not-found`, `validation`, `version-conflict`, `bad-ref`, `forbidden`
(server only) — the message is for the user.
