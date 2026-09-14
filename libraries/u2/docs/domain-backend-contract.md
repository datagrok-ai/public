# The domain backend contract

What every `DomainBackend` must agree on — the platform one (`DgDomainBackend` over
`grok.dapi.domains`, dg) and the in-memory one (`MemoryDomainBackend`, core). A `DomainSource`
and the controls over it are written against this list, not against either backend; the memory
backend is what the headless tests and the gallery run, so a case the platform answers
differently is a bug in one of them, not a difference to code around. Pinned headless in
`tests/memory-domain.test.js`, `tests/domain-source.test.js` and `tests/domain-session.test.js`;
the platform side is pinned live by the U2Demo `U2: domain source` and `U2: domain session`
categories over `apitests.item`.

## Table handle

| Case | Both answer |
|---|---|
| `properties` | the system columns first — `id`, `version`, `created_on`, `updated_on`, `author_id` (semType `User`) — then the declared columns, in declaration order; a `ref` column carries its target address as `semType`; `user`/`group` columns `User`/`Group` |
| `info` | `nameColumn` (the `isName` column, else a `name` string column, else null), `businessKey`, `singularName`, `pluralName`; `searchableColumns` (the `searchable: true` columns, else `[nameColumn]`, else `[]`); `constraints` (the `{expr, message?}` entries of the schema by name — SQL `check`s stay server-side); `refFilters` (`{column: expr}` for every ref column declaring a `filter`); `permissions` (the custom names the schema declares, list or map form); `childTables` (`{schema, table, fkColumn, label}` for every ref column pointing at this table, `label` the column's caption) |
| `access()` | `{can, fields}`: `can.<capability>` for the five plus every custom permission (`can.<name>`, granted by default in memory); `fields[col]` for every column the caller may see — `readonly` for the system columns and autoNumber, `editable` for the rest. Column security only: a row's writability comes from `can.edit`/`can.insert` and the `~can_*` columns |
| `audit(id)` | the row's history, oldest first: `{id, tx_id, op, actor_id, ts, before, after}` per op that touched it — every op of one transaction shares `tx_id`; the memory backend records `actor_id: null` and nothing for seeded rows. Optional on the seam (dg fills it from `client.audit`) |

## Query

| Case | Both answer |
|---|---|
| `filter` | a smart-filter string or the canonical condition tree; the same grammar, evaluated by the same rules (`Filters`); an unknown column or a parse error is a `validation` refusal. A column reference (`end_date >= start_date`) compares row-wise under the six comparators, null on either side false; a `$param` is bound by the caller (`Filters.bind`) before the query — an unbound one is refused |
| `search` | a case-insensitive substring over `info.searchableColumns`, ORed across them and ANDed with `filter`; a table with no searchable column refuses with `validation` |
| `sort` | `'col,!col'`; nulls last ascending, first descending (Postgres) |
| `limit`/`offset` | a page; `limit: 0` answers no rows (an empty or draft source) |
| `columns` | a projection; absent, every column |
| `withAccess` | adds `~can_edit`, `~can_delete`, `~can_share` (and `~can_<name>` per custom permission) to every row: booleans on a row-mode table; off row mode `~can_edit`/`~can_delete` are the table's answer and `~can_share` is **null** in JSON and **absent** from a frame (a bool column holds no null) — `Access.row` treats both as "not carried" |
| `count(filter, search)` | the total under the same filter and search |
| `frame(spec)` | the rows as a frame with the writer attached (`DomainFrameLike`): the frame carries `~state` and the access columns as service columns; `append` adds a page into the same frame; `dispose` detaches the writer |

## Transaction

| Case | Both answer |
|---|---|
| Atomicity | all or nothing across every table and schema the batch touches — one refused op refuses the batch, and nothing before it landed |
| Targets | `op.table` is `'<table>'` in the transaction's schema or `'<schema>.<table>'`; the memory backend resolves a bare name against its own schema |
| Ordering | ops run in a stable topological order by request index: a `$ref` use after the insert that declares it (forward references allowed), a child table's delete before its parent's (a ref column from the child to the parent, none back); anything unconstrained keeps the request order. `results[i]` and every error's op index refer to the REQUEST index |
| `insert` | `id`, `version: 1`, `created_on`, `updated_on`, `author_id` (the caller) are stamped; a supplied `id`/`author_id` is ignored by the server (kept by the memory backend for seeding) |
| `update` | `version + 1`, `updated_on`; with `expectedVersion` a mismatch is a `version-conflict` refusal naming both versions |
| `delete` | by id; a missing row is `not-found`; a row still referenced by a child row is refused (`validation` in memory — the FK veto the child-first ordering beats) |
| Validation | `required` (empty is null, undefined or `''`), `choices`, `min`/`max` — refused as `validation` naming the column; the memory backend and the platform editor report the same messages per cell before the batch is sent (`Value can't be empty`, `Must be one of: …`, `Must be at least N`, `Must be at most N`) |
| `$ref` | an insert may name itself with `ref`; any op's `'$<ref>'` (also inside lists) is that row's id, earlier or later in the batch; `'$$'` escapes a literal `$`; an unknown reference, a duplicate `ref` or a cycle among references is `bad-ref` (the cycle names its first op) |
| Draft ids | the writer stamps `~new:<uuid>` into a draft's `id` cell (`Rows.draftId()`); its `buildOps` names the insert `ref: <draftId>` and rewrites every value equal to a draft id — its own or another writer's — into `'$<draftId>'`, and a literal leading `$` into `'$$…'`; `id` is never sent. A landed batch answers `onSaved({assigned: {draftId: realId}})` |
| Write lock | every participating writer is CLOSED for the whole `saveAll` (`isSaving`): `setValue` and `markDeleted` are refused, adding a row throws — an edit typed while the transaction is in flight is not in the batch being sent and must not read clean afterwards. The platform closes its editors in the js-api `DomainSession`, the memory backend in `saveAll` |
| `saveAll(edits)` | every writer's pending batch concatenated as ONE transaction (each op's `table` qualified by its writer's table), the results sliced back to each writer's `applyResults`; resolves to whether it landed — the memory backend throws its `DomainBackendError`, the platform (the js-api `DomainSession`) answers false after its own conflict and validation dialogs |
| Non-writable columns | a value for a column the caller may not write (`fields` not `editable`, or the row's `~can_edit` false) is dropped by the writer before the batch is sent; the server refuses one that reaches it (server only — the memory backend does not check) |

## Errors

`DomainBackendError.code`: `not-found`, `validation`, `version-conflict`, `bad-ref`, `forbidden`
(server only) — the message is for the user.
