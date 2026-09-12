# The domain backend contract

What every `DomainBackend` must agree on — the platform one (`DgDomainBackend` over
`grok.dapi.domains`, dg) and the in-memory one (`MemoryDomainBackend`, core). A `DomainSource`
and the controls over it are written against this list, not against either backend; the memory
backend is what the headless tests and the gallery run, so a case the platform answers
differently is a bug in one of them, not a difference to code around. Pinned headless in
`tests/memory-domain.test.js` and `tests/domain-source.test.js`; the platform side is pinned
live by the U2Demo `U2: domain source` category over `apitests.item`.

## Table handle

| Case | Both answer |
|---|---|
| `properties` | the system columns first — `id`, `version`, `created_on`, `updated_on`, `author_id` (semType `User`) — then the declared columns, in declaration order; a `ref` column carries its target address as `semType`; `user`/`group` columns `User`/`Group` |
| `info` | `nameColumn` (the `isName` column, else a `name` string column, else null), `businessKey`, `singularName`, `pluralName` |
| `access()` | `{can, fields}`: `can.<capability>` for the five plus every custom permission; `fields[col]` for every column the caller may see — `readonly` for the system columns and autoNumber, `editable` for the rest. Column security only: a row's writability comes from `can.edit`/`can.insert` and the `~can_*` columns |

## Query

| Case | Both answer |
|---|---|
| `filter` | a smart-filter string or the canonical condition tree; the same grammar, evaluated by the same rules (`Filters`); an unknown column or a parse error is a `validation` refusal |
| `sort` | `'col,!col'`; nulls last ascending, first descending (Postgres) |
| `limit`/`offset` | a page; `limit: 0` answers no rows (a draft source) |
| `columns` | a projection; absent, every column |
| `withAccess` | adds `~can_edit`, `~can_delete`, `~can_share` to every row: booleans on a row-mode table; off row mode `~can_edit`/`~can_delete` are the table's answer and `~can_share` is **null** in JSON and **absent** from a frame (a bool column holds no null) — `Access.row` treats both as "not carried" |
| `count(filter)` | the total under the same filter |
| `frame(spec)` | the rows as a frame with the writer attached (`DomainFrameLike`): the frame carries `~state` and the access columns as service columns; `append` adds a page into the same frame; `dispose` detaches the writer |

## Transaction

| Case | Both answer |
|---|---|
| Atomicity | all or nothing — one refused op refuses the batch, and nothing before it landed |
| `insert` | `id`, `version: 1`, `created_on`, `updated_on`, `author_id` (the caller) are stamped; a supplied `id`/`author_id` is ignored by the server (kept by the memory backend for seeding) |
| `update` | `version + 1`, `updated_on`; with `expectedVersion` a mismatch is a `version-conflict` refusal naming both versions |
| `delete` | by id; a missing row is `not-found` |
| Validation | `required` (empty is null, undefined or `''`), `choices`, `min`/`max` — refused as `validation` naming the column; the memory backend and the platform editor report the same messages per cell before the batch is sent (`Value can't be empty`, `Must be one of: …`, `Must be at least N`, `Must be at most N`) |
| `$ref` | an insert may name itself with `ref`; a later op's `'$<ref>'` (also inside lists) is that row's id; `'$$'` escapes a literal `$`; a forward or unknown reference is `bad-ref` |
| Non-writable columns | a value for a column the caller may not write (`fields` not `editable`, or the row's `~can_edit` false) is dropped by the writer before the batch is sent; the server refuses one that reaches it (server only — the memory backend does not check) |

## Errors

`DomainBackendError.code`: `not-found`, `validation`, `version-conflict`, `bad-ref`, `forbidden`
(server only) — the message is for the user.
