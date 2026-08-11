---
title: "Domain schemas"
sidebar_label: "Domain schemas"
keywords:
  - domain schema
  - domain table
  - entity-mapped tables
  - row-level security
  - schema.json
  - CRUD
  - facets
  - default filters
---

:::note

This feature is in Beta

:::

A **domain schema** is a plugin-defined PostgreSQL schema whose tables are mapped to Datagrok
entities. Unlike [plain plugin databases](db-in-plugin.md), where you write your own SQL and
server code, domain tables get the full platform treatment out of the box:

* **Security**: row-level and column-level access control using standard Datagrok
  [groups and permissions](../../../govern/access-control/users-and-groups.md) — no separate
  security model to build.
* **Managed CRUD**: the server builds all queries, validates values, checks permissions, and
  writes an [audit trail](../../../govern/audit/audit.md) in the same transaction. You never
  write SQL or a backend.
* **Standard UI**: browsing, search, filtering, create/edit dialogs, in-grid editing,
  import/export, sharing, watching, and history work automatically — see
  [Domains](../../../govern/catalog/domains.md).
* **Typed API**: a generic JS client (`grok.dapi.domains`) plus generated TypeScript
  interfaces per table.

You declare tables once in a JSON manifest; Datagrok creates the database objects on package
deployment and upgrades them when the manifest changes.

For working examples, see the
[Grit](https://github.com/datagrok-ai/public/tree/master/packages/Grit) (issue tracker),
[Inventory](https://github.com/datagrok-ai/public/tree/master/packages/Inventory)
(batch upload, upsert, aggregation), and
[PlatesFixture](https://github.com/datagrok-ai/public/tree/master/packages/PlatesFixture)
(all security modes) packages.

## Declaring a schema

Add `databases/<schema>/schema.json` to your package. The directory name must match the
manifest's `name`. The manifest is validated against
[domain-schema.schema.json](https://github.com/datagrok-ai/public/blob/master/tools/domain-schema.schema.json).

```json
{
  "name": "grit",
  "version": "1.0.0",
  "tables": {
    "project": {
      "businessKey": ["key"],
      "columns": {
        "key":         {"type": "string", "required": true, "unique": true},
        "name":        {"type": "string", "required": true, "isName": true},
        "description": {"type": "string"}
      }
    },
    "issue": {
      "securityMode": "master",
      "delegate": "project_id",
      "businessKey": ["project_id", "number"],
      "columns": {
        "project_id": {"type": "ref", "ref": "project", "required": true, "onDelete": "cascade"},
        "number":     {"type": "int", "required": true},
        "title":      {"type": "string", "required": true, "isName": true},
        "status":     {"type": "string", "choices": ["open", "in progress", "resolved", "closed"], "default": "open"},
        "assignee":   {"type": "user"}
      }
    }
  }
}
```

Every table automatically gets the system columns `id` (UUID, also the row's entity id),
`version` (optimistic concurrency counter), `created_on`, `updated_on`, `author_id`, and
`is_deleted` (rows are soft-deleted by default). Do not declare them.

### Table options

| Option                 | Default   | Description                                                                                     |
|------------------------|-----------|-------------------------------------------------------------------------------------------------|
| `securityMode`         | `table`   | `table`, `master`, or `row` — see [security modes](#security-modes)                              |
| `delegate`             | —         | Master mode: the `ref` column whose target row's security applies                                |
| `promotion`            | `lazy`    | Row mode: when the row becomes an individually sharable entity (`lazy` = on first share, `eager` = on insert) |
| `defaultRowVisibility` | `table`   | Row/master modes: whether table-level View shows unshared rows (`none` hides them)               |
| `businessKey`          | —         | Natural-key column list: powers deduplication on insert, upsert matching, and search handles     |
| `audit`                | `true`    | In-transaction audit trail with before/after diffs; also enables row history and row-level watch |
| `softDelete`           | `true`    | Deletes mark `is_deleted` instead of removing rows                                               |
| `idempotency`          | `false`   | Adds the `idempotency_key` column for replay-safe creates                                        |
| `extensible`           | `false`   | Lets users add [their own columns](#extending-a-plugin-schema) to this table                     |
| `schemas`              | —         | [Property schemas](#property-schemas-and-column-security) contributing dynamic columns          |
| `filters`              | —         | [Default filters](#default-filters) shown on the table's filter panel                            |
| `friendlyName`, `description` | — | Display metadata                                                                                |

### Column options

Column `type` is one of `string`, `int`, `float`, `bool`, `datetime`, `string_list`,
`ref` (foreign key to another table in the same manifest), `user` or `group`
(foreign keys to platform users/groups, rendered with the matching pickers), or
`file` (a user-uploaded file in platform file storage).

A `file` column stores a `file://<connection>/<path>` string. Row dialogs render a
file input: picking a local file uploads it to the shared `System:DomainFiles`
storage under an unguessable per-value path, and grids render the file with its
name and size, with download and preview available from the **Context Panel**.
Anyone who obtains the stored path can access the file — access control comes from
the row and column security of the column that carries it.

```json
"issue": {"columns": {"attachment": {"type": "file", "friendlyName": "Attachment"}}}
```

| Option         | Description                                                                     |
|----------------|---------------------------------------------------------------------------------|
| `required`     | Rejects null values                                                             |
| `unique`       | Unique among live (not soft-deleted) rows                                       |
| `ref`          | Target table name (for `type: "ref"`); cross-plugin references are not allowed  |
| `onDelete`     | Referential action on delete: `cascade`, `restrict`, or `setnull`               |
| `min`, `max`   | Numeric range validation                                                        |
| `choices`      | Controlled dictionary of allowed values (renders as a combo box)                |
| `default`      | Default value for new rows and dialog inputs                                    |
| `isName`       | Marks the primary display-name column (one per table, `string` only); its value titles cards, tooltips, and entity views. Without it, a string column literally named `name` is used by convention |
| `semType`, `friendlyName`, `description`, `format` | Display and semantic metadata                     |

Validation runs twice with the same code: client-side in dialogs (instant feedback) and
server-side on every write (authoritative). Integrity constraints (`unique`, foreign keys,
`required`) are additionally enforced by the database.

### Property schemas and column security

Besides relational columns, a table can carry dynamic, strongly-typed values grouped into
**property schemas** — the mechanism behind column-level security:

```json
{
  "tables": {
    "item": {
      "columns": {"sku": {"type": "string", "required": true}},
      "schemas": ["chemistry", "procurement"]
    }
  },
  "propertySchemas": {
    "chemistry":   {"cas_number": {"type": "string"}, "hazard_class": {"type": "string", "choices": ["1", "2", "3"]}},
    "procurement": {"unit_cost": {"type": "float", "min": 0}, "reorder_point": {"type": "int"}}
  }
}
```

A user sees (and can edit) a column only when a group they belong to holds View (or Edit) on
one of the column's property schemas. In the example, chemists granted the `chemistry` schema
see `cas_number` and `hazard_class`, while procurement sees `unit_cost` — on the same rows.
Hidden columns never leave the server: they are absent from query results, exports, and
filters. Relational columns belong to the table's built-in "core" schema, which is granted to
all users on deployment.

### Default filters

Without configuration, the table's
[filter panel](../../../govern/catalog/domains.md#filtering) is constructed automatically
from the column types and value cardinality. To control it, declare a per-table `filters`
section — it replaces the automatic selection for that table, and the declared filters appear
pre-opened, in order (users can still add more from the panel):

```json
"filters": [
  {"column": "sample_state"},
  {"column": "volume", "type": "histogram", "bins": 24},
  {"column": "plate_id.barcode", "label": "Plate barcode"},
  {"column": "created_on", "type": "range"}
]
```

| Key      | Description                                                                                       |
|----------|----------------------------------------------------------------------------------------------------|
| `column` | A declared column, a property-schema column, a system column (such as `created_on`), or a dotted reference path like `plate_id.barcode` |
| `type`   | `categories`, `histogram`, `range`, `text`, or `bool`. Omit to pick automatically from the column type |
| `bins`   | Histogram bucket count (1–200). Requires an explicit `histogram` or `range` type                    |
| `label`  | Display caption for the filter                                                                      |

The section is validated on deployment, and every problem fails the publish with a distinct
error:

* Unknown columns, unknown descriptor keys, and duplicate columns are rejected.
* `type` must match the column type: `histogram` and `range` apply to `int`, `float`, and
  `datetime` columns, `text` to `string`, `bool` to `bool` (`categories` applies to any).
* `bins` without an explicit `histogram` or `range` type is rejected.
* In a dotted path, every segment but the last must be a `ref` column, and paths are capped
  at three hops. Paths point forward only — from the table to the tables it references — and
  support the `categories` type only.

## Security modes

Each table declares how its rows are protected:

| Mode     | Typical use                              | Behavior                                                                                       |
|----------|------------------------------------------|-------------------------------------------------------------------------------------------------|
| `table`  | Lookup and reference tables (default)    | One permission check against the table itself: a View grant shows all rows, Edit allows writes  |
| `master` | Detail tables (issue → project, well → plate) | Each row inherits the security of the row it references through the `delegate` column; chains up to two hops deep |
| `row`    | Registration masters (studies, plates)   | Individual rows can be shared with users and groups; unshared rows follow the table-level grant (or stay hidden with `defaultRowVisibility: "none"`) |

Grants use the standard permissions (View, Edit, Delete, Share) on the schema, table, and
property-schema entities. Grant them from the UI (the table's **Sharing** pane) or
programmatically — `grok.dapi.domains.table('s.t').grants()/grant()/revoke()` for a table,
`grok.dapi.domains.schema('s').grants()/grant()/revoke()` for the schema entity. Note that
schema-level grants gate schema *operations* (apply requires Edit, delete requires Delete,
sharing requires Share) — they do not grant access to row data; grant per table for that.
In row mode, sharing
a row for the first time transparently *promotes* it to a full platform entity — after that it
behaves like any entity: sharing dialog, favorites, comments, global search.

On deployment, the publishing user receives the full permission set on the schema, its tables,
and property schemas.

## Deployment and upgrades

The schema deploys when the package is published (`grok publish`) and updates on every
subsequent publish:

* **Additive changes** (new tables, new columns, new property schemas, metadata edits) are
  applied automatically.
* **Destructive changes** (dropping or retyping columns) are refused with an error — declare
  explicit statements under the manifest's `migrations` key to apply them deliberately.
* **Uninstalling** the package keeps the data (the schema is orphaned); reinstalling re-adopts
  it. A full purge is available to administrators.

## Extending a plugin schema

Every deployment is different, and users routinely need one more field. A schema can invite
that instead of forcing a fork: opt in from the manifest, and users you trust add their own
tables and columns to your database at runtime.

```json
{
  "name": "grit", "version": "1.2.0",
  "extensible": {"tables": true},
  "tables": {
    "issue": {"extensible": true, "columns": {"title": {"type": "string"}}}
  }
}
```

* `"extensible": {"tables": true}` at the root lets users add **their own tables**.
* `"extensible": true` on a table lets them add **their own columns** to it.

Both are off by default, and both are yours to revoke: turning a flag off blocks new
extensions while everything already added keeps working.

Users also need the **Extend** permission on the schema entity — grant it from the schema's
`Share...` dialog, or with
`grok.dapi.domains.schema('grit').grant(groupId, 'Extend')`. Extend is a schema-level
permission; it does not by itself grant access to any row data.

What users may do is deliberately narrow:

* Their own **tables** get the full manifest vocabulary, and they manage what they created.
* Their own **columns on your tables** stay optional and non-unique, never join the business
  key, never become the display-name column, and may only `restrict` or `setnull` on delete —
  your existing rows can never be invalidated or deleted by someone else's column.
* Everything **you** declared is immutable to them: your columns, your table metadata, your
  property schemas.

Republishing your package is safe. The deploy diff sees only your objects: user tables and
columns are retained untouched, and a plugin column or table that would collide with one of
theirs is refused with a named error rather than silently adopting their data. Uninstalling
orphans the schema and keeps their work; reinstalling re-adopts it.

Two things to know about how extension columns are stored:

* They live under an `x_` **physical prefix** in PostgreSQL (`x_customer_id`) so they can
  never collide with a column you add later. Every API surface — insert, query, filter, patch,
  batch, d42, the audit trail — uses the declared name (`customer_id`). The prefix does surface
  in raw PostgreSQL messages that quote an identifier verbatim, such as a constraint-violation
  error naming `fk_issue_x_customer_id`. Because of the prefix, **no column you declare may
  start with `x_`**; the manifest validator rejects it.
* They exist only in the server registry, so `grok api` **never generates them** into your
  typed clients — your generated code keeps describing exactly what your manifest declares.
  Users reach their columns through the generic client,
  `grok.dapi.domains.table('grit.issue')`.

Extension applies use the same endpoint as any other schema change
([`schema.apply`](#schema-lifecycle-grants-and-watching)) with an `extend` section, and the
same dry-run/confirm flow for anything destructive:

```ts
await grok.dapi.domains.schema('grit').apply({
  tables: {customers: {columns: {name: {type: 'string', isName: true}}}},
  extend: {issue: {columns: {customer_id: {type: 'ref', ref: 'customers', onDelete: 'setnull'}}}},
});
```

The `extend` section is full state for the user's own columns of that table: omitting one they
previously added proposes its drop. Optimistic concurrency runs on the schema's extension
counter (`ifVersion` ↔ `extVersion`, echoed in the dry-run plan) — your package's own version
is never touched by an extension. See `dapi/domains/extend-schema.js` in ApiSamples.

## Working with data from JS

The generic client is `grok.dapi.domains`. Address tables as `'<schema>.<table>'`:

```ts
const issues = grok.dapi.domains.table('grit.issue');

// Filtered, sorted query (same grammar as entity search; 10k-row cap)
const open = await issues.query({
  filter: 'status = "open" and title starts "Crash"',
  sort: '!created_on',
  limit: 100,
});

// The same query as a typed DataFrame (10M-row cap); columns carry semantic types,
// choices, and property tags, so grids render them correctly out of the box
const df = await issues.queryDf({filter: 'assignee = @current'});

// Single row by id; null when absent or not visible
const issue = await issues.get(id);

// Insert (per-row report; business-key duplicates are reported, not duplicated)
const [r] = await issues.insert({project_id: projectId, number: 42, title: 'Crash on save'});

// Update with optimistic concurrency: fails with a version conflict if the row
// changed since you read it
await issues.update(r.id, {status: 'resolved'}, {version: issue.version});

// Soft delete (declared referential actions are enforced)
await issues.delete(r.id);

// Aggregate over the rows and columns the caller can see
const counts = await issues.aggregate({
  groupBy: ['status'],
  measures: [{fn: 'count'}],
});

// Row history (audit-enabled tables)
const trail = await issues.audit(r.id);
```

For a typed DataFrame of aggregation results, use `DomainQuery` in aggregate mode (below):
one row per group, columns typed from the resolved measures. A boolean `groupBy` column
comes back as strings (`'true'`/`'false'`), and aggregate outputs carry no property tags —
they are not registry columns.

The same select and aggregate surface is packaged as the `DomainQuery`
[function](../../../datagrok/concepts/functions/functions.md) — the reproducible query the
platform records behind **Open in Table View** and data-synced domain dashboards (see
[Queries and dashboards](../../../govern/catalog/domains.md#queries-and-dashboards)). It is
callable from JS like any function:

```ts
// Select mode: typed DataFrame, ref-column captions applied (10M-row ceiling)
const open = await grok.functions.call('DomainQuery', {
  schema: 'grit', table: 'issue',
  filters: ['status = "open"'], orderBy: ['!created_on'], limit: 100,
});

// Aggregate mode: one row per group
const byStatus = await grok.functions.call('DomainQuery', {
  schema: 'grit', table: 'issue', groupBy: ['status'], aggregations: ['count'],
});
```

Note that only frames opened through the UI (**Open in Table View**, or a `DomainQuery` run
from the console) record their generation script for project data sync; a frame fetched
programmatically via `grok.functions.call` carries no creation script — save it as static
data, or re-run the query in the UI to make it data-sync-ready.

Reads only ever return rows and columns the current user can see — there is no way to opt out
client-side. Writes are validated, permission-checked, and audited server-side.

### Batch upload and upsert

`batch` loads large payloads efficiently (a million-row CSV registers in under a minute) and
merges by the table's business key in upsert mode:

```ts
const report = await items.batch(csvString, {mode: 'upsert', allOrNothing: false});
// {inserted, updated, skipped, errorCount, rows: [{index, id, status, errors?}, ...]}
```

Accepted payloads: a `DG.DataFrame`, a CSV string, an array of row objects, or raw bytes
(`d42`, or `parquet` converted client-side via the Arrow package). With `allOrNothing: true`
(the default) any bad row aborts the whole batch; with `false`, good rows are applied and bad
ones are reported per row.

### Transactions

Multiple operations — including across tables of one schema — commit or roll back atomically.
An op can name its new row's id with `ref` for later ops to reference:

```ts
await grok.dapi.domains.transaction('grit', [
  {op: 'insert', table: 'project', ref: 'p', values: {key: 'GRIT', name: 'Grit'}},
  {op: 'insert', table: 'issue', values: {project_id: '$p', number: 1, title: 'First issue'}},
]);
```

### Expansions

Fetch related rows in one query: `expand: ['project_id']` adds the master row's columns
prefixed `project_id.<name>`; `expand: ['details:comment']` adds capped child-row arrays
(JSON queries only). The expanded table's own row and column security applies.

### Facets and saved filters

`facets` powers filter panels: one batched request computes category counts, histogram
buckets, value ranges, the row count, and column profiles in a single round trip. Category
counts, buckets, and `count` are computed under the passed `filter` with the conditions on
each facet's own column stripped, so a filter control shows counts under all *other* filters
(classic faceted search). The exception is the stable-axis rule: `minMax`, `plan`, and
histogram *bounds* are computed under your access predicate only, ignoring `filter` — a
narrowing filter never re-derives a filter's axis. All results respect row-level access and
column security, so two users can legitimately get different counts for the same data.

```ts
const wells = grok.dapi.domains.table('plates.plate_well');

const res = await wells.facets({
  filter: [{property: 'sample_state', operator: '=', value: ['filled', 'dosed']}],
  facets: [
    // → {categories: [{value, display?, total, filtered}], hasMore}
    {id: 'state', kind: 'categories', column: 'sample_state', limit: 100},
    // Dotted reference path (up to 3 hops); groups by id, `display` carries the name
    {id: 'plate', kind: 'categories', column: 'plate_id.barcode', search: 'P-1'},
    // → {min, max, buckets, nulls}; datetime bounds are ISO-8601 strings
    {id: 'vol', kind: 'histogram', column: 'volume', bins: 24},
    // → {count}
    {id: 'n', kind: 'count'},
    // → {columns: [{name, distinct, min?, max?}]} — for choosing filter types
    {id: 'plan', kind: 'plan', columns: ['sample_state', 'volume']},
  ],
});
grok.shell.info(`${res.facets['n'].count} rows, buckets: ${res.facets['vol'].buckets}`);
```

At most 32 facets per request. Category lists are capped (default 100, `hasMore` set when
more remain) — narrow with `search` (compiled server-side as a bound substring match)
instead of raising the cap.

Saved filter presets are small shareable entities carrying the filter panel's state maps
verbatim (the shape `DG.FilterGroup` saves):

```ts
const preset = await wells.filters.save('Filled wells',
  {'Sample state': {type: 'categorical', column: 'sample_state', selected: ['filled']}});
const presets = await wells.filters.list();  // presets visible to the caller, by name
await wells.filters.delete(preset.id);
```

Saving with `{id}` updates a preset in place, preserving its original author. Share a preset
with users or groups like any other entity.

## Typed clients

`grok api` (already part of standard package build scripts) detects
`databases/*/schema.json` and generates `src/generated/db.ts` with per-table row and insert
interfaces, column-name unions, expand maps, a typed transaction union, and a lazy
per-schema client:

```ts
import {gritDb, IssueStatus} from './generated/db';

const projects = await gritDb.project.query({sort: 'name'});   // ProjectRow[]
await gritDb.issue.insert({project_id: projects[0].id, number: 7, title: 'Typed!'});
// gritDb.issue.insert({}) — compile error: required columns are enforced
```

The generated surface is truthful about the wire:

* **Datetimes are dayjs.** Declared datetime columns and `created_on`/`updated_on` are typed
  `Dayjs` and materialize as dayjs objects on JSON reads, including expanded master fields
  and detail child rows (inserts also accept ISO strings). Untyped `table('s.t')` clients
  keep ISO strings. Regenerating db.ts across this change is breaking — fix call sites that
  treated datetimes as strings (`a.created_on.valueOf()` instead of `localeCompare`).
* **`choices` columns are literal unions** (`IssueStatus = 'open' | 'in progress' | ...`)
  used in both row and insert types — a typo in a status value no longer compiles.
* **Column names and expand keys are compile-checked** through the client generics: filter
  conditions, `columns`, `groupBy`, and `expand` reject unknown names.

### Fluent queries and bound conditions

Bare `query()` returns an awaitable builder; `query(spec)` is unchanged:

```ts
const top = await gritDb.issue.query()
  .where('project_id', '=', projectId)
  .where({status: 'open'})                        // equality map, AND-combined
  .orderBy('number', true)
  .top(20);

const one = await gritDb.issue.query().where('number', '=', 7).first();  // row | null
const df = await gritDb.issue.query().where({status: 'open'}).df();      // typed DataFrame
const n = await gritDb.issue.query().where({status: 'open'}).count();
```

Condition values are **bound server-side, never interpolated** — any string value works,
including apostrophes that the filter-string grammar cannot express:

```ts
await gritDb.project.query().where('name', '=', "O'Brien's project");
await gritDb.issue.query({filter: DG.or(
  DG.cond('status', '=', 'open'), DG.cond('priority', '=', 'critical'))});
```

`.expand('details:comment')` types the child arrays into the awaited rows. `.select(...)`
narrows the projection — system columns always ride along, and call it before `.expand()`.

### Typed errors

Failures are `DG.DomainError` subclasses discriminated by class and `code` — never match
message text:

```ts
try {
  await gritDb.issue.update(id, {status: 'resolved'}, {version});
} catch (e) {
  if (e instanceof DG.DomainVersionConflictError)
    grok.shell.info(`expected v${e.expectedVersion}, current v${e.currentVersion}`);
}
```

The family: `DomainValidationError` (per-row `.rows`, `.isDuplicate`),
`DomainVersionConflictError`, `DomainRestrictError`, `DomainFilterError`,
`DomainForbiddenError`, `DomainNotFoundError`, `DomainManifestValidationError`. A failed
transaction carries `.opIndex` — the index of the failing op.

### Optimistic concurrency

```ts
const saved = await gritDb.project.save({key: 'GRIT', name: 'Grit'}); // insert-or-update
await gritDb.issue.updateWithRetry(id, (fresh) =>
  fresh.status === 'open' ? {priority: 'high'} : null);   // null skips the write
await DG.retryOnVersionConflict(async () => {/* fresh read + transaction write */});
```

`save` addresses rows by identity — a business-key duplicate applies your values to the
existing row under a versioned update. An idempotency-key replay applies nothing — the
original insert already did — and resolves the existing row's fresh version.
`updateWithRetry` re-reads and retries on conflict
(default five retries after the initial attempt). Typed transactions get per-op result
types from a tuple ops literal: `const [upd, ins] = await gritDb.transaction([...]);`.

### Bulk delete

```ts
const stamp = `test-${Date.now()}`;
while ((await gritDb.issue.deleteWhere(DG.cond('title', 'like', stamp + '%'))).hasMore);
```

Soft-deletes up to 1000 matching rows you may delete per call, oldest first, in one
transaction. Declared referential actions run per row, and a `restrict` reference rejects
the whole call with a `DomainRestrictError`. Never write per-row delete loops.

### Schema lifecycle, grants, and watching

```ts
await grok.dapi.domains.createSchema('inv', {friendlyName: 'Inventory'});
const handle = grok.dapi.domains.schema('inv');
await handle.apply({tables: {/* manifest fragment */}}, {dryRun: true}); // change plan
await handle.apply({tables: {/* manifest fragment */}});
const events = await handle.audit({limit: 50});  // row + ddl history, newest first
await handle.delete();                           // full purge
```

Table-scoped access control lives on the table client: `grants()`,
`grant(group, permission)`, `revoke(group)`, and column security via
`shareColumn`/`restrictColumn`/`restoreColumnVisibility`. Schema-entity grants
(`handle.grant(...)`) gate schema *operations* — apply requires Edit, delete requires
Delete, sharing requires Share — and do not grant access to row data. Grant per table for
that.

`watch()`/`unwatch()`/`isWatching()` subscribe the current user to change notifications for
a table or one row (row watch requires the table's audit trail). `audit(id)` reads a row's
history, `auditLog({limit})` the table-wide trail.

### Samples

Runnable in the platform's samples gallery:
[crud](https://public.datagrok.ai/js/samples/dapi/domains/crud),
[typed-client](https://public.datagrok.ai/js/samples/dapi/domains/typed-client),
[aggregate](https://public.datagrok.ai/js/samples/dapi/domains/aggregate),
[batch](https://public.datagrok.ai/js/samples/dapi/domains/batch),
[transaction](https://public.datagrok.ai/js/samples/dapi/domains/transaction),
[filters](https://public.datagrok.ai/js/samples/dapi/domains/filters),
[dataframe](https://public.datagrok.ai/js/samples/dapi/domains/dataframe),
[idempotency](https://public.datagrok.ai/js/samples/dapi/domains/idempotency),
[schema](https://public.datagrok.ai/js/samples/dapi/domains/schema).

## Customizing the UI

The default UI (cards, tooltips, context panel, entity view) works for every table with no
code. Plugins customize it through the standard mechanisms only — rows of each table carry
the semantic type `<schema>.<table>`:

* **Custom rendering and views**: register an
  [ObjectHandler](register-identifiers.md) for the semantic type — it overrides the
  default cards, tooltips, context panel, and entity view for that table only:

  ```ts
  class IssueHandler extends DG.ObjectHandler {
    get type() { return 'grit.issue'; }
    // renderCard, renderTooltip, renderProperties, renderView, ...
  }
  DG.ObjectHandler.register(new IssueHandler());
  ```

* **Context actions**: any package function with an input of semantic type
  `<schema>.<table>` appears in the row's **Actions** pane automatically.
* **Info panels**: `#panel`-tagged functions with the same input appear in the
  [Context Panel](../../../datagrok/navigation/panels/panels.md#context-panel).
* **Search patterns**: claim [identifier patterns](register-identifiers.md) (like
  `GRIT-123`) in the handler, and they resolve from global search.

See also:

* [Domains](../../../govern/catalog/domains.md) — the user-facing guide
* [Plugin Postgres databases](db-in-plugin.md) — raw SQL storage without entity mapping
* [Access data](access-data.md) — connections and queries
