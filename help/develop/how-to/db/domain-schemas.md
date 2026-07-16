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
* **Standard UI**: browsing, search, create/edit dialogs, in-grid editing, import/export,
  sharing, watching, and history work automatically — see [Domains](../../../govern/catalog/domains.md).
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
| `schemas`              | —         | [Property schemas](#property-schemas-and-column-security) contributing dynamic columns          |
| `friendlyName`, `description` | — | Display metadata                                                                                |

### Column options

Column `type` is one of `string`, `int`, `float`, `bool`, `datetime`, `string_list`,
`ref` (foreign key to another table in the same manifest), `user`, or `group`
(foreign keys to platform users/groups, rendered with the matching pickers).

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

## Security modes

Each table declares how its rows are protected:

| Mode     | Typical use                              | Behavior                                                                                       |
|----------|------------------------------------------|-------------------------------------------------------------------------------------------------|
| `table`  | Lookup and reference tables (default)    | One permission check against the table itself: a View grant shows all rows, Edit allows writes  |
| `master` | Detail tables (issue → project, well → plate) | Each row inherits the security of the row it references through the `delegate` column; chains up to two hops deep |
| `row`    | Registration masters (studies, plates)   | Individual rows can be shared with users and groups; unshared rows follow the table-level grant (or stay hidden with `defaultRowVisibility: "none"`) |

Grants use the standard permissions (View, Edit, Delete, Share) on the schema, table, and
property-schema entities. Grant them from the UI (the table's **Sharing** pane) or
programmatically via `grok.dapi2.domains` (grant/revoke/list endpoints). In row mode, sharing
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

// The same query as a typed DataFrame (1M-row cap); columns carry semantic types,
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

## Typed clients

`grok api` (already part of standard package build scripts) detects
`databases/*/schema.json` and generates `src/generated/db.ts` with per-table row and insert
interfaces, column-name unions, and a ready-made client:

```ts
import {gritDb} from './generated/db';

const projects = await gritDb.project.query({sort: 'name'});   // ProjectRow[]
await gritDb.issue.insert({project_id: projects[0].id, number: 7, title: 'Typed!'});
// gritDb.issue.insert({}) — compile error: required columns are enforced
```

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
