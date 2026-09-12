# Managing a Datagrok Server with `grok s`

The `grok server` command (alias `grok s`) talks to a running Datagrok instance over
its public REST API. It is the right tool for **any server-management task** that
does not require rendering the UI: creating users and groups, sharing entities,
running functions, browsing files, hitting raw API endpoints, or scripting bulk
imports.

Use this instead of UI automation (Playwright, Selenium) whenever the goal is data
management — it is faster, idempotent, scriptable, and does not depend on the
browser or a logged-in session.

## When to use

| Task                                                        | Command                                                |
|-------------------------------------------------------------|--------------------------------------------------------|
| Create / update a user                                      | `grok s users save --json user.json`                   |
| Block / unblock a user                                      | `grok s users block <login>` / `users unblock <login>` |
| Create / update a group                                     | `grok s groups save --json group.json`                 |
| Add or remove users in a group                              | `grok s groups add-members <group> <user>...`          |
| List members of a group                                     | `grok s groups list-members <group>`                   |
| List the groups a user belongs to                           | `grok s groups list-memberships <user>`                |
| Share a connection / query / script / project with a group  | `grok s shares add <entity> <group> [--access View\|Edit]` |
| See who an entity is shared with                            | `grok s shares list <entity-id-or-name>`               |
| Create / update a connection                                | `grok s connections save --json conn.json`             |
| Test a connection                                           | `grok s connections test <id-or-name>`                 |
| Run a registered function                                   | `grok s functions run 'Pkg:fn(arg1,arg2)'`             |
| List functions with rich filters                            | `grok s functions list --type script --language python --package Chem` |
| List / browse files in a file share                         | `grok s files list "System:AppData" -r`                |
| Upload / download a table (CSV or d42)                      | `grok s tables upload <name> file.csv\|file.d42` / `tables download <name> -O out.csv` |
| Check whether a package is deployed                         | `grok s packages list --filter "MyPlugin"`             |
| Install / update / uninstall server plugins                 | `grok s packages install Chem Bio` / `packages update --all` / `packages uninstall Chem` |
| See installed vs latest plugin versions                     | `grok s packages outdated` / `packages versions Chem`  |
| Count entities                                              | `grok s users count --filter 'status = "active"'`      |
| See what fields an entity type has                          | `grok s describe connections`                          |
| Hit any undocumented endpoint                               | `grok s raw GET /users/current` / `raw POST <path> --data '{...}'` |
| Check server + per-module health                            | `grok s healthcheck [--module <name>]`                 |
| Bulk operations in one round-trip                           | `grok s batch <entity> <verb> --json items.json`       |
| Move entities dev to prod (bundle, or instance to instance) | `grok s pull ... --out ./bundle` / `grok s migrate ... --from dev --to prod` |
| Browse / query / edit domain-table rows                    | `grok s domains query grit.issue --filter 'status = "open"'` / `domains insert` / `domains upload` |
| Manage domain schemas (manifest, apply, grants)             | `grok s domains get grit` / `domains apply` / `domains grant grit.issue Chemists --access Edit` |

## Configuration

Servers and credentials live in `~/.grok/config.yaml`:

```yaml
default: local
servers:
  local:
    url: http://localhost:8888/api
    key: admin
  dev:
    url: https://dev.datagrok.ai/api
    key: <developer-key>
```

- `grok config add --alias <name> --server <url> --key <key>` writes a new entry.
- Add `--default` to make it the active server.
- Every `grok s ...` command accepts `--host <alias-or-url>` to override the default. The URL
  is the API base (`https://host/api`, or `http://host:8082` for a bare Datlas).

## Entity operations

### List / count / get / delete

```bash
grok s users list                              # default: table, 50 rows
grok s users list --filter "login = 'admin'"   # smart filter
grok s users list --limit 20 --offset 40       # third page of 20
grok s users count --filter 'status = "active"'
grok s connections list --output json
grok s packages list --filter "name:MyPlugin"  # table shows friendlyName
grok s groups get <id-or-name>
grok s users get alice.mendel                  # users: id or login
grok s connections delete "Admin:MyConnection" # id or namespace:name
```

Options that work on every entity:

| Flag                  | Meaning                                           |
|-----------------------|---------------------------------------------------|
| `--output table\|json\|csv\|quiet` | Output format; `quiet` prints ids only |
| `--filter "<expr>"`   | Server-side smart filter                          |
| `--limit <n>`         | Page size (default 50)                            |
| `--offset <n>`        | Skip first `n` rows                               |
| `--host <alias\|url>` | Target a specific server from your config         |

A missing entity, a rejected save (duplicate login or group name), or any other server-side
failure prints one line on stderr and exits 1 — never an error object on stdout.

`delete` is not available for every entity: **users cannot be deleted** (Datagrok has no user
deletion; `users delete` refuses and points at `users block`), and `functions delete` covers
scripts and queries only (package functions go away with their package). `count` has no server
endpoint for `reports`.

### Save from JSON

`grok s users save` and `grok s groups save` accept the same JSON shape that the
REST API returns on `get`. The simplest valid bodies:

```json
// user.json
{ "#type": "User", "login": "alice.mendel", "firstName": "Alice", "lastName": "Mendel", "status": "active" }
```

```json
// group.json
{ "#type": "UserGroup", "name": "Chemists", "friendlyName": "Chemists" }
```

```bash
grok s users save --json user.json
grok s groups save --json group.json --save-relations
```

Include an `id` to update an existing entity; omit it to create. To introspect the
shape of an existing entity, use `grok s <entity> get <id-or-name> --output json`.

### Connections

```bash
grok s connections save --json conn.json --save-credentials
grok s connections test "MyUser:MyConnection"         # by "author:name"
grok s connections test --json conn.json              # test before saving
```

See `public/packages/Chembl/connections/` for connection JSON examples.

## Group membership

`grok s groups add-members` is idempotent — it resolves each member (by login, name,
or UUID), compares against the group's current children, and only writes when
something changes. Results are returned per member with statuses `added`, `updated`,
`noop`, `not-member`, or `error`.

```bash
grok s groups add-members Chemists alice.mendel bob.curie
grok s groups add-members Chemists alice.mendel --admin        # promote / add as admin
grok s groups add-members Admins analysts                      # nest a group inside a group
grok s groups add-members Chemists alice.mendel --user         # force personal-group lookup
grok s groups remove-members Chemists alice.mendel
grok s groups list-members Chemists                            # all members
grok s groups list-members Chemists --admin                    # admin members only
grok s groups list-members Chemists --no-admin                 # non-admin members only
grok s groups list-memberships alice.mendel                    # groups a user belongs to
grok s groups list-memberships admin --user                    # the personal group, not "Administrators"
```

Names resolve by substring, an exact `name`/`friendlyName` match wins. When a name is still
ambiguous the command prints every matching group and exits non-zero — pass a UUID to
disambiguate, or add `--user` (accepted by all four membership verbs) to restrict lookup to
personal groups.

## User administration

```bash
grok s users block   alice.mendel        # accept login, UUID, or namespace:name
grok s users unblock alice.mendel
```

Both commands resolve the argument to a full user record first, so the server receives
the entire object (matching the Python client's `grok.users.block(user)`). Blocking
flips `status` to `Blocked` and terminates active sessions; unblocking restores `Active`.

## Sharing entities

```bash
grok s shares add "MyUser:MyConnection" Chemists,Biologists --access Edit
grok s shares list "MyUser:MyConnection"
```

The entity argument accepts either a UUID or an `"author:name"` pair. `--access`
defaults to `View`. `shares list` prints one row per group and permission, including the
grants inherited through the entity's project links (`inherited: true`).

## Running functions

```bash
grok s functions run 'Chem:smilesToMw("CCO")'                  # positional args
grok s functions run 'Pkg:fn({smiles:"CCO", radius:2})'        # named args
grok s functions run Pkg:fn --json params.json                 # big input from a file
```

The server binds arguments by parameter name, so positional values are mapped onto the
function's inputs in declared order (one extra lookup of the function); more values than inputs
is an error. Named arguments are sent as they are.

## Listing functions with filters

`grok s functions list` composes a smart filter from flags so you don't have to hand-roll
`text=` expressions. All flags are optional and combine with `and`; `--filter` lets you
add any other smart-filter clause.

```bash
grok s functions list --type script --language python          # all python scripts
grok s functions list --type query --package Chembl            # data-queries in one package
grok s functions list --package PowerPack --type package       # package-bundled funcs only
grok s functions list --language r --limit 20 --output json
grok s functions list --type script --filter 'name contains "similarity"'
```

Flag cheat sheet:

| Flag         | Values                                                     | Notes                                                             |
|--------------|------------------------------------------------------------|-------------------------------------------------------------------|
| `--type`     | `script`, `query` (alias for `data-query`), `function`, `package` | `function` matches both standalone and package-bundled funcs     |
| `--language` | `python`, `r`, `julia`, `nodejs`, `octave`, `grok`, `javascript`, `pyodide` | Filters in the CLI (the server doesn't index Script.language); implies `--type script` when used alone |
| `--package`  | package short name (e.g. `Chem`, `PowerPack`)              | Maps to `package.shortName` in the smart filter                  |
| `--filter`   | any smart-filter expression                                | Combined with the other flags via `and`                          |
| `--limit` / `--offset` | integers                                         | Applied client-side — the public list endpoint returns the full match set |

## Files

Remote paths are `<connector>/<file-path>`, where `<connector>` is the connection's
full name — including any namespace — e.g. `System:DemoFiles/smiles.csv`. The
separator between connector and path is the first `/`; the connector part may
contain colons (the namespace separator).

```bash
grok s files list "System:AppData" -r                          # recursive
grok s files list "System:AppData/MyPlugin"                    # path, kind (file|dir), size, updatedOn
grok s files list "System:DemoFiles" --output quiet            # paths only, pipe-friendly
grok s files get  "System:AppData/MyPlugin/config.json"
grok s files put  ./smiles.csv "System:DemoFiles/smiles.csv"   # upload local file
grok s files delete "System:AppData/MyPlugin/old.csv"
```

`files put` streams the file as raw bytes (no base64), so it handles GB-scale uploads
without blowing up memory. Use `batch files.put` only when you want to bundle an
upload with other operations in a single round-trip (the batch path base64-encodes
the `source` before sending).

## Tables

Table I/O — the shell counterpart to Python's `grok.tables.upload/download`.
Unlike `files put`, this registers a proper Datagrok table entity (returns `{ID, ...}`).

```bash
grok s tables upload MyTable ./data.csv                         # CSV → table
grok s tables upload MyTable ./data.d42                         # d42 binary → table
grok s tables upload MyTable ./data.csv --output json           # get ID and markup back
grok s tables download MyTable                                  # CSV to stdout (pipe-friendly)
grok s tables download MyTable -O ./data.csv                    # CSV to a local file
grok s tables download "Admin:MyTable:MyTable"                  # full name, or the UUID
grok s tables list --filter MyTable                             # registered tables
grok s tables get MyTable                                       # the TableInfo record
grok s tables delete MyTable
```

A table lives under a project namespace (the upload prints its full name,
`Admin:MyTable:MyTable`); the bare name works while only one table carries it, otherwise the
CLI lists the candidates and asks for the full name or the id.

Upload streams raw bytes — `Content-Type: text/csv` for `.csv` and
`application/octet-stream` for `.d42` (auto-detected from the file extension). Both
formats handle large tables without loading the whole file into a JSON envelope.
`-O` / `--output-file` avoids colliding with the format flag (`--output
table|json|csv|quiet`), which still controls how the upload result is printed.

Download is CSV-only — the server reads the stored d42 blob and converts. If you
need the raw d42 bytes, hit `/tables/data/<id>` directly via `grok s raw`.

## Packages

Manage server plugins the same way the Package Manager UI does — the server pulls
published packages from its configured package repository (npm). This is for
installing **released** packages on a running server; publishing your own local
package sources is still `grok publish`.

```bash
grok s packages install Chem Bio PowerGrid       # install latest of each, one line
grok s packages install Chem --version 1.14.0    # pin a specific version (single package only)
grok s packages outdated                         # installed vs registry-latest
grok s packages update Chem Bio                  # upgrade named packages to latest
grok s packages update --all                     # upgrade everything outdated
grok s packages versions Chem                    # published versions w/ current/latest/debug flags
grok s packages set-version Chem 1.13.0          # activate a specific version (pulls it if needed)
grok s packages uninstall Chem                   # remove; repository entry stays installable
grok s packages share Chem Chemists --access View
```

Semantics worth knowing:

- **`install` without `--version` sets the package to `latest`**, which is an
  auto-update intent: the server re-resolves it against the registry every ~15
  minutes, so the package keeps itself current. A pinned `--version` (and `update`,
  which pins the concrete latest version) stays put until you change it.
- **Install is synchronous** — the call returns once the version is downloaded,
  published, and made current. A heavy first-time install (e.g. Chem) can take
  minutes; if the HTTP call times out, just re-run it — install is idempotent, and
  a version that finished installing server-side is activated instantly on retry.
- Multi-package `install`/`update` runs sequentially and prints a per-package
  status table (`installed` / `noop` / `error`); a partial failure sets exit code 1.
- **`update` preserves auto-update tracking**: a package whose desired version is
  `latest` is updated by re-resolving `latest` (it keeps auto-updating); a pinned
  package gets pinned to the concrete registry-latest version. Note that a
  locally-published dev build counts as outdated once the registry moves past its
  base version, so `update --all` will replace it with the registry release.
- Status semantics: `noop` means "already at the requested latest"; installing an
  explicit `--version` reports `installed` (with the published-version id) even
  when that version was already current — the server re-activates it.
- **`uninstall` of a repository-backed package keeps the package entry** (it shows
  up as installable again); only locally-published packages are deleted outright.
  It also clears the desired version, so the package drops out of `outdated` until
  reinstalled. Neither path cleans up package credentials or package DB schemas.
- `versions`, `outdated`, and name resolution read the repository catalog, so they
  need a configured package repository — check with `grok s raw GET /api/packages/repos`.
- `share` targets the package's project under the hood (same as sharing from the UI),
  so the whole package — functions, queries, connections — is shared at once. For the
  same reason `grok s shares list <package-uuid>` won't show the grant — it lives on
  the package's project, not the package entity itself.

## Domain schemas and rows

Entity-mapped domain tables — the schemas plugins declare in `databases/<schema>/schema.json`
and the user-managed ones created at runtime — are reachable through one entity, `domains`.
Every verb takes an address: a bare `<schema>` names a schema, `<schema>.<table>` names a
table, mirroring `grok.dapi.domains.schema('grit')` and `grok.dapi.domains.table('grit.issue')`.
Reads return only the rows and columns the key's user may see; writes are validated,
permission-checked and audited by the server exactly as they are from the UI.

### Browsing

```bash
grok s domains list                                 # schemas: name, managedBy, version, table count
grok s domains list grit                            # tables of one schema: security mode, business key, ...
grok s domains get grit                             # the manifest, as JSON (doubles as an export)
grok s domains get grit.issue                       # the table's columns (--output json: its manifest section)
grok s domains get grit.issue <row-id>              # one row
grok s domains access grit.issue                    # can.view/insert/edit/delete/share + editable/readonly column lists
```

### Querying

```bash
grok s domains query grit.issue --filter 'status = "open"' --sort '!created_on' --limit 20
grok s domains query grit.issue --columns title,status --expand project_id --offset 100
grok s domains count grit.issue --filter 'status = "open"'
grok s domains aggregate grit.issue --measures 'count,avg(estimate) as mean' --group-by status
grok s domains aggregate grit.issue --json spec.json             # any DomainAggregateSpec
grok s domains download grit.issue -O ./issues.csv --filter 'status = "open"'
grok s domains download grit.issue -O ./issues.d42                # typed DataFrame, 10M-row cap
```

`--filter` is the domain smart-filter grammar (`status = "open"`, `title contains "crash"`,
`quantity > 10`); values are bound server-side. `--sort` is a comma list with `!` for
descending; `--expand` takes `<fk_column>`, `details:<child>` or a relation name. JSON output
and CSV downloads are capped at 10k rows by the server; a `.d42` download uses the DataFrame
path and goes up to 10M. `--limit` defaults to 50 for `query` and is unbounded for `download`.

### Writing rows

```bash
grok s domains insert grit.issue title="Crash on save" status=open project_id=<uuid>
grok s domains insert grit.issue --json rows.json                # one object or an array
grok s domains update grit.issue <row-id> status=closed --version 3   # optimistic concurrency
grok s domains delete grit.issue <row-id>
grok s domains delete grit.issue --filter 'status = "closed"' --limit 500   # bulk, oldest first
grok s domains transaction grit --json ops.json                  # ordered ops, one transaction
```

Inline `col=value` pairs are typed when the value parses as JSON (`quantity=5`,
`done=true`, `note=null`, `tags=["a","b"]`) and sent as strings otherwise. A business-key
duplicate is reported as `status: duplicate` with the existing id; add `--error-on-duplicate`
to fail instead. A validation failure prints one line per offending column and exits 1.
`delete --filter` removes at most 1000 rows per call and says when more remain.
`ops.json` is the `DomainsDataSource.transaction` ops list — `{op, table, ref?, values?, id?,
expectedVersion?}`, with `"$ref"` placeholders for earlier ops' ids.

### Bulk upload

```bash
grok s domains upload grit.issue ./issues.csv                    # insert
grok s domains upload grit.issue ./issues.csv --upsert           # merge by business key
grok s domains upload grit.issue ./issues.d42                    # d42 DataFrame
grok s domains upload grit.issue ./issues.json --no-all-or-nothing --error-on-duplicate
```

The format follows the extension (`.csv`, `.d42`, `.json` — a row array, bare or under
`rows`). The whole batch is one transaction by default; `--no-all-or-nothing` applies the good
rows and reports the bad ones per row. The report prints `inserted / updated / skipped /
errors` plus a row table for the failures; any error sets exit code 1.

### Schema lifecycle (user-managed schemas)

```bash
grok s domains create inventory --friendly-name "Inventory" --description "Lab stock"
grok s domains apply inventory --json schema.json --dry-run      # the change plan, no writes
grok s domains apply inventory --json schema.json                # create / alter tables
grok s domains apply inventory --json schema.json --confirm-destructive --if-version 3
grok s domains audit inventory --limit 50                        # DDL and row events, newest first
grok s domains delete inventory --force                          # purge: data, audit, registry, grants
```

`schema.json` may be a full manifest or a partial apply body — only `tables`, `extend`,
`propertySchemas` and `dropTables` are sent; `name`, `version` and `description` are dropped.
Named tables replace their definition wholesale; untouched tables stay as the registry has
them. A plan that drops or narrows anything is refused until `--confirm-destructive` is
passed, and the plan is printed with the refusal. `--if-version` fails the apply when the
schema's apply counter (or `ext_version` on a package schema) has moved. Package-deployed
schemas cannot be applied to or deleted here — `apply` on one is the user-extension path
(needs `Extend`), and the manifest is owned by `grok publish`. `delete <schema>` requires
`--force` because it takes every row with it.

### Grants

```bash
grok s domains grants grit.issue                                 # direct permission rows
grok s domains grant grit.issue Chemists,Biologists --access Edit
grok s domains revoke grit.issue Chemists --access Edit          # omit --access to revoke all
grok s domains grant grit Chemists --access Extend               # schema-level: may add own tables/columns
```

Groups resolve by name the same way `groups add-members` does (a login resolves to the
personal group; ambiguous names list the candidates). Table grants gate row data; schema
grants gate schema operations (`Edit` to apply, `Delete` to purge, `Share`, `Extend`) and do
not reach rows. A grant on a table restricted by per-column sharing does not un-hide the
column — that path (`/domains/grants/column`) is `grok s raw` territory for now.

### Not covered here

Watch/subscribe, facets, saved filters, per-column sharing and row promotion have no verb;
`grok s raw <METHOD> /api/domains/...` reaches them. Domain schemas and their data do not
take part in `pull` / `push` / `migrate` — install the package (or re-`apply` the manifest)
on the target and `upload` the rows.

## Server health

```bash
grok s healthcheck                             # full per-module health
grok s healthcheck --module scripting          # filter to one module
grok s healthcheck --output json               # machine-readable
```

Hits `GET /public/v1/healthcheck`. Response:

```json
{
  "status":  "ok",
  "server":  "https://public.datagrok.ai",
  "version": "1.27",
  "time":    "2026-04-19T19:20:00.000Z",
  "services": [
    { "key": "scripting", "type": "Service", "name": "...", "status": "Running", "started": true, "enabled": true, "time": "..." }
  ]
}
```

`services` is the same payload as `/admin/health` (per-`GrokServiceInfo` records). A server
that reports no services (a dev stack) prints `(no services reported)`; `--module` for a
module the server does not report exits 1. Requires a valid dev key (standard `grok s` auth).
For an anonymous liveness probe — load balancer, k8s readiness — hit `/admin/health` directly;
it's on the server's unauthenticated allowlist.

## Describing an entity type

```bash
grok s describe connections                    # fields of a DataConnection, with types and examples
grok s describe Project --output json          # by type name; JSON carries the registry record and a sample
grok s describe users --output quiet           # field names only
```

The server publishes no JSON schema, so `describe` combines the entity-type registry record
with the top-level fields of one existing entity of that type (`field`, `type`, `example`).
Aliases: `users groups connections queries scripts functions packages reports tables projects
files`; any other registered type name (`Project`, `ViewLayout`) works when at least one
entity exists.

## Raw API access

When no dedicated subcommand exists, fall through to `grok s raw`:

```bash
grok s raw GET  /users/current
grok s raw GET  /packages/dev/MyPlugin
grok s raw POST /admin/reload-settings
grok s raw POST /public/v1/functions/Sin/call --data '{"x": 1}'
grok s raw POST /domains/grants/<id> --json grant.json
```

Paths are relative to the API base of the target server (`/users/current` becomes
`https://host/api/users/current`, or `http://host:8082/users/current` on a bare Datlas); a
leading `/api` is accepted and dropped, so old `/api/...` paths keep working on both host
shapes. A body comes from `--json <file>` or `--data '<json>'`. A non-2xx answer, or a
200 carrying an `ApiError`, prints the message with the HTTP status on stderr and exits 1;
under `--output json` the stderr line is `{"error", "errorCode", "body"}`.

On **Windows Git Bash**, prefix raw paths with `MSYS_NO_PATHCONV=1` to stop the shell
from rewriting POSIX paths into Windows paths:

```bash
MSYS_NO_PATHCONV=1 grok s raw GET /users/current
```

## Batch operations

Apply the same verb to many items in one round-trip:

```bash
# Inline args
grok s batch files delete "System:AppData/old1.txt" "System:AppData/old2.txt"

# From a JSON array
grok s batch users save --json users.json      # [{...user1}, {...user2}, ...]

# Full workflow manifest — mixed actions, optional transaction / stopOnError
grok s batch manifest.json
```

Manifest shape:

```json
{
  "stopOnError": true,
  "transaction": false,
  "operations": [
    { "id": "op1", "action": "users.create",  "params": {"login": "alice", "firstName": "Alice"} },
    { "id": "op2", "action": "groups.create", "params": {"name": "Chemists"} }
  ]
}
```

Actions the server accepts: `create | get | delete` for `users`, `groups`, `connections`,
`functions`, `queries`, `scripts` (`get | delete` for `reports`), `functions.run`
(`{name, params}`), and `files.list | get | put | delete`. For `files.put`, add
`"source": "<local-path>"` and the CLI base64-encodes the file into `content` before sending.
`users.delete` removes the entity record only (see "List / count / get / delete").

## Scripting pattern

`--output quiet` and `--output json` make `grok s` safe to compose with standard
shell tooling:

```bash
# Pipe IDs
grok s users list --filter "status = 'active'" --output quiet \
  | xargs -I{} grok s users get {}

# Filter with jq
grok s connections list --output json \
  | jq '.[] | select(.dataSource=="Postgres") | .name'

# Generate JSON on the fly
for login in alice bob carol; do
  printf '{"#type":"User","login":"%s","firstName":"%s","status":"active"}\n' \
    "$login" "${login^}" > user.json
  grok s users save --json user.json
done
```

## Worked example: seed users and populate groups

Create a group, bulk-create users, and drop each into the right group without touching
the UI:

```bash
# 1. Create the group (idempotent if you include the existing id)
cat > /tmp/g.json <<'EOF'
{ "#type": "UserGroup", "name": "Chemists", "friendlyName": "Chemists" }
EOF
grok s groups save --json /tmp/g.json

# 2. Create 8 users
for row in \
  "alice.mendeleev:Alice:Mendeleev" \
  "bob.curie:Bob:Curie" \
  "carol.pauling:Carol:Pauling" ; do
  IFS=: read -r login first last <<<"$row"
  cat > /tmp/u.json <<EOF
{ "#type": "User", "login": "$login", "firstName": "$first", "lastName": "$last", "status": "active" }
EOF
  grok s users save --json /tmp/u.json --output quiet
done

# 3. Add them all to the group in one call (resolves logins to personal groups)
grok s groups add-members Chemists alice.mendeleev bob.curie carol.pauling --user

# 4. Verify
grok s groups list-members Chemists --no-admin
```

Steps 1 and 2 create; re-running them fails with "already exists" (exit 1) unless the JSON
carries the existing `id`. Step 3 is idempotent (`noop` on a re-run).

## Migrating entities between instances (pull / push / migrate)

Entities built in the UI on one instance — connections, queries, scripts, dashboards,
spaces, layouts, tables, files, jobs, notebooks, models, and the groups and grants they
need — are promoted to another instance through a **bundle**: a directory of one JSON file
per entity that travels by any means (a commit, a PR, a USB stick).

```bash
grok s pull   <selection> --out ./bundle --host dev      # instance → directory
grok s bundle ls ./bundle                                # what is in it
grok s diff   ./bundle --host prod                       # what a push would change
grok s push   ./bundle --host prod [--dry-run]           # directory → instance
grok s migrate <selection> --from dev --to prod          # pull + push, temp dir in between
```

Every entity keeps the **same UUID** on both instances, so a push is idempotent: pushing an
unchanged bundle a second time writes nothing, and 1.28's built-in server-to-server sync
recognises what the CLI pushed as its own.

### Bundle layout

```
bundle/
  manifest.json                     # source url + version, pull history, FK-safe order
  DataConnection/Chem.Chembl.json   # <Type>/<nqName with ':' and '/' as '.'>.json
  DataQuery/Chem.CompoundsByTarget.json
  Project/Chem.Dashboard.json
  UserGroup/Chemists.json           # bare group + member logins, never member ids
  FileInfo/reports.readme.md.json   # a file without a namespace is named by its path
  tables/<id>.d42                   # table data (on by default for pulled tables)
  files/<id>                        # file bytes (only with --include-files)
  idmap.json                        # sourceId -> targetId, written by --on-conflict adopt and
                                    # whenever a FileInfo save answers with an existing row's id
```

**Pulls accumulate.** Pulling into an existing bundle merges: entities already there are
overwritten by id, new ones are added, nothing is removed, and `manifest.pulls[]` keeps one
entry per invocation. `--replace` clears the directory first.

```bash
grok s pull Chem:TargetDashboard --out ./release --host dev
grok s pull --type script --author alice --no-deps --out ./release --host dev
grok s pull Chemists --out ./release --host dev
grok s push ./release --host prod
```

### Selection

| Flag | Selects |
|---|---|
| positional `Chem:Dashboard <uuid> ...` | exactly these entities, by nqName or id |
| `--type conn,query,script,project,dashboard,space,view,layout,table,file,group,job,notebook,model` | which types to list (default: all of them) |
| `--name <glob>` | free-text search, then a client-side glob on name and friendly name (`Cereal*`, `*demo*`) |
| `--namespace Chem` | everything **under** the namespace, recursively — the space `Chem` itself is not one of them |
| `--space Chem:Reports` | the space itself **plus** everything under it |
| `--author alice` | entities authored by a login |
| `--tag demo` | entities carrying a tag (types whose router has no `tags` param are skipped with a warning) |
| `--since 2w` / `--since=-30d` / `--since 2026-08-01` | updated since (bare `2w` means `-2w`; the shell eats a leading `-` unless you use `=`) |
| `--filter "<expr>"` | a smart-filter expression, ANDed with the rest |
| `--no-deps` | do not follow dependencies (what never travels is still excluded — see below) |
| `--no-include-data` | do not pull the `.d42` data of the tables that were pulled (data is on by default) |
| `--include-files` | also pull the bytes of the files that were pulled (off by default) |
| `--verbose` | print the stack of a runtime failure instead of the one-line message |

`--type space` and `--type dashboard` are both `Project` with a different listing rule, so
they cannot be combined in one command (neither can `--type project,space`).

### What the walker adds, and what never travels

A selected entity brings its dependencies with it: a project brings its relation children
(recursively), their views, layouts and tables; a query brings its connection; a datasync
table brings the connection its creation script opens; a notebook and a model bring their
tables; every entity brings the non-personal groups that hold grants on it, and a group
brings its parents and members. A **job** brings nothing — 1.27 does not persist the link
between a job and the queries it runs, so pull those explicitly.

The table below holds however the entity was chosen: `--no-deps` skips the walk, not these
rules, and `push` refuses the same connections even if a bundle was hand-edited to carry one
(`skip(platform_connection)` / `skip(personal_storage)` / `skip(space_files_connection)`).

| Never travels | Row on the report |
|---|---|
| Passwords and other password-class connection parameters | `needs-credentials` (see `--creds`) |
| Users — memberships are replayed by login and group name | `warn(member_not_found)` when a login is not on the target |
| Package-owned entities — install the package instead | `warn(package_entity)`, `warn(package_not_installed)` |
| `System:` connections, the personal `Home` share, a space's own `Files` connection | `info(platform_connection)` / `warn(personal_storage)` / `info(space_files_connection)` |
| A trained model blob | `info(model_blob_skipped)` |

### Credentials: `--creds`

A pushed connection arrives without its secrets. Author a YAML file **for the target** and
pass it to `push` or `migrate`; the values are merged into the connection's parameters before
the save, and the server encrypts and masks them itself. Nothing is ever read back into the
bundle.

```yaml
# creds.yaml — keys are connection nqNames as the bundle spells them
Chem:Chembl:
  password: ${CHEMBL_PROD_PASSWORD}
Admin:Northwind:
  password: ${NORTHWIND_PASSWORD}
```

```bash
CHEMBL_PROD_PASSWORD=... grok s push ./release --host prod --creds ./creds.yaml
```

`${VAR}` is resolved from the environment exactly as `grok publish` resolves it in
`connections/*.json` (the file is parsed as YAML first, so a `${VAR}` written in flow style
needs quoting: `{password: "${VAR}"}`); a variable that is not set aborts the push before the
first write. A connection the file covers gets no `needs-credentials` row, and is **always**
written — a secret is invisible in the payload, so a connection that would otherwise be
`identical` is planned as `update` with reason `credentials`. That is how a password is
rotated: re-run the push with a new value.

### Conflicts

An entity whose id is absent on the target but whose name is taken by another id is a
conflict. `--on-conflict` decides:

| Policy | What happens |
|---|---|
| `fail` (default) | nothing is written; every conflict is listed and the command exits 1 |
| `skip` | the twin is left alone; anything in the bundle that points at it is reported `failed(dependency_skipped)` |
| `adopt` | the bundle entity is written **into** the twin, and `idmap.json` records `sourceId -> targetId` so every later reference and every later push follows it |
| `duplicate` | the bundle entity is created under its own id next to the twin (the server renames it `Name_1`). A `UserGroup` cannot be duplicated — group names are unique on 1.27, so that row fails |

### Reading the plan

`diff` and `--dry-run` print the plan without writing (`diff` plans with `skip`, so a
conflict does not abort it). Actions: `create`, `update` (with the changed top-level keys as
the detail), `identical` (no write), `skip`, `failed`, plus the `warn` / `info` /
`needs-credentials` notes. `--output json` emits the 1.28-compatible shape — the same one for
`diff`, `push` and `migrate`, with `detail` always present (`""` when there is nothing to
say):

```json
{
  "items": [{"name": "Chem:Chembl", "entityType": "DataConnection", "action": "create", "reason": "", "detail": ""}],
  "counts": {"create": 1, "identical": 4},
  "status": "ok",
  "remoteUrl": "https://prod.datagrok.ai/api"
}
```

A push exits 1 if any row is `failed`.

### 1.27 limits worth knowing

- A **FileInfo that lives in a share** is a dead row on the target — only stand-alone blobs
  migrate (`info(file_in_share_not_migratable)`), the same rule the 1.28 sync applies.
- `metaParams` do not survive a save on connections and queries, so metadata attached there
  does not travel.
- `GET /projects/relations` fails on a project that links domain-table rows; the walker falls
  back to the project's own `relations[]` and reports `warn(relations_degraded)`.
- Relations are **merged**, never replaced: a link that exists only on the target is kept, so
  removing a relation from a bundle never unlinks it there (`info(relation_not_removed)`).
  Unlink it in the UI on the target instead.
- A space delete leaves its `…:Files` connection behind, and it cannot be deleted through the
  API. Deleting a group that still holds grants fails — revoke the grants first.

### dev to prod, end to end

```bash
# 1. see what would come over, from dev, without touching prod
grok s migrate Chem:TargetDashboard --from dev --to prod --dry-run

# 2. promote it, filling in the target's secrets
CHEMBL_PROD_PASSWORD=... grok s migrate Chem:TargetDashboard \
  --from dev --to prod --creds ./creds.yaml --keep

# 3. re-run: everything is `identical`, nothing is written
grok s migrate Chem:TargetDashboard --from dev --to prod

# 4. or keep the bundle under version control instead
grok s pull Chem:TargetDashboard --out ./release --host dev
git add release && git commit -m "release: target dashboard"
grok s push ./release --host prod --creds ./creds.yaml
```

`migrate` pulls into a temporary directory and deletes it afterwards; `--keep` keeps it and
prints its path on **stderr**, so `--output json` stays one parseable document. `--from` is
only ever read from.

## Implementation notes

- Source: `public/tools/bin/commands/server.ts`, `public/tools/bin/utils/node-dapi.ts`;
  pull / push / migrate in `bin/commands/server-migrate.ts` + `bin/utils/migrate/`.
  Domain schemas and rows in `bin/commands/server-domains.ts` (`NodeDomainsDataSource` in
  `node-dapi.ts`); they call the internal `/domains/` router the browser uses.
- The Node client talks directly to `/public/v1/` — no Dart interop, no browser, no
  logged-in session required. Authentication uses the developer key from the config.
- If `grok s` is not working, start by running `grok s healthcheck` — it verifies the
  URL, the key, and basic connectivity, and returns per-module status if the server is
  reachable. Fall back to `grok s raw GET /users/current` to isolate auth issues. A `--host`
  URL that is not the API base fails at login with the reason (`should end with /api`).
- Cross-instance sync (1.28 servers): `grok s sync pairs list`, `sync setups list --pair <id>`,
  `sync setup get <id>`, `sync run <id>`; on 1.27 the routes do not exist and the commands
  answer "not found".
