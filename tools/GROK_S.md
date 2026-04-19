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
| See who an entity is shared with                            | `grok s shares list <entity-id>`                       |
| Create / update a connection                                | `grok s connections save --json conn.json`             |
| Test a connection                                           | `grok s connections test <id-or-name>`                 |
| Run a registered function                                   | `grok s functions run 'Pkg:fn(arg1,arg2)'`             |
| List functions with rich filters                            | `grok s functions list --type script --language python --package Chem` |
| List / browse files in a file share                         | `grok s files list "System:AppData" -r`                |
| Upload / download a table (CSV)                             | `grok s tables upload <name> file.csv` / `tables download <name> -O out.csv` |
| Check whether a package is deployed                         | `grok s packages list --filter "MyPlugin"`             |
| Hit any undocumented endpoint                               | `grok s raw GET /api/users/current`                    |
| Check server + per-module health                            | `grok s healthcheck [--module <name>]`                 |
| Bulk operations in one round-trip                           | `grok s batch <entity> <verb> --json items.json`       |

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

- `grok config --server --alias <name> --server <url> --key <key>` writes a new entry.
- Add `--default` to make it the active server.
- Every `grok s ...` command accepts `--host <alias-or-url>` to override the default.

## Entity operations

### List / get / delete

```bash
grok s users list                              # default: table, 50 rows
grok s users list --filter "login = 'admin'"   # smart filter
grok s connections list --output json
grok s packages list --filter "name:MyPlugin"  # table shows friendlyName
grok s groups get <id-or-name>
grok s connections delete <id>
```

Options that work on every entity:

| Flag                  | Meaning                                           |
|-----------------------|---------------------------------------------------|
| `--output table\|json\|csv\|quiet` | Output format; `quiet` prints ids only |
| `--filter "<expr>"`   | Server-side smart filter                          |
| `--limit <n>`         | Page size (default 50)                            |
| `--offset <n>`        | Skip first `n` rows                               |
| `--host <alias\|url>` | Target a specific server from your config         |

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
```

When a name is ambiguous the command prints every matching group and exits non-zero —
pass a UUID to disambiguate, or add `--user` to restrict lookup to personal groups.

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
grok s shares list <entity-uuid>
```

The entity argument accepts either a UUID or an `"author:name"` pair. `--access`
defaults to `View`.

## Running functions

```bash
grok s functions run 'Chem:smilesToMw("CCO")'                  # positional args
grok s functions run 'Pkg:fn({smiles:"CCO", radius:2})'        # named args
grok s functions run Pkg:fn --json params.json                 # big input from a file
```

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
grok s files list "System:AppData/MyPlugin"
grok s files get  "System:AppData/MyPlugin/config.json"
grok s files put  ./smiles.csv "System:DemoFiles/smiles.csv"   # upload local file
grok s files delete "System:AppData/MyPlugin/old.csv"
```

`files put` streams the file as raw bytes (no base64), so it handles GB-scale uploads
without blowing up memory. Use `batch files.put` only when you want to bundle an
upload with other operations in a single round-trip (the batch path base64-encodes
the `source` before sending).

## Tables

CSV-level table I/O — the shell counterpart to Python's `grok.tables.upload/download`.
Unlike `files put`, this registers a proper Datagrok table entity (returns `{ID, ...}`).

```bash
grok s tables upload MyTable ./data.csv                         # CSV → table
grok s tables upload MyTable ./data.csv --output json           # get ID and markup back
grok s tables download MyTable                                  # CSV to stdout (pipe-friendly)
grok s tables download MyTable -O ./data.csv                    # CSV to a local file
grok s tables download <uuid>                                   # UUID or namespace:name both work
```

Upload streams raw bytes with `Content-Type: text/csv` so it handles large tables
without loading the whole file into a JSON envelope. `-O` / `--output-file` avoids
colliding with the format flag (`--output table|json|csv|quiet`), which still controls
how the upload result is printed.

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

`services` is the same payload as `/admin/health` (per-`GrokServiceInfo` records).
Requires a valid dev key (standard `grok s` auth). For an anonymous liveness probe —
load balancer, k8s readiness — hit `/admin/health` directly; it's on the server's
unauthenticated allowlist.

## Raw API access

When no dedicated subcommand exists, fall through to `grok s raw`:

```bash
grok s raw GET  /api/users/current
grok s raw GET  /api/packages/dev/MyPlugin
grok s raw POST /api/admin/reload-settings
```

On **Windows Git Bash**, prefix raw paths with `MSYS_NO_PATHCONV=1` to stop the shell
from rewriting POSIX paths into Windows paths:

```bash
MSYS_NO_PATHCONV=1 grok s raw GET /api/users/current
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
    { "id": "op1", "action": "users.save",  "params": {"login": "alice", "firstName": "Alice"} },
    { "id": "op2", "action": "groups.save", "params": {"name": "Chemists"} }
  ]
}
```

For `files.put`, add `"source": "<local-path>"` and the CLI base64-encodes the file
into `content` before sending.

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

Every subcommand above is idempotent — re-running the whole block is safe.

## Implementation notes

- Source: `public/tools/bin/commands/server.ts`, `public/tools/bin/utils/node-dapi.ts`.
- The Node client talks directly to `/public/v1/` — no Dart interop, no browser, no
  logged-in session required. Authentication uses the developer key from the config.
- If `grok s` is not working, start by running `grok s healthcheck` — it verifies the
  URL, the key, and basic connectivity, and returns per-module status if the server is
  reachable. Fall back to `grok s raw GET /api/users/current` to isolate auth issues.
