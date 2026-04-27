---
title: "Server management with grok s"
sidebar_position: 6
keywords: [grok, cli, server management, automation, integration, ldap, active directory]
---

`grok s` (alias of `grok server`) is a command-line tool for managing a running
Datagrok server over its [public REST API](packages/rest-api.md). Use it for
any task that creates, updates, shares, or inspects entities on the server —
users, groups, connections, queries, scripts, packages, files, function calls —
and for bulk operations or integrations that should not depend on a browser.

It is the recommended alternative to UI automation (Playwright, Selenium,
hand-written HTTP) for data and configuration management: faster, idempotent,
scriptable, and free of any logged-in session.

`grok s` ships with [`datagrok-tools`](https://www.npmjs.com/package/datagrok-tools).
Install with `npm install -g datagrok-tools`, then configure a server.

## When to use grok s

| Task                                                       | Command                                                    |
|------------------------------------------------------------|------------------------------------------------------------|
| Create or update a user                                    | `grok s users save --json user.json`                       |
| Block or unblock a user                                    | `grok s users block <login>` / `users unblock <login>`     |
| Create or update a group                                   | `grok s groups save --json group.json`                     |
| Add or remove members in a group                           | `grok s groups add-members <group> <user>...`              |
| List members of a group                                    | `grok s groups list-members <group>`                       |
| List the groups a user belongs to                          | `grok s groups list-memberships <user>`                    |
| Share a connection, query, script, or project with a group | `grok s shares add <entity> <group> [--access View\|Edit]` |
| Create or update a connection                              | `grok s connections save --json conn.json`                 |
| Test a connection                                          | `grok s connections test <id-or-name>`                     |
| Run a registered function                                  | `grok s functions run 'Pkg:fn(arg1,arg2)'`                 |
| List or browse files in a file share                       | `grok s files list "System:AppData" -r`                    |
| Upload or download a table (CSV)                           | `grok s tables upload <name> file.csv`                     |
| Check whether a package is deployed                        | `grok s packages list --filter "MyPlugin"`                 |
| Hit any undocumented endpoint                              | `grok s raw GET /api/users/current`                        |
| Check server and per-module health                         | `grok s healthcheck [--module <name>]`                     |
| Bulk operations in one round-trip                          | `grok s batch <entity> <verb> --json items.json`           |

## Configuration

Servers and developer keys live in `~/.grok/config.yaml`:

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

Add a new entry with `grok config --server --alias <name> --server <url> --key <key>`,
or pass `--host <alias-or-url>` to any `grok s` command to override the default.
The developer key is the one shown on your **User Profile** page next to **API Key**.

## Common workflows

### Manage users and groups

```bash
grok s users list                                    # default: table, 50 rows
grok s users list --filter "status = 'active'"       # smart filter
grok s users save --json user.json                   # create or update
grok s users block alice.mendel                      # accept login, UUID, or namespace:name

grok s groups save --json group.json
grok s groups add-members Chemists alice.mendel bob.curie
grok s groups add-members Admins analysts            # nest a group inside a group
grok s groups list-members Chemists --no-admin
grok s groups list-memberships alice.mendel
```

`add-members` is idempotent: each member is resolved (by login, name, or UUID),
compared against the group's current children, and only written when something
changes. Pass `--user` to disambiguate a name as a personal group.

### Share entities

```bash
grok s shares add "MyUser:MyConnection" Chemists,Biologists --access Edit
grok s shares list <entity-uuid>
```

The entity argument accepts either a UUID or an `"author:name"` pair.
`--access` defaults to `View`.

### Run functions

```bash
grok s functions run 'Chem:smilesToMw("CCO")'                  # positional args
grok s functions run 'Pkg:fn({smiles:"CCO", radius:2})'        # named args
grok s functions run Pkg:fn --json params.json                 # big input from a file
```

### Browse files and tables

```bash
grok s files list "System:AppData" -r
grok s files put ./smiles.csv "System:DemoFiles/smiles.csv"

grok s tables upload MyTable ./data.csv
grok s tables download MyTable -O ./data.csv
```

`files put` streams raw bytes — it handles GB-scale uploads. `tables upload`
registers a proper Datagrok [table entity](../datagrok/concepts/objects.md)
and returns its ID.

### Server health

```bash
grok s healthcheck                             # full per-module health
grok s healthcheck --module scripting          # filter to one module
grok s healthcheck --output json               # machine-readable
```

Hits `GET /public/v1/healthcheck`. For an anonymous liveness probe (load
balancer, Kubernetes readiness), use `/admin/health` directly.

### Raw API access

```bash
grok s raw GET  /api/users/current
grok s raw POST /api/admin/reload-settings
```

On Windows Git Bash, prefix with `MSYS_NO_PATHCONV=1` so the shell does not
rewrite POSIX paths.

## Scripting patterns

`--output quiet` and `--output json` make `grok s` safe to compose with shell
tooling:

```bash
grok s users list --filter "status = 'active'" --output quiet \
  | xargs -I{} grok s users get {}

grok s connections list --output json \
  | jq '.[] | select(.dataSource=="Postgres") | .name'
```

Every subcommand on this page is idempotent — re-running a script that already
ran is safe.

## Sync an AD group with Datagrok

A common integration scenario: keep a Datagrok group in sync with the membership
of an Active Directory group. Pull the AD members from your directory of choice,
then drive Datagrok with `grok s`:

```bash
# 1. Make sure the Datagrok-side group exists.
cat > /tmp/g.json <<'EOF'
{ "#type": "UserGroup", "name": "Chemists", "friendlyName": "Chemists" }
EOF
grok s groups save --json /tmp/g.json

# 2. Create any users that don't exist yet (one user per logon-name from AD).
for row in alice.mendeleev:Alice:Mendeleev bob.curie:Bob:Curie ; do
  IFS=: read -r login first last <<<"$row"
  cat > /tmp/u.json <<EOF
{ "#type": "User", "login": "$login", "firstName": "$first", "lastName": "$last", "status": "active" }
EOF
  grok s users save --json /tmp/u.json --output quiet
done

# 3. Reconcile membership in one call (idempotent: existing members are noop).
grok s groups add-members Chemists alice.mendeleev bob.curie --user

# 4. Verify.
grok s groups list-members Chemists --no-admin
```

Combine with the [LDAP authentication setup](../deploy/complete-setup/configure-auth.md#ldap-authentication)
to let the same users sign in with their AD credentials.

## Full reference

The page above covers the most common tasks. For the exhaustive command
reference — JSON shapes, batch manifests, `connections save --save-credentials`,
all output flags — see
[`tools/GROK_S.md`](https://github.com/datagrok-ai/public/blob/master/tools/GROK_S.md)
in the public repository.

## See also

* [REST API](packages/rest-api.md) — the underlying HTTP surface
* [JavaScript API](packages/js-api.md) — for in-platform plugin code
* [Users and groups](../govern/access-control/users-and-groups.md)
* [Configure authentication](../deploy/complete-setup/configure-auth.md)
* [Manage credentials](how-to/packages/manage-credentials.md)
