---
name: db-in-docker
version: 0.3.0
description: |
  Ship a database engine (Postgres with custom extensions, MariaDB,
  ClickHouse, Oracle, ...) inside a sibling container alongside a Datagrok
  package and wire the platform to it. For plugin authors whose data layer
  needs an engine the shared RDS doesn't provide or DBA-level isolation.
  Produces `dockerfiles/Dockerfile` + `dockerfiles/container.json` +
  `connections/<name>.json` whose `server` placeholder resolves to the
  running container at deploy time, so queries against
  `<Pkg>:<Connection>` route transparently through it.
  Use when asked to "ship a bundled database with my plugin", "wire a
  connection to a sibling container database", or "run postgres with
  custom extensions alongside the package".
triggers:
  - bundled database alongside the package
  - sibling container database connection
  - postgres extensions not in shared rds
  - wire connection to containerized database
  - custom postgres engine for a plugin
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
harness-authored: true
---

# db-in-docker

## When to use

Pick this path when your package needs a database engine the shared
platform RDS won't run:

- a non-Postgres DBMS (MariaDB, ClickHouse, Oracle);
- a Postgres image with extensions outside the shared instance's
  allowlist;
- required DBA-level isolation from other tenants (`DG-FACT-418`).

For ordinary CRUD storage that fits in the shared Postgres, prefer
[[db-in-plugin]] — it's simpler, the platform manages migrations and
lifecycle, and it uses far fewer resources. If none of the three
bullets above apply, stop and use that path instead.

## Prerequisites

- Package scaffold (`grok create <Name>`); paths below are relative to
  the package root. Use [[create-package]] if absent.
- API key from `<GROK_HOST>/u` for the credentials POST in step 5.
- Familiarity with the connection-and-query spine from [[access-data]] —
  this skill only replaces step 1 of that flow.

## Steps

1. **Write `dockerfiles/Dockerfile`.** One file directly under
   `dockerfiles/` (single-container packages); the friendly name then
   defaults to the package name. Bake user/password/database into the
   image via env vars — the connection JSON in step 3 must match these
   byte-for-byte. Declare exactly one `EXPOSE <internal-port>`; the
   platform proxies it dynamically and forbids multiple ports per
   image (`DG-FACT-416`). For multi-container packages, use sub-folders
   `dockerfiles/<friendly-suffix>/Dockerfile`; the friendly name
   becomes `<package>-<friendly-suffix>` (`DG-FACT-414`).
   ```dockerfile
   FROM postgres
   ENV POSTGRES_PASSWORD datagrok
   ENV POSTGRES_USER datagrok
   ENV POSTGRES_DB world
   COPY world.sql /docker-entrypoint-initdb.d/
   EXPOSE 5432
   HEALTHCHECK --interval=5s --timeout=3s --retries=10 --start-period=10s \
     CMD pg_isready -U "$POSTGRES_USER" -d "$POSTGRES_DB" || exit 1
   ```
   Expected: `docker build dockerfiles/` succeeds locally. Reference:
   `packages/DBTests/dockerfiles/Dockerfile`.

2. **Add `dockerfiles/container.json` for resource allocation.** The
   article omits this file, but every production db-in-docker package
   ships one — without it the platform applies defaults (`memory: 512`,
   `cpu: 0.25`, `on_demand: false`) that are usually too small for a
   real DB (`DG-FACT-DRIFT-DBDOCKER-001`). Use the documented keys
   only — `shutdown_timeout` is in minutes; the variant
   `timeout_minutes` seen in `SureChembl/dockerfiles/container.json`
   is unsupported and silently ignored, do not propagate it.
   ```json
   {"cpu": 1, "memory": 2048, "on_demand": true, "shutdown_timeout": 30}
   ```
   Expected: container starts on first query and idles out after
   `shutdown_timeout` minutes. Schema:
   `{{ lattice.harness.help_develop_root }}/how-to/packages/docker-containers.md:41-86`.
   Reference: `packages/DBTests/dockerfiles/container.json`.

3. **Write `connections/<name>.json` with the `<DockerContainer>`
   placeholder.** The platform recognizes the connection as
   container-bound by the `parameters.server` value. Emit the
   two-segment form `${<Package>:<Friendly><DockerContainer>}`
   (`DG-FACT-412`). Normalize BOTH the package and friendly segments
   to first-letter capitalization with the rest lowercased — package
   `DBTests` becomes `Dbtests`; `SureChembl` becomes `Surechembl`.
   Literal `${DBTests:DBTests<DockerContainer>}` will NOT resolve
   (`DG-FACT-413`, `DG-FACT-DRIFT-DBDOCKER-002`). `db` pins the
   schema/database inside the image and must equal the Dockerfile's
   `ENV POSTGRES_DB`. Credentials go under `credentials.parameters`
   and must match the Dockerfile's `ENV POSTGRES_USER` /
   `ENV POSTGRES_PASSWORD`. NEVER add `port` — the platform resolves
   it from the running container (`DG-FACT-415`).
   ```json
   {
     "#type": "DataConnection",
     "name": "PostgresDocker",
     "friendlyName": "PostgresDocker",
     "parameters": {
       "server": "${Dbtests:Dbtests<DockerContainer>}",
       "db": "world"
     },
     "credentials": {
       "parameters": {"login": "datagrok", "password": "datagrok"}
     },
     "dataSource": "Postgres"
   }
   ```
   Expected: `grok publish` accepts the JSON without "missing
   dataSource" or unknown-server errors. Source:
   `{{ lattice.harness.help_develop_root }}/how-to/db/db-in-docker.md:13-65`.
   Reference: `packages/DBTests/connections/postgres-test-docker.json`.

4. **Publish the package.** Push the Dockerfile + container.json +
   connection JSON together so the platform builds the image and
   registers the connection in one round-trip.
   ```bash
   grok publish
   ```
   Expected: the connection appears under **Browse > Databases** as
   `<Package>:<Name>` (`DG-FACT-417`).

5. **(Production only) Move credentials out of the JSON.** For a
   throwaway DB the JSON-embedded credentials in step 3 suffice. For
   production, leave the JSON's `credentials.parameters` empty and
   POST the real login/password once after deploy so they route
   through the credentials store (`DG-FACT-409`).
   ```bash
   curl -X POST "$GROK_HOST/api/credentials/for/$PACKAGE.$CONNECTION" \
     -H "Authorization: $API_KEY" -H "Content-Type: application/json" \
     -d '{"login":"datagrok","password":"datagrok"}'
   ```
   Expected: the connection still resolves and now reads credentials
   from the secret store, not the world-readable JSON.

## Common failure modes

- **Connection deploys but queries fail with "unknown server".** The
  `server` placeholder used the original mixed-case package or friendly
  name (e.g. `${DBTests:DBTests<DockerContainer>}`). Normalize to
  first-letter-capitalized form: `Dbtests`, `Surechembl`
  (`DG-FACT-413`, `DG-FACT-DRIFT-DBDOCKER-002`).
- **Container never accepts traffic / build queue stalls.** The
  Dockerfile has zero or multiple `EXPOSE` directives. Exactly one
  internal port is allowed; remove extras (`DG-FACT-416`).
- **Container is OOM-killed or evicted on first real query.** No
  `container.json` was shipped, so the platform applied 512 MB / 0.25
  CPU defaults. Add `dockerfiles/container.json` with realistic
  `memory` / `cpu` (`DG-FACT-DRIFT-DBDOCKER-001`).
- **Connection JSON includes `parameters.port`.** Datagrok ignores the
  static port and the proxy still binds dynamically, but the platform
  rejects the JSON as malformed. Remove the field entirely
  (`DG-FACT-415`).
- **Login/password sit inside `parameters` or in `connString`.** That
  block is world-readable; the credentials store is bypassed. Move
  them under `credentials.parameters` and POST per step 5
  (`DG-FACT-409`).
- **Friendly-name segment is wrong.** With one Dockerfile directly in
  `dockerfiles/`, the friendly name equals the (normalized) package
  name. With sub-folders, it becomes `<package>-<folder>`. Mismatch →
  "unknown server" at query time (`DG-FACT-414`).

## Verification

- `grok publish` returns success and lists the new connection.
- In Datagrok, **Browse > Databases > `<Package>:<Name>`** opens; its
  schema tree populates after the first query.
- Fetch the connection object via `grok.functions.eval` (preserves
  the original package case — normalization applies ONLY inside the
  `${...}` placeholder) and call `.test()`. Pass: the dispatcher
  treats it as a normal Postgres/MariaDB connection (`DG-FACT-417`).
  ```typescript
  const conn = await grok.functions.eval('DbTests:PostgresDocker');
  await conn.test();
  ```
  Reference pattern: `packages/DBTests/src/connections/queries-test.ts:43-50`.

## See also

- Source: `{{ lattice.harness.help_develop_root }}/how-to/db/db-in-docker.md`
  (mirror: `docs/_internal/articles-mirror/how-to/db/db-in-docker.md`).
- Related: `{{ lattice.harness.help_develop_root }}/how-to/packages/docker-containers.md`
  (container.json schema), `db-in-plugin.md` (simpler alternative),
  `access-data.md` (connection-and-query spine this plugs into).
- Knowledge (`docs/_internal/knowledge/knowledge-graph.md`):
  `DG-FACT-412` – `DG-FACT-418`,
  `DG-FACT-DRIFT-DBDOCKER-001`, `DG-FACT-DRIFT-DBDOCKER-002`,
  `DG-FACT-409` (credentials POST).
- Related skills: [[create-package]], [[docker-containers]] (general
  container packaging), [[db-in-plugin]] (managed Postgres alternative),
  [[access-data]] (TS dispatcher), [[manage-credentials]] (step 5 POST).
