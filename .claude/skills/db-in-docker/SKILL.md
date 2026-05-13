---
name: db-in-docker
version: 0.3.1
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
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# db-in-docker

## When to use

Pick this path when your package needs a database engine the shared
platform RDS won't run:

- a non-Postgres DBMS (MariaDB, ClickHouse, Oracle);
- a Postgres image with extensions outside the shared instance's
  allowlist;
- required DBA-level isolation from other tenants (`DG-FACT-418`).

For ordinary CRUD storage that fits in the shared Postgres, prefer
[[db-in-plugin]] — simpler, platform-managed, far fewer resources. If
none of the three bullets above apply, stop and use that path instead.

## Prerequisites

- API key from `<GROK_HOST>/u` for the credentials POST in step 5.
- Familiarity with the connection-and-query spine from [[access-data]] —
  this skill only replaces step 1 of that flow.

## Steps

1. **Write `dockerfiles/Dockerfile`.** One file directly under `dockerfiles/` (single-container packages) — the friendly name then defaults to the package name. For multi-container layouts use sub-folders (see DG-FACT-414). Declare exactly one `EXPOSE <internal-port>` (see DG-FACT-416). Bake user/password/database into env vars — the connection JSON in step 3 must match them byte-for-byte.
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

2. **Add `dockerfiles/container.json` for resource allocation.** Without it the platform applies defaults too small for a real DB (see DG-FACT-DRIFT-DBDOCKER-001). Use documented keys only — `shutdown_timeout` is in minutes; do NOT use `timeout_minutes` (silently ignored).
   ```json
   {"cpu": 1, "memory": 2048, "on_demand": true, "shutdown_timeout": 30}
   ```
   Schema: `help/develop/how-to/packages/docker-containers.md:41-86`.

3. **Write `connections/<name>.json` with the `<DockerContainer>` placeholder.** Emit the two-segment form `${<Package>:<Friendly><DockerContainer>}` (see DG-FACT-412) and normalize BOTH segments to first-letter capitalization (`DBTests` → `Dbtests`, `SureChembl` → `Surechembl` — see DG-FACT-413, DG-FACT-DRIFT-DBDOCKER-002). `db` must match `ENV POSTGRES_DB`; credentials must match `ENV POSTGRES_USER` / `ENV POSTGRES_PASSWORD`. NEVER add `port` (see DG-FACT-415).
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

4. **Publish the package.** Push the Dockerfile + container.json + connection JSON together.
   ```bash
   grok publish
   ```
   Expected: the connection appears under **Browse > Databases** as `<Package>:<Name>` (see DG-FACT-417).

5. **(Production only) Move credentials out of the JSON.** Leave `credentials.parameters` empty and POST the real login/password once after deploy (see DG-FACT-409).
   ```bash
   curl -X POST "$GROK_HOST/api/credentials/for/$PACKAGE.$CONNECTION" \
     -H "Authorization: $API_KEY" -H "Content-Type: application/json" \
     -d '{"login":"datagrok","password":"datagrok"}'
   ```

## Common failure modes

- **"Unknown server" at query time.** `server` placeholder kept mixed-case (must be normalized — see DG-FACT-413, DG-FACT-DRIFT-DBDOCKER-002), or wrong friendly-name segment for multi-container layout (see DG-FACT-414).
- **Container never accepts traffic / build queue stalls.** Zero or multiple `EXPOSE` directives — exactly one internal port allowed (see DG-FACT-416).
- **Container OOM-killed on first real query.** No `container.json` shipped — see DG-FACT-DRIFT-DBDOCKER-001.
- **Connection JSON includes `parameters.port`.** Remove it (see DG-FACT-415).
- **Login/password sit inside `parameters` or `connString`.** Move them under `credentials.parameters` and POST per step 5 (see DG-FACT-409).

## See also

- Source: `help/develop/how-to/db/db-in-docker.md`
  (mirror: `docs/_internal/articles-mirror/how-to/db/db-in-docker.md`);
  related article `how-to/packages/docker-containers.md` (container.json schema).
- Knowledge (`docs/_internal/knowledge/knowledge-graph.md`):
  `DG-FACT-412`–`DG-FACT-418`, `DG-FACT-DRIFT-DBDOCKER-001`,
  `DG-FACT-DRIFT-DBDOCKER-002`, `DG-FACT-409`.
- Related skills: [[create-package]], [[docker-containers]],
  [[db-in-plugin]], [[access-data]], [[manage-credentials]].
