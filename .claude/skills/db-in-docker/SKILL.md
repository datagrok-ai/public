---
name: db-in-docker
description: Ship a database inside the package's Docker container and connect to it from Datagrok
---

# db-in-docker

## When to use

Your package needs to ship a database the platform-managed Postgres can't
host — a non-Postgres engine, a Postgres plugin RDS doesn't support, or a
preloaded dataset baked into the image. Triggers: "bundle SQLite/MySQL with
my package", "use pgvector / postgis", "ship a demo DB with seed data", or
"my plugin needs its own DB and I don't want a plugin Postgres". For an
empty platform-managed Postgres, prefer the simpler `db-in-plugin` flow.

## Prerequisites

- A package scaffold (`grok create <Name>`). Paths below are relative to
  the package root.
- Familiarity with the `connections/*.json` shape — see the `access-data`
  skill (knowledge `DG-FACT-033`, `DG-FACT-040`).
- Docker CLI installed locally for `docker build` smoke tests; `grok-spawner`
  must be running on the target Datagrok instance for image build/run.
- A user API key from `<GROK_HOST>/u` for transferring credentials post-deploy.

## Steps

1. **Decide the dockerfile layout.**
   One container per package — put the Dockerfile at
   `dockerfiles/Dockerfile`. For multiple containers, create one
   subfolder per image: `dockerfiles/<folder>/Dockerfile`. Container
   "friendly name" follows from the layout: single → `<PackageName>`;
   multi → `<PackageName>-<folder>` (knowledge `DG-FACT-043`).
   ```bash
   mkdir -p dockerfiles
   ```
   Expected: a `dockerfiles/` directory at the package root.

2. **Author the Dockerfile.**
   Set the DB engine, seed credentials via env vars, COPY init scripts,
   and EXPOSE exactly one port. EXPOSE is mandatory and must appear once
   only — multiple EXPOSE directives are unsupported (knowledge
   `DG-FACT-046`). A `HEALTHCHECK` is recommended so the platform sees
   when the DB is ready (knowledge `DG-FACT-048`).
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
   Expected: `dockerfiles/Dockerfile` (or `dockerfiles/<folder>/Dockerfile`)
   builds locally with `docker build dockerfiles -t pkg-db`.

3. **(Optional) Add `container.json` next to the Dockerfile.**
   Databases are good candidates for on-demand startup with an idle
   timeout (knowledge `DG-FACT-049`). The article on `db-in-docker.md`
   does NOT mention this file — see `docker-containers.md` for the full
   schema (drift `DG-FACT-DRIFT-014`).
   ```json
   {"on_demand": true, "shutdown_timeout": 30}
   ```
   Expected: `dockerfiles/container.json` exists; container starts on
   first DB request and shuts down after 30 idle minutes. Omit the file
   to keep the container running continuously.

4. **Write the connection JSON — only programmatic creation works.**
   Drop a file under `connections/<name>.json`. There is NO UI form for
   Docker-bound connections; the platform identifies one by the
   `<DockerContainer>` token in `parameters.server` (knowledge
   `DG-FACT-042`). The `port` parameter MUST be omitted — Datagrok
   resolves it dynamically from the container address (knowledge
   `DG-FACT-044`). Add `db` to select a schema/database name; it must
   match the engine's `POSTGRES_DB` (or equivalent) from step 2
   (knowledge `DG-FACT-045`).
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
   Expected: a JSON file at `connections/postgres-docker.json`. For
   `dataSource` values see knowledge `DG-FACT-034`.

5. **Get the placeholder casing right.**
   `parameters.server` follows
   `${<PackageName>:<ContainerFriendlyName><DockerContainer>}`. Both
   segments use first-letter-only-capitalized form — package `DBTests`
   becomes `Dbtests:Dbtests`, NOT `DBTests:DBTests` (knowledge
   `DG-FACT-043`, drift `DG-FACT-DRIFT-012`). The article example
   (`${Test:Test<DockerContainer>}`) hides this rule because `Test` is
   already first-letter-capitalized.
   Expected: `${Dbtests:Dbtests<DockerContainer>}` for a `DBTests`
   package; `${Foo:Foo-worker<DockerContainer>}` for `dockerfiles/worker/`
   inside a `Foo` package.

6. **Match Dockerfile credentials to the connection JSON.**
   Datagrok routes `credentials.parameters.{login,password}` to the
   container at runtime; the Dockerfile's `ENV POSTGRES_USER` /
   `ENV POSTGRES_PASSWORD` seed the DB at build time. The two MUST be
   identical or the connection will be rejected (knowledge `DG-FACT-047`).
   Expected: `login` / `password` in the JSON exactly equal
   `POSTGRES_USER` / `POSTGRES_PASSWORD` in the Dockerfile.

7. **Publish and transfer the credentials.**
   Same flow as any other connection — see `access-data` step 2.
   ```bash
   webpack
   grok publish dev
   curl -X POST "$GROK_HOST/api/credentials/for/$PACKAGE_NAME.PostgresDocker" \
     -H "Authorization: $API_KEY" -H "Content-Type: application/json" \
     -d '{"login":"datagrok","password":"datagrok"}'
   ```
   Expected: publish exits `0`; credentials POST returns `200 OK`. The
   image lands in **Manage → Docker** as building (grey blink) → ready
   (green).

8. **Verify and run a query.**
   The connection appears in **Browse → Databases**. Author a
   parameterised query (`access-data` step 3) and call it from JS:
   ```typescript
   import * as grok from 'datagrok-api/grok';
   const df = await grok.data.query(`${_package.name}:Cities`, {limit: 10});
   grok.shell.addTableView(df);
   ```
   Expected: a TableView opens with rows from inside the containerized DB.

## Common failure modes

- **Server placeholder rejected.** Wrong casing (`${DBTests:DBTests…}`)
  or wrong segment count fails resolution; canonical form is
  first-letter-only-capitalized for both segments (knowledge
  `DG-FACT-043`, drift `DG-FACT-DRIFT-012`). Fix: rename to
  `${Dbtests:Dbtests<DockerContainer>}`. Mirror the casing used by
  `packages/DBTests/connections/postgres-test-docker.json`.
- **Multi-container friendly-name confusion.** With several
  `dockerfiles/<folder>/Dockerfile`, friendly name becomes
  `<PackageName>-<folder>` — not the package name alone (knowledge
  `DG-FACT-043`). Fix: use `${Foo:Foo-worker<DockerContainer>}` for
  `dockerfiles/worker/`.
- **Connection times out / silently never resolves.** `parameters.port`
  is set in the JSON; Datagrok cannot map the dynamic container port
  (knowledge `DG-FACT-044`). Fix: remove the `port` line entirely.
- **Image build fails with "EXPOSE" error or container has no exposed
  port.** Multiple `EXPOSE` directives, or no `EXPOSE` at all — the
  platform allows exactly one (knowledge `DG-FACT-046`). Fix: collapse
  to a single `EXPOSE <port>`.
- **Auth fails on first query but Dockerfile builds fine.** Login or
  password in `credentials.parameters` does not match
  `POSTGRES_USER`/`POSTGRES_PASSWORD` (knowledge `DG-FACT-047`). Fix:
  copy the values from the Dockerfile into the JSON verbatim.
- **`db` property points at a database the container never created.**
  Postgres images only auto-create the database named by `POSTGRES_DB`
  (knowledge `DG-FACT-045`). Fix: align `parameters.db` with
  `ENV POSTGRES_DB`, or COPY init SQL that runs `CREATE DATABASE` into
  `/docker-entrypoint-initdb.d/`.

## Verification

- `docker build dockerfiles -t pkg-db` exits `0` (local smoke test).
- `grok publish dev` exits `0` and prints no `<DockerContainer>`-related
  warnings.
- **Manage → Docker** shows a green dot for the container card.
- **Browse → Databases** lists the new connection; clicking it opens
  the schema browser.
- A test query (`grok.data.query` or the in-platform query editor)
  returns rows.

## See also

- Source articles:
  - `help/develop/how-to/db/db-in-docker.md`
  - `help/develop/how-to/packages/docker-containers.md` (container.json,
    EXPOSE rule)
- Knowledge:
  - `docs/_internal/knowledge/knowledge-graph.md` — facts `DG-FACT-042`
    through `DG-FACT-049` and drifts `DG-FACT-DRIFT-012..015`.
- Related skills:
  - `access-data` (sibling — covers connection JSON shape, query
    authoring, credentials POST).
  - `db-in-plugin` (preferred when a managed empty Postgres is enough).
