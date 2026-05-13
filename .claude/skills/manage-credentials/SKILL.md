---
name: manage-credentials
version: 0.1.0
description: |
  Store and read secrets — API keys, database logins, service tokens —
  that a Datagrok package needs at runtime, without committing them to
  source. For plugin authors whose code calls a paid third-party API or
  a credentialed database. Produces a key-value record in the platform's
  secured credentials store (UI dialog or REST POST), plus a code-side
  reader using `_package.getCredentials()`.
  Use when asked to "store an API key for my plugin",
  "keep the DB password out of source", or
  "rotate a service-user token without redeploying".
triggers:
  - api key for plugin runtime
  - keep db password out of source
  - rotate service-user token
  - read my package secret
  - third-party token for plugin
  - secrets store for plugin
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# manage-credentials

## When to use

Your package calls an external resource (REST API, JDBC database, AWS
service) that needs a key or login at runtime, and the secret must NOT
live in the repo. Typical phrasings: "where do I keep the AWS access
key for my package", "the DB password belongs in the secured store, not
the connection JSON", "let a cron script rotate our service token".

## Prerequisites

- A published package on the target server — the store keys on package
  name, so the package must exist before POSTs to `for/<Pkg>` succeed.
- Personal API key from the user profile page (e.g.
  `https://public.datagrok.ai/u`). Sent raw as `Authorization` — NO
  `Bearer ` prefix (`DG-FACT-349`).
- For automation: a service user created via
  **Manage → Users → Actions → Add Service User** (`DG-FACT-350`).
- For connection credentials: a `DataConnection` JSON at
  `packages/<Pkg>/connections/<name>.json` (see `access-data`).

## Steps

Two routes share the same secured store but differ in REST path and read API (see DG-FACT-346): free-form package secrets at `for/<PackageName>`, connection credentials at `for/<PackageName>.<ConnectionName>`.

1. **Set package credentials in the UI** (one-off case).
   ```
   Manage → Packages → right-click <Pkg> → Credentials...
   ```
   Set **Credentials owner** to a user or group (see DG-FACT-350).

2. **Set credentials programmatically** (rotation scripts, CI). Payload keys are entity-specific (see DG-FACT-349).
   ```bash
   curl -X POST \
     "$GROK_HOST/api/credentials/for/$PACKAGE_NAME" \
     -H "Authorization: $API_KEY" \
     -H "Content-Type: application/json" \
     -d '{"login":"datagrok","password":"s3cr3t"}'
   ```
   For a connection, target `for/$PACKAGE_NAME.$CONNECTION_NAME`. Same body POSTed twice overwrites in place.

3. **Read the package's own credentials from code.** Use `.parameters[key]` (owner-only real values — see DG-FACT-347, DG-FACT-348).
   ```typescript
   import * as grok from 'datagrok-api/grok';
   import * as DG from 'datagrok-api/dg';

   export const _package = new DG.Package();

   async function getApiKey(): Promise<string> {
     const credentials = await _package.getCredentials();
     if (credentials == null) {
       grok.shell.info('Package credentials are not set.');
       return '';
     }
     return credentials.parameters['apiKey'];
   }
   ```
4. **Use `openParameters` when surfacing what's stored to a user dialog.** Redacted display form, safe for non-owners (see DG-FACT-348).
   ```javascript
   let p = await grok.dapi.packages.filter('<Pkg>').first();
   let c = await p.getCredentials();
   grok.shell.info(c ? c.openParameters : 'Credentials are not set.');
   ```

5. **Declare connection credentials in JSON at deploy time.** Add the `credentials.parameters` block to `packages/<Pkg>/connections/<name>.json` — at deploy, Datagrok strips it and stores under `<Pkg>.<Conn>`. Literal values only — there is NO `${VAR_NAME}` substitution (see DG-FACT-351):
   ```json
   {
     "#type": "DataConnection",
     "name": "PostgreSQLTest",
     "parameters": {
       "server": "db.example.com",
       "port": 5432,
       "db": "northwind"
     },
     "credentials": {
       "parameters": {
         "login": "datagrok",
         "password": "literal-password-here"
       }
     },
     "dataSource": "PostgresDart"
   }
   ```
   For a non-private repo, omit the `credentials` block entirely and POST after deploy with step 2.

## Common failure modes

- **HTTP 401 from REST POST.** `Authorization` has a `Bearer ` prefix or quotes — send raw key (see DG-FACT-349).
- **`_package.getCredentials()` returns `null` after a successful set.** Package not republished on the same host, OR caller isn't in the **Credentials owner** group (see DG-FACT-350).
- **Connection credentials missing after deploy.** Used `${DB_PASSWORD}` expecting env substitution — not supported (see DG-FACT-351). Inline literal or POST after deploy.
- **Rotation script breaks when its author leaves.** Used a personal API key. Create a service user via **Manage → Users → Actions → Add Service User** (see DG-FACT-350).
- **Secrets leak through `connString`.** Always put login/password in `credentials.parameters`.

## See also

- Source articles: `help/develop/how-to/packages/manage-credentials.md`.
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-346` (two entity types, read APIs), `DG-FACT-347` (return
  shape + null check), `DG-FACT-348` (`parameters` vs `openParameters`),
  `DG-FACT-349` (REST POST surface, Authorization rules), `DG-FACT-350`
  (owner-scoped, service users), `DG-FACT-351` (connection JSON block,
  no env substitution).
- Reference packages: `packages/NLP/src/package.ts:104-115` +
  `aws/nlp-user.py:40-53` (full set-and-read cycle);
  `packages/Chemspace/src/package.ts:405-413` (single-key read);
  `packages/KnimeLink/src/credentials.ts:6-22` (multi-key read with
  explicit error on missing fields);
  `packages/DBTests/connections/postgresql-test.json` (connection JSON).
- Related skills: `access-data`, `publish-packages`.
