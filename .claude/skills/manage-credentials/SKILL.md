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
harness-authored: true
---

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

Two routes share the same secured store but differ in REST path and
read API (`DG-FACT-346`). Free-form package secrets → **package
credentials**, path `for/<PackageName>`, read via
`_package.getCredentials()`. Login/password for a declared
`DataConnection` → **connection credentials**, path
`for/<PackageName>.<ConnectionName>`, consumed implicitly by the
connection.

1. **Set package credentials in the UI** (one-off case).
   ```
   Manage → Packages → right-click <Pkg> → Credentials...
   ```
   Add key-value pairs; set **Credentials owner** to a user or group
   (`DG-FACT-350`).
   Expected: only members of the owner group see real values from
   `_package.getCredentials()`; non-members get `null`.

2. **Set credentials programmatically** (rotation scripts, CI).
   Payload keys are entity-specific — `{login, password}` for a DB,
   `{accessKeyId, secretAccessKey}` for AWS, etc. (`DG-FACT-349`).
   ```bash
   curl -X POST \
     "$GROK_HOST/api/credentials/for/$PACKAGE_NAME" \
     -H "Authorization: $API_KEY" \
     -H "Content-Type: application/json" \
     -d '{"login":"datagrok","password":"s3cr3t"}'
   ```
   For a connection, target `for/$PACKAGE_NAME.$CONNECTION_NAME`.
   Canonical Python reference: `packages/NLP/aws/nlp-user.py:47-53`.
   Expected: HTTP `200 OK`. Same body POSTed twice overwrites in place.

3. **Read the package's own credentials from code.** Lenient `== null`
   check, access via `.parameters[key]` (not `.openParameters`)
   (`DG-FACT-347`, `DG-FACT-348`).
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
   Reference: `packages/NLP/src/package.ts:104-115`,
   `packages/Chemspace/src/package.ts:405-413`.
   Expected: returns the value set in step 1 or 2 when the caller is
   in the owner group; otherwise `null`.

4. **Use `openParameters` when surfacing what's stored to a user
   dialog.** Different visibility — redacted display form, safe for
   non-owners (`DG-FACT-348`).
   ```javascript
   let p = await grok.dapi.packages.filter('<Pkg>').first();
   let c = await p.getCredentials();
   grok.shell.info(c ? c.openParameters : 'Credentials are not set.');
   ```
   Reference: `packages/ApiSamples/scripts/misc/package-credentials.js:6-8`.
   Expected: displays the redacted key set without leaking values.

5. **Declare connection credentials in JSON at deploy time.**
   In `packages/<Pkg>/connections/<name>.json`, add the
   `credentials.parameters` block alongside `parameters`
   (`DG-FACT-351`):
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
   At deploy, Datagrok strips the `credentials` block off the connection
   record and stores it under the entity `<Pkg>.<Conn>`. Values must be
   literal — there is NO `${VAR_NAME}` env-substitution. Canonical:
   `packages/DBTests/connections/postgresql-test.json:10-15`. For a
   non-private repo, omit the `credentials` block entirely and POST
   after deploy with step 2.

## Common failure modes

- **HTTP 401 from the REST POST.** The `Authorization` header has a
  `Bearer ` prefix or stray quotes around the key. Send the raw API
  key string from `/u` with NO prefix (`DG-FACT-349`).
- **`_package.getCredentials()` returns `null` after a successful
  set.** Either the package was not republished on the same host, or
  the caller isn't in the **Credentials owner** group — `parameters`
  is owner-gated (`DG-FACT-350`). Fix: re-run `grok publish dev`
  against the same `$GROK_HOST`; for non-owners, change owner to a
  shared group.
- **Connection credentials missing after deploy.** Connection JSON
  used `${DB_PASSWORD}` expecting env substitution. The platform does
  NOT expand `${...}` (`DG-FACT-351`). Fix: inline the literal value
  (private repo) or omit the block and POST to `for/<Pkg>.<Conn>`
  after deploy.
- **Rotation script breaks when its author leaves.** Script used a
  personal API key from `/u`. Fix: create a service user via
  **Manage → Users → Actions → Add Service User** and embed *its* key
  (`DG-FACT-350`; ref `packages/NLP/aws/nlp-user.py:50-52`).
- **Secrets leak through the connection string.** When `connString` is
  supplied, it persists in the connection record itself — not the
  secured store. Always put login/password in the
  `credentials.parameters` block.

## Verification

- POST the same body in step 2 a second time — endpoint is idempotent,
  HTTP `200 OK` both times.
- In the Datagrok console:
  ```javascript
  (await (await grok.dapi.packages.filter('<Pkg>').first()).getCredentials()).openParameters
  ```
  Returns the redacted key set you stored (`DG-FACT-348`).
- For a connection: open
  **Browse → Connections → \<Pkg\>:\<Conn\>**, click **Test connection** —
  it succeeds without re-prompting for a password.

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
