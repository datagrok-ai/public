---
name: manage-credentials
description: Store and read secured credentials for a Datagrok package or one of its data connections
harness-authored: true
---

# manage-credentials

## When to use

A package needs an API key, login/password, or token at runtime and you
don't want it in the repo. Triggers: "where do I keep the AWS key,"
"hide the DB password," "rotate the service-user token," "read my own
package's secret from `package.ts`."

## Prerequisites

- A published package on the target Datagrok server (the credentials
  store keys on the package name; the package must exist before you
  POST to it).
- Your personal **API key** from the user-profile page (e.g.
  `https://public.datagrok.ai/u`). Sent raw as the `Authorization`
  header — NO `Bearer` prefix (knowledge `DG-FACT-349`).
- For programmatic rotation: a **service user** created via
  **Manage | Users | Actions | Add Service User** so the script doesn't
  depend on a human account (knowledge `DG-FACT-350`).
- For connection credentials: the connection already declared in
  `packages/<Pkg>/connections/<name>.json` (see `access-data` skill).

## Steps

1. **Pick the entity type.**
   - Free-form package secrets (API keys, tokens, anything not tied to
     a connection) → **package credentials**, REST path
     `for/<PackageName>`.
   - Login/password for a declared `DataConnection` → **connection
     credentials**, REST path `for/<PackageName>.<ConnectionName>`.

   Both routes funnel through the same credentials store; only the
   READ API differs (knowledge `DG-FACT-346`).

2. **(UI path) set package credentials by hand.**
   In Datagrok: **Manage → Packages**, right-click the package,
   choose **Credentials...**, add key-value pairs, set
   **Credentials owner** to a user or group.
   Expected: only members of the owner group will see real values back
   from `_package.getCredentials()`; others get `null` (knowledge
   `DG-FACT-350`).

3. **(REST path) set credentials programmatically.**
   ```bash
   curl -X POST \
     "$GROK_HOST/api/credentials/for/$PACKAGE_NAME" \
     -H "Authorization: $API_KEY" \
     -H "Content-Type: application/json" \
     -d '{"login":"...","password":"..."}'
   ```
   For a connection credential, target the dotted form:
   `$GROK_HOST/api/credentials/for/$PACKAGE_NAME.$CONNECTION_NAME`
   (knowledge `DG-FACT-349`).
   Expected: HTTP `200 OK`. Body keys are entity-specific (`login` /
   `password` for a DB; `accessKeyId` / `secretAccessKey` for AWS —
   canonical: `packages/NLP/aws/nlp-user.py:40-53`).

4. **Read your own package's credentials from code.**
   Use the lenient `== null` form and name the local `credentials`, not
   `credentialsResponse` (drift `DG-FACT-DRIFT-CREDS-002`):
   ```typescript
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
   Reference: `packages/NLP/src/package.ts:104-115`.
   Expected: returns the value set in step 2 or 3 when the caller is a
   member of the **Credentials owner** group; otherwise `null`.

5. **Pick `parameters` vs `openParameters` for the right reason.**
   - `credentials.parameters[key]` — the real secret. Use this when
     your package consumes the value (auth header, AWS SDK init).
     The map is hidden from non-owners.
   - `credentials.openParameters` — the redacted display form. Use this
     only when you intentionally surface "credentials exist" to a user
     dialog. Reference: `packages/ApiSamples/scripts/misc/package-credentials.js:6-8`.

   Rule (knowledge `DG-FACT-348`): own-code reading own secret →
   `parameters`; UI listing → `openParameters`.

6. **(Connection JSON) declare credentials at deploy time.**
   In `packages/<Pkg>/connections/<name>.json`, add the
   `credentials.parameters` block alongside the `parameters` block:
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
   At deploy time, Datagrok strips the `credentials` block from the
   connection record and moves it into the secured store under
   `<Pkg>.<Conn>` (knowledge `DG-FACT-351`). Reference:
   `packages/DBTests/connections/postgresql-test.json`.

   **Values must be literal.** There is no `${VAR_NAME}` substitution
   at deploy time (drift `DG-FACT-DRIFT-CREDS-001`). For a repo you
   cannot keep private, OMIT the `credentials` block entirely and POST
   credentials after deploy with step 3 (path b).

7. **Never put secrets in `connString`.**
   When the `connString` parameter is supplied, all other connector
   parameters are ignored, but `connString` itself ends up in the
   connection record — not the secured store. Keep credentials in the
   `credentials.parameters` block (knowledge `DG-FACT-507`).

## Common failure modes

- **`_package.getCredentials()` returns `null` and you ARE the owner.**
  Package wasn't republished after the credentials were set, OR the
  POST in step 3 targeted a stale instance. Fix: re-run
  `grok publish dev` against the same `$GROK_HOST` you POSTed to and
  reload the package.
- **HTTP 401 from the REST call.** The `Authorization` header has a
  `Bearer ` prefix or quoted value. Send the raw API key string from
  `/u` with NO prefix (knowledge `DG-FACT-349`).
- **Credentials silently missing after deploy.** The connection JSON
  used `${DB_PASSWORD}` (or similar) expecting env substitution. The
  platform does NOT expand `${...}` in connection JSONs (drift
  `DG-FACT-DRIFT-CREDS-001`). Fix: either inline literal values (if
  the repo is private) or omit the `credentials` block and POST to
  `/api/credentials/for/<Pkg>.<Conn>` after deploy.
- **`creds.parameters['x']` is `undefined` for a teammate but works
  for you.** The teammate isn't in the **Credentials owner** group;
  `parameters` is owner-gated (knowledge `DG-FACT-350`). Fix: change
  owner to a group both of you are in (e.g. `All users` for
  truly-public secrets, or a shared team group).
- **Rotation script breaks when its author leaves.** Script uses the
  author's personal API key from `/u`. Fix: create a service user via
  **Manage | Users | Actions | Add Service User** and embed *its* API
  key (knowledge `DG-FACT-350`; reference
  `packages/NLP/aws/nlp-user.py:50-52`).
- **Showing `parameters` to a non-owner user dialog leaks structure.**
  Use `openParameters` when surfacing "what's stored" in UI; reserve
  `parameters` for the package's own consumption (knowledge
  `DG-FACT-348`).

## Verification

- After step 3, re-run the POST with the same body; it should still
  return 200 (the endpoint is idempotent — same keys overwrite).
- In the Datagrok console:
  ```javascript
  (await (await grok.dapi.packages.filter('<Pkg>').first()).getCredentials()).openParameters
  ```
  Returns the redacted key set you stored.
- For connection credentials: open the connection in
  **Browse → Connections → <Pkg>:<Conn>**, click **Test connection**
  — it should succeed without the dialog re-prompting for a password.

## See also

- Source articles:
  - `help/develop/how-to/packages/manage-credentials.md`
- Knowledge:
  - `docs/_internal/knowledge/knowledge-graph.md` — facts
    `DG-FACT-346..351`, related `DG-FACT-499` (REST POST surface),
    `DG-FACT-507` (connString hazard); drifts
    `DG-FACT-DRIFT-CREDS-001` (no `${VAR}` substitution) and
    `DG-FACT-DRIFT-CREDS-002` (`== null`, name `credentials`).
- Related skills:
  - `access-data` (connection JSON shape and `dataSource` choice).
  - `db-in-plugin` (declaring a package-owned connection that needs
    credentials).
  - `publish-packages` (the `grok publish dev` step that makes a newly
    created package visible to `for/<Pkg>` POSTs).
