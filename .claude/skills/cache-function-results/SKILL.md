---
name: cache-function-results
version: 0.1.0
description: |
  Make a Datagrok function (script, TypeScript wrapper, or SQL query)
  return the previously-computed answer for identical inputs instead of
  re-executing. For plugin authors whose functions are deterministic
  but slow — remote model calls, RDKit similarity passes, parameterized
  SQL queries hitting a stable backend. Produces a `meta.cache` /
  `meta.cache.invalidateOn` annotation on the function header (or the
  equivalent decorator field), plus the regenerated `package.g.ts`
  the platform reads at startup. Connection-wide caching is a separate
  switch on the connection JSON.
  Use when asked to "remember a previous answer", "skip re-running an
  expensive computation", or "stop re-paying for identical model calls".
triggers:
  - remember a previous answer
  - skip re-running an expensive computation
  - memoize results for identical inputs
  - avoid repeating an expensive api call
  - reuse last answer across users
  - persist query results between sessions
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# cache-function-results

## When to use

A Datagrok function is called repeatedly with the same arguments (LLM
completion, RDKit similarity, parameterized SQL) and you want the
second call to return the first call's result without re-executing.
The body must be effectively pure over the cache's lifetime — if the
backing data mutates, pick a coarse `invalidateOn` or do not cache.

## Steps

1. **Pick a cache mode and add the header annotation.**
   `//meta.cache: <mode>` accepts EXACTLY four validator-checked
   values — `'client'`, `'server'`, `'all'`, `'true'` (knowledge
   `DG-FACT-340`; `tools/bin/utils/utils.ts:191`,
   `tools/bin/commands/check.ts:448-449`). The article marks `true`
   as a legacy synonym for `all`
   (`cache-function-results.md:40-44`); prefer the explicit forms.
   ```python
   #name: Example
   #language: python
   #meta.cache: client
   #input: string table [Data table]
   #output: int result
   ```
   Expected: `grok check` reports no `unsupported cache variable`
   error. For TS packages, the emitted `src/package.g.ts` carries the
   same line — `packages/Chem/src/package.g.ts:719-720` (`client`),
   `packages/SureChembl/src/package.g.ts:7,17` (`all`).

2. **Add `meta.cache.invalidateOn` with a 5-field cron if the answer
   can change.** Without it the cache never auto-invalidates; entries
   only age out under the size/record-count limits. Validated by
   `utils.isValidCron` (`tools/bin/commands/check.ts:450-451`).
   Production cadences (knowledge `DG-FACT-342`):
   ```text
   //meta.cache.invalidateOn: 0 * * * *        # hourly  — Chem synthonSearchFunc
   //meta.cache.invalidateOn: 0 0 * * *        # daily   — SureChembl
   //meta.cache.invalidateOn: 0 0 1 * *        # monthly — Reinvent4, MolTrack
   ```
   Never write `meta.invalidateOn` without `meta.cache:` — the
   validator emits `Can't use invalidateOn without cache`
   (`check.ts:446-447`). The dotted `meta.cache.invalidateOn` is
   the only correct form.

3. **For a TypeScript function, the header form is the lowest-friction
   path; the decorator is an equivalent alternative.** Either lands the
   same annotation in `package.g.ts` (knowledge `DG-FACT-343`,
   `DG-FACT-465`).

   **Header form (preferred when adding cache to an existing
   function).** JSDoc-style line comments directly above the
   `export`ed function — `grok api` reads them and emits the same
   `//meta.cache:` into `src/package.g.ts`. Used in production by
   `packages/RevvitySignalsLink/src/package.ts:198-211` (`getUsers`,
   `getLibraries`, `getTerms` — all daily-invalidated remote-API
   wrappers, the same shape as a remote-model call):
   ```typescript
   //name: askRemoteModel
   //input: string prompt
   //output: string answer
   //meta.cache: all
   //meta.cache.invalidateOn: 0 0 * * *
   export async function askRemoteModel(prompt: string): Promise<string> {
     /* … */
   }
   ```

   **Decorator form (when the function already carries other
   `@grok.decorators.func` metadata).** Nest under `meta:`
   (article-documented; codegens to `//meta.cache:`):
   ```typescript
   @grok.decorators.func({
     name: 'Search Synthons',
     meta: { cache: 'client', cacheInvalidateOn: '0 * * * *' },
   })
   static async synthonSearchFunc(/* … */) { /* … */ }
   ```
   Reference `packages/Chem/src/package.ts:1669-1675`. The top-level
   form (`{cache: 'all', cacheInvalidateOn: '…'}` without `meta:`) is
   also accepted and codegens to `//cache:` —
   `packages/MolTrack/src/package.ts:155-161`.

4. **Stay inside the client-cache budget if you pick `client` or
   `all`** (knowledge `DG-FACT-341`, `cache-function-results.md:7-13`):
   IndexedDB-backed; per-function total ≤ 100 MB across all input
   variants; total records ≤ 100 000; output type must be scalar
   (`int`, `double`, `string`), `dataframe`, `graphics`, or
   `datetime`. Anything outside that set (`list`, `column`, …) cannot
   use client cache — switch to `server`. Server cache has no
   documented size or type restrictions
   (`cache-function-results.md:23-26`); `all` stacks both.

5. **For SQL queries, cache at the connection level instead of
   per-query.** `connections/<name>.json` exposes three
   `DataConnectionCacheProperties` keys (knowledge `DG-FACT-344`,
   `js-api/src/entities/types.ts:33-40`): `cacheResults: true`
   activates BOTH client and server caches for every query under
   the connection (not toggleable per-side here), `cacheSchema: true`
   also caches the introspected DB schema, and
   `cacheInvalidateSchedule` is the cron. Canonical
   `packages/DBTests/connections/cached-postgresql-test.json:5-11`:
   ```json
   "parameters": { "server": "…", "port": 54327, "db": "test",
                   "cacheResults": true,
                   "cacheInvalidateSchedule": "0 1 * * *" }
   ```
   UI equivalent: right-click the connection > **Edit…** > **Cache
   Results** / **Invalidate On** writes the same JSON.

6. **(Optional) Clear the client cache at runtime via
   `grok.functions.clientCache`** when a downstream source changed
   and the cron is too coarse (knowledge `DG-FACT-345`,
   `js-api/src/functions.ts:134,193-212`):
   ```typescript
   await grok.functions.clientCache.clear();         // drop all
   await grok.functions.clientCache.clear('myFn');   // one function id
   ```
   Also: `start()`, `stop()`, `cleanup()`, `getRecordCount()`,
   `isRunning`. The article documents only the Settings UI path.

7. **Build, check, publish.**
   ```bash
   npm install && grok check && grok publish <host>   # add --release once stable
   ```
   The user-side switches in **Settings > Cache** must also be on for
   entries to be stored — they default on, but the annotations alone
   do not force them.

## Common failure modes

- **`unsupported cache variable: <mode>`** — value outside the
  allowed set. Use `'all'`, `'server'`, `'client'`, or legacy
  `'true'` (knowledge `DG-FACT-340`).
- **`Can't use invalidateOn without cache`** — header has the
  invalidator but no `meta.cache:` line; the invalidator is parented
  on cache (knowledge `DG-FACT-342`).
- **`unsupported invalidateOn: <expr>`** — failed `isValidCron`. Use
  a 5-field unix cron (`0 0 * * *`), not the 7-field Quartz form
  (`0 0 * ? * * *`).
- **Cached output never refreshes.** `meta.cache.invalidateOn` was
  omitted. Add a cron matching the data's refresh cadence, or call
  `grok.functions.clientCache.clear(<id>)` on the known event.
- **`meta.cache: client` silently degrades.** Output type is outside
  IndexedDB eligibility (`list`, `column`, etc. — knowledge
  `DG-FACT-341`). Narrow the output, or switch to `server` / `all`.
- **Connection-level cache on but queries still recompute.**
  Platform-wide **Settings > Cache** switches still gate storage
  (`cache-function-results.md:17-19,30-32`).
- **TS function annotated via decorator but writer wanted minimal
  touch.** Use the JSDoc header form instead (`//meta.cache: …` lines
  directly above the `export function`) — `grok api` codegens it
  identically to the decorator form (knowledge `DG-FACT-465`,
  `packages/RevvitySignalsLink/src/package.ts:198-211`).

## See also

- Source articles: `help/develop/how-to/functions/cache-function-results.md`.
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-340` (mode enum + validator), `DG-FACT-341` (client cache
  size/record/type limits), `DG-FACT-342` (`invalidateOn` cron +
  parent-cache constraint), `DG-FACT-343` (header vs decorator
  forms), `DG-FACT-344` (connection-level keys), `DG-FACT-345`
  (`grok.functions.clientCache` runtime surface), `DG-FACT-465`
  (TS header form is preferred for minimal-touch caching).
- Reference packages: `packages/RevvitySignalsLink/src/package.ts:198-211`
  (TS header-form, `all` + daily — the canonical "add cache to an
  existing async TS function" shape; same pattern repeats at lines
  213-231, 243-260, 287-303, 322-340);
  `packages/Chem/src/package.ts:1669-1675` +
  `package.g.ts:719-720` (decorator-`meta`, `client` + hourly);
  `packages/MolTrack/src/package.ts:155-161` + `package.g.ts:32-35`
  (decorator-top-level, `all` + monthly);
  `packages/SureChembl/src/package.g.ts:7,17` (script-header,
  `all` + daily);
  `packages/DBTests/connections/cached-postgresql-test.json`
  (connection-level);
  `packages/ApiSamples/scripts/functions/caching-results.js`
  (dynamic `grok.functions.register({…, options: {cache,
  'cache.invalidateOn'}})`).
- Related skills: `access-data` (SQL queries / connections that
  benefit from connection-level caching); `python-functions`
  (caching annotations apply identically in script-header form).
