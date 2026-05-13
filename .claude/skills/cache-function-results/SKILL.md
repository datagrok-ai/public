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
   `//meta.cache: <mode>` accepts exactly `'client'`, `'server'`,
   `'all'`, or `'true'` (legacy synonym for `all`; prefer explicit
   forms) (`DG-FACT-340`).
   ```python
   #name: Example
   #language: python
   #meta.cache: client
   #input: string table [Data table]
   #output: int result
   ```

2. **Add `meta.cache.invalidateOn` with a 5-field cron if the answer
   can change.** Without it the cache never auto-invalidates. Use the
   dotted `meta.cache.invalidateOn` form — `meta.invalidateOn` without
   `meta.cache:` is a validator error (`DG-FACT-342`).
   ```text
   //meta.cache.invalidateOn: 0 * * * *        # hourly
   //meta.cache.invalidateOn: 0 0 * * *        # daily
   //meta.cache.invalidateOn: 0 0 1 * *        # monthly
   ```

3. **For a TypeScript function, prefer the JSDoc header form.** Place
   `//meta.cache:` / `//meta.cache.invalidateOn:` line comments
   directly above the `export`ed function — `grok api` codegens them
   to `package.g.ts` (`DG-FACT-465`):
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
   Equivalent decorator forms (use when the function already carries
   `@grok.decorators.func` metadata): nested `meta: {cache,
   cacheInvalidateOn}` or top-level `{cache, cacheInvalidateOn}`
   (`DG-FACT-343`).

4. **Stay inside the client-cache budget if you pick `client` or
   `all`.** IndexedDB-backed; ≤ 100 MB per function, ≤ 100 000 total
   records, output must be scalar / `dataframe` / `graphics` /
   `datetime`. Anything else must switch to `server` (`DG-FACT-341`).

5. **For SQL queries, cache at the connection level** via
   `cacheResults`, `cacheSchema`, `cacheInvalidateSchedule` keys in
   `connections/<name>.json` `parameters` (`DG-FACT-344`):
   ```json
   "parameters": { "server": "…", "port": 54327, "db": "test",
                   "cacheResults": true,
                   "cacheInvalidateSchedule": "0 1 * * *" }
   ```
   UI equivalent: right-click connection > **Edit…** > **Cache
   Results** / **Invalidate On**.

6. **(Optional) Clear the client cache at runtime** when a downstream
   source changed and the cron is too coarse (`DG-FACT-345`):
   ```typescript
   await grok.functions.clientCache.clear();         // drop all
   await grok.functions.clientCache.clear('myFn');   // one function id
   ```

## Common failure modes

- **`unsupported cache variable: <mode>`** — value outside the allowed
  set (`DG-FACT-340`).
- **`Can't use invalidateOn without cache`** — header has the
  invalidator but no `meta.cache:` line (`DG-FACT-342`).
- **`unsupported invalidateOn: <expr>`** — use a 5-field unix cron, not
  the 7-field Quartz form.
- **Cached output never refreshes.** `meta.cache.invalidateOn` was
  omitted. Add a cron or call `grok.functions.clientCache.clear(<id>)`.
- **`meta.cache: client` silently degrades.** Output type is outside
  IndexedDB eligibility — switch to `server` / `all` (`DG-FACT-341`).
- **Connection-level cache on but queries still recompute.**
  Platform-wide **Settings > Cache** switches still gate storage.

## See also

- Source: `help/develop/how-to/functions/cache-function-results.md`.
- Knowledge: `DG-FACT-340`–`345`, `DG-FACT-465`.
- Related skills: `access-data`, `python-functions`.
