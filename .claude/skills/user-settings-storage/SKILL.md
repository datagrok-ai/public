---
name: user-settings-storage
version: 0.2.1
description: |
  Persist small key-value blobs from a Datagrok package via
  `grok.userSettings` — last-used inputs, the row a user had selected,
  a theme choice, a per-vault search payload — so state survives page
  reload and (optionally) is shared across users. For plugin authors
  who need lightweight, platform-managed per-user state without
  standing up a database or hitting `localStorage`.
  Use when asked to "persist a filter selection across sessions",
  "remember the last compound a user looked at", or "stash a small
  per-user JSON blob the platform syncs".
triggers:
  - persist state between sessions
  - remember last user selection
  - stash per-user json blob
  - save preferences across reloads
  - share a value across packages
  - keep settings after refresh
allowed-tools:
  - Read
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# user-settings-storage

## When to use

Your Datagrok package needs a tiny piece of state — a saved filter, the
last-opened compound, a theme choice — that must outlive page reloads
and follow the user across devices, with no separate database.

## Prerequisites

- A package with `import * as grok from 'datagrok-api/grok';` —
  `grok.userSettings` is the singleton entry point (`DG-FACT-394`).
- Decide upfront whether the record is private to the user (default)
  or shared; `isPrivate` must match on every read and write of it.
- Each `value` is a `string` ≤ 5000 characters; pass structured data
  through `JSON.stringify` / `JSON.parse` (`DG-FACT-395`, `DG-FACT-398`).

## Steps

1. **Pick a stable, package-scoped storage name.** Namespace by package
   so two plugins don't clobber each other's keys. The singleton lives
   at `js-api/grok.ts:29`, typed in `js-api/src/user_settings_storage.ts:12-75`.
   ```typescript
   import * as grok from 'datagrok-api/grok';
   const STORAGE_NAME = 'MyPackage.lastSearch';
   ```
   Expected: a constant used by every read/write call site below.

2. **Write a single value with `add` (string only).** Synchronous,
   returns `void`, appends to (does not replace) the record at
   `STORAGE_NAME`. `JSON.stringify` structured data first. Signature at
   `js-api/src/user_settings_storage.ts:23`. Knowledge: `DG-FACT-398`.
   ```typescript
   const payload = {protocol: 'A12', similarity: 0.7};
   grok.userSettings.add(STORAGE_NAME, 'vault-42', JSON.stringify(payload));
   ```
   Expected: no exception. The platform flushes to the server within
   ~10 seconds (knowledge `DG-FACT-396`). Pattern:
   `packages/CddVaultLink/src/search-function-editor.ts:144-154`.

3. **Read a single value back with `getValue`.** Returns
   `string | undefined` — always handle `undefined` (key may not exist
   yet, or the sync window hasn't elapsed on a fresh device). Signature
   at `js-api/src/user_settings_storage.ts:62`.
   ```typescript
   const raw = grok.userSettings.getValue(STORAGE_NAME, 'vault-42');
   const stored = raw ? JSON.parse(raw) as {protocol: string; similarity: number} : null;
   ```
   Expected: `stored` is the parsed object, or `null` on first run.
   Pattern: `packages/CddVaultLink/src/search-function-editor.ts:113-119`.

4. **Overwrite the whole record atomically with `put`.** Takes a plain
   `{[key: string]: string}` object — not a JS `Map` (knowledge
   `DG-FACT-DRIFT-USS-002`). Every existing key at `STORAGE_NAME` is
   discarded. Signature at `js-api/src/user_settings_storage.ts:43`.
   ```typescript
   grok.userSettings.put(STORAGE_NAME, {
     date: '2026-05-12',
     count: '42',
   });
   ```
   Expected: previous record at `STORAGE_NAME` is gone; only the two
   keys above remain. Pattern: `packages/KnimeLink/src/function-cache.ts:29-34`.

5. **Merge keys into the record with `addAll`.** Same `{[key: string]:
   string}` shape as `put`, but existing keys are preserved; only the
   argument's keys are added or overwritten. Signature at
   `js-api/src/user_settings_storage.ts:33`. Knowledge: `DG-FACT-397`.
   ```typescript
   grok.userSettings.addAll(STORAGE_NAME, {
     lastVault: 'vault-42',
     lastQuery: 'CCO',
   });
   ```
   Expected: `get(STORAGE_NAME)` returns prior keys plus `lastVault`
   and `lastQuery`.

6. **Iterate the whole record.** There is no `keys()` / `has()` /
   `list()`; call `get(name)` then `Object.keys(...)`. `get` returns
   `{[key: string]: string} | undefined`. An empty namespace returns
   `{}` (the article's `entries !== null` check misleads — knowledge
   `DG-FACT-DRIFT-USS-001`).
   ```typescript
   const entries = grok.userSettings.get(STORAGE_NAME) ?? {};
   for (const [key, raw] of Object.entries(entries)) {
     const value = JSON.parse(raw);
     // ... do something with key, value
   }
   ```
   Expected: empty object `{}` when nothing is stored; never `null`.
   Pattern: `packages/PowerPack/src/utils.ts:24-39`.

7. **Choose private vs shared explicitly.** `isPrivate` is the last
   argument on every method, defaults to `true`, and is part of the
   lookup key — pass `false` on every read AND every write for a shared
   store. Knowledge: `DG-FACT-394`, `DG-FACT-397`.
   ```typescript
   grok.userSettings.add(STORAGE_NAME, 'theme', 'dark', false);
   const theme = grok.userSettings.getValue(STORAGE_NAME, 'theme', false);
   ```
   Expected: a different signed-in user sees the same `theme` value.

8. **Delete a single key.** Leaves the rest of the record intact.
   Signature at `js-api/src/user_settings_storage.ts:72`.
   ```typescript
   grok.userSettings.delete(STORAGE_NAME, 'vault-42');
   ```
   Expected: `getValue(STORAGE_NAME, 'vault-42')` returns `undefined`;
   sibling keys still resolve.

9. **Wipe the whole namespace.** The signature declares `key: string`,
   but the runtime accepts `null` — cast it, since the article's
   `delete(name, null)` form fails strict TS (knowledge
   `DG-FACT-DRIFT-USS-003`).
   ```typescript
   grok.userSettings.delete(STORAGE_NAME, null as unknown as string);
   ```
   Expected: subsequent `get(STORAGE_NAME)` returns `{}`.

10. **(Optional) Chunk values that would exceed 5000 chars.** A single
    `value` over the cap is silently rejected (knowledge `DG-FACT-395`).
    Split into N parts plus a metadata key; reconstruct on read.
    Reference pattern: `packages/ClinicalCase/src/utils/layout-utils.ts:14-72`.
    ```typescript
    const MAX = 5000;
    const totalParts = Math.ceil(blob.length / MAX);
    for (let i = 0, n = 0; i < blob.length; i += MAX, n++)
      grok.userSettings.add(STORAGE_NAME, `${id}|${n}`, blob.slice(i, i + MAX));
    if (totalParts > 1)
      grok.userSettings.add(STORAGE_NAME, `${id}|meta`, JSON.stringify({parts: totalParts}));
    ```
    Expected: the load path reads `|meta`, then concatenates the N parts
    back into the original blob.

## Common failure modes

- **`JSON.parse: unexpected token` on a `getValue` result.** Object
  stored without `JSON.stringify` — platform serialised it as
  `"[object Object]"`. Fix: `JSON.stringify` non-string values on write
  (`DG-FACT-398`).
- **Page reload right after `add()` shows old or missing data.** Writes
  flush every ~10 seconds (`DG-FACT-396`). Fix: don't force a reload
  immediately after a write, or mirror in component state.
- **TypeScript error on `delete(name, null)`.** Signature is
  `key: string` (`DG-FACT-DRIFT-USS-003`). Fix: cast to
  `null as unknown as string`, or wrap the wipe in a helper.
- **Shared value appears empty for other users.** The read used default
  `isPrivate: true`. Fix: pass `false` on the read — the flag is part
  of the storage key (`DG-FACT-397`).
- **Silent truncation of a long value.** Value exceeded 5000 chars
  (`DG-FACT-395`). Fix: chunk per step 10.
- **Used `grok.dapi.userDataStorage` instead.** Class is `@deprecated`
  (`js-api/src/dapi.ts:168,712,716`), async, and does not share data
  with `grok.userSettings`. Fix: migrate to `grok.userSettings`
  (`DG-FACT-394`).

## See also

- Source articles: `help/develop/how-to/data/user-settings-storage.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-394`..`DG-FACT-398` and drifts `DG-FACT-DRIFT-USS-001`..`003`.
- Related skills: `build-an-app` (the most common consumer of saved
  settings); `custom-package-settings-editors` (when the saved value
  should surface as a tunable knob in the Settings panel).
