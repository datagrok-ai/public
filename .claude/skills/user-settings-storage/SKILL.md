---
name: user-settings-storage
description: Persist key-value settings per user (or shared across users) via grok.userSettings
---

# user-settings-storage

## When to use

Your package needs to remember a small piece of state across sessions —
last-used inputs, a saved layout, a per-vault search payload, a daily usage
counter. Triggers: "save user preferences", "remember last selection",
"share a setting across packages", "store JSON in user storage".

## Prerequisites

- A package importing `datagrok-api/grok` (`grok.userSettings` is the entry
  point — knowledge `DG-FACT-014`).
- Each value serialises to ≤5000 characters; longer payloads must be chunked
  (knowledge `DG-FACT-015`).
- Decide up front whether the data is private to the current user or shared
  across all users (knowledge `DG-FACT-016`). The flag is the same on every
  read and write — if you write with `isPrivate: false` you must read with
  `isPrivate: false`.

## Steps

1. **Pick a stable storage name and serialise the value.**
   `grok.userSettings` only stores strings; structured data must be
   `JSON.stringify`-ed on write and `JSON.parse`-ed on read
   (knowledge `DG-FACT-019`).
   ```typescript
   import * as grok from 'datagrok-api/grok';

   const STORAGE_NAME = 'MyPackage.lastSearch';   // namespace by package
   const inputObj = {x: 1, y: 2, z: 3};
   grok.userSettings.add(STORAGE_NAME, 'coords', JSON.stringify(inputObj));
   ```
   Expected: no exception; the entry is in memory and will be flushed to the
   server within ~10 seconds (knowledge `DG-FACT-018`).

2. **Read a single value back.**
   `getValue` returns `string | undefined`. Always handle `undefined` — the
   key may not exist yet, or you may be on a fresh device before the first
   server sync. Knowledge `DG-FACT-017` (synchronous, no `await` needed).
   ```typescript
   const raw = grok.userSettings.getValue(STORAGE_NAME, 'coords');
   const coords = raw ? JSON.parse(raw) as {x: number; y: number; z: number} : null;
   ```
   Expected: `coords` is your object, or `null` on first run.

3. **Replace a whole record atomically with `put`.**
   `put` overwrites the record at `name`; `addAll` merges into it. Both take
   a plain `{[key: string]: string}` object — NOT a JS `Map`
   (knowledge `DG-FACT-020`, drift `DG-FACT-DRIFT-006`).
   ```typescript
   grok.userSettings.put(STORAGE_NAME, {
     date: '2025-01-15',
     count: '42',
   });
   ```
   Expected: the previous record at `STORAGE_NAME` is gone; only the two keys
   above remain.

4. **Iterate the whole record.**
   You cannot list keys without first fetching the record. `get` returns
   `{[key: string]: string} | undefined` (article calls it a "Map" — that is
   the drift in `DG-FACT-DRIFT-006`).
   ```typescript
   const entries = grok.userSettings.get(STORAGE_NAME) ?? {};
   for (const [key, raw] of Object.entries(entries)) {
     const value = JSON.parse(raw);
     // do something with key, value
   }
   ```
   Expected: empty record (`{}`) when nothing is stored, never `null`.

5. **Choose private vs shared explicitly.**
   `isPrivate` is the third (or fourth) argument and defaults to `true`
   (knowledge `DG-FACT-016`, drift `DG-FACT-DRIFT-005`). If you need a
   shared store, pass `false` on every read AND every write.
   ```typescript
   // Shared across all users:
   grok.userSettings.add(STORAGE_NAME, 'theme', 'dark', false);
   const theme = grok.userSettings.getValue(STORAGE_NAME, 'theme', false);
   ```
   Expected: a different user logging in sees the same `theme` value.

6. **Delete a single key, or wipe the whole storage.**
   The TypeScript signature is `delete(name, key, isPrivate?)`. Passing
   `null` for `key` wipes the storage at runtime, but TS rejects `null` —
   cast or use a no-key call site (drift `DG-FACT-DRIFT-007`).
   ```typescript
   grok.userSettings.delete(STORAGE_NAME, 'coords');           // one key
   grok.userSettings.delete(STORAGE_NAME, null as any);         // wipe all
   ```
   Expected: subsequent `get(STORAGE_NAME)` returns `{}` after the wipe.

7. **(Optional) Chunk values over 5000 characters.**
   When a single value would exceed the limit (knowledge `DG-FACT-015`),
   split into N parts plus a metadata key. Reference pattern:
   `packages/ClinicalCase/src/utils/layout-utils.ts:18-35`.
   ```typescript
   const MAX = 5000;
   for (let i = 0, n = 0; i < blob.length; i += MAX, n++)
     grok.userSettings.add(STORAGE_NAME, `${id}|${n}`, blob.slice(i, i + MAX));
   grok.userSettings.add(STORAGE_NAME, `${id}|meta`, JSON.stringify({parts: Math.ceil(blob.length / MAX)}));
   ```
   Expected: load path reads `meta`, then concatenates the N parts.

## Common failure modes

- **Reload immediately after `add()` shows old/missing data.** Writes are
  in-memory and synced to the server every ~10 seconds (knowledge
  `DG-FACT-018`, drift `DG-FACT-DRIFT-008`). Fix: avoid forcing a reload
  right after a write, or store the value in component state too.
- **`JSON.parse: unexpected token` on `getValue` result.** You stored an
  object without `JSON.stringify`. The platform stores `[object Object]`
  literally. Fix: always `JSON.stringify` non-string values on write
  (knowledge `DG-FACT-019`).
- **TypeScript error on `delete(name, null)`.** The signature is
  `key: string`; `null` is not assignable (drift `DG-FACT-DRIFT-007`). Fix:
  cast (`null as any`) or migrate the wipe call into a helper that does
  the cast in one place.
- **Value lost silently when longer than 5000 chars.** The platform
  truncates / rejects (knowledge `DG-FACT-015`). Fix: chunk per step 7.
- **Shared data appears empty for other users.** A read used the default
  `isPrivate: true`. Pass `false` on the read too — the flag must match
  the write (knowledge `DG-FACT-016`).
- **Used `grok.dapi.userDataStorage` instead.** That class is
  `@deprecated` (knowledge `DG-FACT-021`); it is async and does not share
  data with `userSettings`. Fix: migrate to `grok.userSettings`.

## Verification

- After step 1, run `JSON.parse(grok.userSettings.getValue(STORAGE_NAME,
  'coords')!)` in the Datagrok console — it returns your object.
- Reload the page after waiting >10 seconds, repeat the read — value
  persists across the reload.
- Open the same package as a different user (or in an incognito session
  signed in as another account): private values are absent; shared values
  (written with `isPrivate: false`) are present.

## See also

- Source articles:
  - `help/develop/how-to/data/user-settings-storage.md`
- Knowledge:
  - `docs/_internal/knowledge/knowledge-graph.md` — facts `DG-FACT-014`
    through `DG-FACT-021` and drifts `DG-FACT-DRIFT-005..008`.
- Related skills:
  - `build-an-app` (apps are the most common consumer of saved settings).
