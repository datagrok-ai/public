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

1. **Pick a stable, package-scoped storage name** to avoid cross-plugin clobbering.
   ```typescript
   import * as grok from 'datagrok-api/grok';
   const STORAGE_NAME = 'MyPackage.lastSearch';
   ```

2. **Write a single value with `add` (string only).** Sync; appends without replacing. `JSON.stringify` structured data first (`DG-FACT-398`).
   ```typescript
   const payload = {protocol: 'A12', similarity: 0.7};
   grok.userSettings.add(STORAGE_NAME, 'vault-42', JSON.stringify(payload));
   ```
   Flushes to server within ~10s (`DG-FACT-396`).

3. **Read a single value back with `getValue`.** Returns `string | undefined`; always handle `undefined`.
   ```typescript
   const raw = grok.userSettings.getValue(STORAGE_NAME, 'vault-42');
   const stored = raw ? JSON.parse(raw) as {protocol: string; similarity: number} : null;
   ```

4. **Overwrite the whole record atomically with `put`.** Takes a plain `{[key: string]: string}` (not a `Map` — `DG-FACT-DRIFT-USS-002`); existing keys are discarded.
   ```typescript
   grok.userSettings.put(STORAGE_NAME, {
     date: '2026-05-12',
     count: '42',
   });
   ```

5. **Merge keys into the record with `addAll`** — same shape as `put`, but existing keys are preserved (`DG-FACT-397`).
   ```typescript
   grok.userSettings.addAll(STORAGE_NAME, {
     lastVault: 'vault-42',
     lastQuery: 'CCO',
   });
   ```

6. **Iterate the whole record** via `get(name)` + `Object.keys`. Empty namespace returns `{}`, never `null` (`DG-FACT-DRIFT-USS-001`).
   ```typescript
   const entries = grok.userSettings.get(STORAGE_NAME) ?? {};
   for (const [key, raw] of Object.entries(entries)) {
     const value = JSON.parse(raw);
     // ... do something with key, value
   }
   ```

7. **Choose private vs shared explicitly.** `isPrivate` defaults to `true` and is part of the lookup key — pass `false` on EVERY read AND write for a shared store (`DG-FACT-394`, `DG-FACT-397`).
   ```typescript
   grok.userSettings.add(STORAGE_NAME, 'theme', 'dark', false);
   const theme = grok.userSettings.getValue(STORAGE_NAME, 'theme', false);
   ```

8. **Delete a single key** (leaves siblings intact).
   ```typescript
   grok.userSettings.delete(STORAGE_NAME, 'vault-42');
   ```

9. **Wipe the whole namespace** — runtime accepts `null`, but TS signature declares `key: string`; cast (`DG-FACT-DRIFT-USS-003`).
   ```typescript
   grok.userSettings.delete(STORAGE_NAME, null as unknown as string);
   ```

10. **(Optional) Chunk values that would exceed 5000 chars** — single value over the cap is silently rejected (`DG-FACT-395`). Split with a metadata key.
    ```typescript
    const MAX = 5000;
    const totalParts = Math.ceil(blob.length / MAX);
    for (let i = 0, n = 0; i < blob.length; i += MAX, n++)
      grok.userSettings.add(STORAGE_NAME, `${id}|${n}`, blob.slice(i, i + MAX));
    if (totalParts > 1)
      grok.userSettings.add(STORAGE_NAME, `${id}|meta`, JSON.stringify({parts: totalParts}));
    ```

## Common failure modes

- **`JSON.parse: unexpected token`** — wrote a non-string without `JSON.stringify` (`DG-FACT-398`).
- **Reload shows stale data** — writes flush every ~10s (`DG-FACT-396`); mirror in component state.
- **TS error on `delete(name, null)`** — cast `null as unknown as string` (`DG-FACT-DRIFT-USS-003`).
- **Shared value empty for other users** — `isPrivate` defaults to `true`; pass `false` on read too (`DG-FACT-397`).
- **Silent truncation** — value exceeded 5000 chars; chunk (`DG-FACT-395`).
- **Used `grok.dapi.userDataStorage`** — deprecated and async; migrate to `grok.userSettings` (`DG-FACT-394`).

## See also

- Source articles: `help/develop/how-to/data/user-settings-storage.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-394`..`DG-FACT-398` and drifts `DG-FACT-DRIFT-USS-001`..`003`.
- Related skills: `build-an-app` (the most common consumer of saved
  settings); `custom-package-settings-editors` (when the saved value
  should surface as a tunable knob in the Settings panel).
