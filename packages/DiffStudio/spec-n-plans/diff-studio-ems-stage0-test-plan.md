# Diff Studio — EMS Stage 0 Test Plan (schema + typed client)

Verifies the **foundation** of the EMS migration end-to-end on a real server: the `diffstudio`
domain schema deploys, the generated typed client (`src/generated/db.ts`) round-trips against the
backend, server-side validation holds, and row-level privacy reproduces MyFiles semantics.

This is the **exit criterion for Stage 0**. It does not touch UI or user-facing behavior (the
`diffStudioEms` flag stays `false`); it exercises only the data plane.

## Preconditions

1. **Deployed to dev.** `npx grok publish dev` completed (local `datagrok-tools` 6.5.7, so
   `databases/` ships and the server applies the schema). Give the server a few seconds after
   publish to apply the schema.
2. **Dev stand has EMS enabled** (domain databases). Confirmed available.
3. **Two user accounts** for the privacy test (TC8): user **A** (the publisher) and user **B**
   (any other account, or a second browser / incognito session).
4. **Where to run:** the dev stand's **JS console** (browser platform console). All snippets use
   the generic client `grok.dapi.domains.table('diffstudio.model')` — no package import needed.
   (The typed `diffstudioDb` from `src/generated/db.ts` is equivalent and used from package code.)

Set up a shared handle once per console session:

```js
const db = grok.dapi.domains.table('diffstudio.model');
const src = '#name: CRUD test\n#equations:\n  dx/dt = -x\n#argument: t\n  t0 = 0\n  t1 = 1\n  h = 0.1\n#inits:\n  x = 1';
```

---

## TC0 — Schema and table exist

**Goal:** the server accepted `schema.json` and created `diffstudio.model` (with system columns,
`securityMode: row`, `defaultRowVisibility: none`, and the `method` choices constraint).

```js
await db.query({});           // → [] (empty array), NO exception
```

**Pass:** returns an array (empty on a fresh schema) without throwing `DomainNotFoundError` /
manifest errors. Optionally confirm in the UI: **Data → Domains → diffstudio → model** shows the
declared columns plus `id/version/created_on/updated_on/author_id/is_deleted`.

**Fail signals:** "schema not found", "table not found", or a manifest-validation error → the
`row` + `defaultRowVisibility: none` combo or a column option was rejected server-side; capture
the exact message.

---

## TC1 — Insert

**Goal:** a row is created; system columns are populated; datetimes materialize as dayjs.

```js
const [r] = await db.insert({name: 'CRUD test', source: src, method: 'ros34prw'});
console.log(r.id, r.version, r.author_id);
console.log('created_on is dayjs:', r.created_on && typeof r.created_on.format === 'function');
```

**Pass:** `r.id` is a UUID, `r.version` is `1` (or the initial counter), `r.author_id` is the
current user, `created_on.format(...)` works (dayjs, not a string).

---

## TC2 — Get by id

```js
const got = await db.get(r.id);
console.log(got.name, got.source === src, got.method);
```

**Pass:** returns the row; `name === 'CRUD test'`, `source` round-trips byte-for-byte,
`method === 'ros34prw'`. `await db.get('00000000-0000-0000-0000-000000000000')` → `null`.

---

## TC3 — Query / filter

```js
await db.query({filter: 'name = "CRUD test"'});          // → [row]
await db.query({filter: 'method = "ros34prw"', sort: '!created_on', limit: 10});
```

**Pass:** the filter returns the inserted row; sort/limit are honored.

---

## TC4 — Update + optimistic concurrency

**Goal:** update in place bumps the version; a stale version is rejected (this is the mechanism
Stage 1's save-in-place relies on).

```js
const before = await db.get(r.id);
await db.update(r.id, {description: 'updated'}, {version: before.version});
const after = await db.get(r.id);
console.log('version bumped:', after.version === before.version + 1, 'desc:', after.description);

// Stale version must be rejected:
try {
  await db.update(r.id, {description: 'stale write'}, {version: before.version}); // old version
  console.log('UNEXPECTED: stale update succeeded');
} catch (e) {
  console.log('version conflict OK:', e instanceof DG.DomainVersionConflictError,
              e.expectedVersion, e.currentVersion);
}
```

**Pass:** first update succeeds and `version` increments; the second throws
`DG.DomainVersionConflictError`.

---

## TC5 — Soft-delete

```js
await db.delete(r.id);
console.log('gone from get:', await db.get(r.id));       // → null
console.log('gone from query:', await db.query({filter: 'name = "CRUD test"'})); // → []
```

**Pass:** after delete, `get` returns `null` and the row no longer appears in queries (soft-delete
hides it; the row is retained with `is_deleted = true` server-side).

---

## TC6 — Validation: `choices` constraint on `method`

```js
try {
  await db.insert({name: 'bad method', source: src, method: 'bogus'});
  console.log('UNEXPECTED: invalid method accepted');
} catch (e) {
  console.log('choices rejected OK:', e instanceof DG.DomainValidationError, e.message);
}
```

**Pass:** throws `DG.DomainValidationError` (server rejects a value outside the declared choices).

---

## TC7 — Validation: required fields

```js
try {
  await db.insert({name: 'no source'});                  // source is required
  console.log('UNEXPECTED: missing source accepted');
} catch (e) {
  console.log('required rejected OK:', e instanceof DG.DomainValidationError, e.message);
}
try {
  await db.insert({source: src});                        // name is required (isName)
  console.log('UNEXPECTED: missing name accepted');
} catch (e) {
  console.log('required name rejected OK:', e instanceof DG.DomainValidationError, e.message);
}
```

**Pass:** both throw `DG.DomainValidationError`.

---

## TC8 — Privacy (row mode + `defaultRowVisibility: none`) — the key check

**Goal:** a model is private to its author until shared — reproduces MyFiles semantics. Requires
two users.

**As user A** (JS console on the dev stand):
```js
const [priv] = await db.insert({name: 'private A', source: src});
console.log('A sees own row:', (await db.get(priv.id)) !== null);   // true
console.log('A row id:', priv.id);
```

**As user B** (second account / incognito, JS console):
```js
const dbB = grok.dapi.domains.table('diffstudio.model');
console.log('B query excludes A row:', (await dbB.query({})).every(m => m.name !== 'private A')); // true
console.log('B get returns null:', await dbB.get('<paste A row id>'));                            // null
```

**Back as user A — share the row with B**, then B re-checks:
```js
// A: grant View on the specific row to B (or B's group), or use the row's Sharing dialog in the UI.
await db.grant('<B user or group id>', 'View', priv.id);   // per-row share (promotes the row to an entity)
```
```js
// B: now the row is visible
console.log('B sees shared row:', (await dbB.get('<A row id>')) !== null);   // true
```

**Pass:** before sharing, B cannot see or `get` A's row; after sharing, B can. (If `db.grant`
with a row id is not the exact signature on this build, use the row's **Sharing** pane in the UI —
the behavior under test is the same.)

**Fail signals:** B sees A's unshared row in TC8 step 2 → `defaultRowVisibility: none` did not take
effect; this is the highest-value finding to catch here.

---

## TC9 — Audit trail (optional; audit is on by default)

```js
const [a] = await db.insert({name: 'audit test', source: src});
await db.update(a.id, {description: 'v2'}, {version: a.version});
console.log(await db.audit(a.id));    // → history entries (insert, then update with before/after)
await db.delete(a.id);
```

**Pass:** `audit(id)` returns the change history with before/after diffs.

---

## Cleanup (teardown — always run)

Per the repo testing rules, remove all server-side state created above:

```js
for (const m of await db.query({filter: 'name in ("CRUD test", "bad method", "no source", "private A", "audit test")'}))
  await db.delete(m.id);
```

Also delete any leftover rows from failed asserts. For the two-user test, A should delete
`private A` (or unshare + delete).

---

## Exit criterion (Stage 0 done)

| Check | Proves |
|-------|--------|
| TC0 | schema deploys server-side (row + `defaultRowVisibility: none` accepted) |
| TC1–TC5 | typed client ↔ backend: insert/get/query/update/soft-delete, versions, dayjs |
| TC4 (stale) | optimistic concurrency — the basis for Stage 1 save-in-place |
| TC6–TC7 | server-side validation (choices, required) |
| TC8 | privacy = MyFiles semantics (prerequisite #3) |

All green → proceed to **Stage 1** (wire save/load into `app.ts` behind the flag). Any failure in
TC0 or TC8 is a schema-design issue to fix before writing UI code.

## Notes

- **Error types** are `DG.DomainError` subclasses — match by class, never by message text:
  `DomainValidationError`, `DomainVersionConflictError`, `DomainNotFoundError`,
  `DomainForbiddenError`, `DomainFilterError`.
- **Datetimes** (`created_on`, `updated_on`) come back as **dayjs** objects, not strings.
- The `diffStudioEms` flag is irrelevant to these tests — they hit the data plane directly. The
  flag only gates the app's choice of DB vs files (Stage 1+).
