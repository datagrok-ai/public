---
feature: general
target_layer: manual-only
coverage_type: regression
produced_from: split
original_path: public/packages/UsageAnalysis/files/TestTrack/General/files-cache.md
split_date: 2026-06-16
related_bugs: []
manual_only_reason: |
  Connection-level file caching — the "Cache" checkbox in the connection Edit
  dialog, the "Cache..." mapping dialog with its cron invalidation
  (`*/2 * * * *`), and the assertion that a mapping is removed together with its
  folder — is a Dart-only connection configuration surface with no public JS
  API. The folder/file CRUD lifecycle is automated in files-cache-spec.ts; these
  cache-configuration steps remain manual. (Distinct from LocalCashing/
  local-cashing-ui.md, which is about function/query result caching, not file
  storage caching.)
---

# Files cache — connection cache configuration (manual)

Manual companion to `files-cache.md`. The folder/file CRUD lifecycle (create
folder, create/write file, rename file, rename folder, delete folder) is
automated by `files-cache-spec.ts`. This file covers the connection-cache
configuration, which has no JS API.

## Preconditions

- A file connection with caching support (e.g. **Files > Demo**).

## Steps

1. **Enable caching on the connection.**
   - Open the context menu on the Demo connection → **Edit...** → tick the
     **Cache** checkbox → save.

2. **Create a cache mapping with cron invalidation.**
   - Create a folder in the connection's root to serve as the mapping target,
     e.g. `Folder cache test` (folder/file creation itself is covered by
     `files-cache-spec.ts` — here it is only the mapping target).
   - Open the context menu on the Demo connection → **Cache...** → choose that
     folder → add cron invalidation `*/2 * * * *` (every 2 minutes). The mapping
     is saved.

3. **Verify the mapping is removed together with its folder.**
   - Delete the mapped folder (folder deletion itself is covered by the spec).
   - Open the context menu → **Cache...**; the cache mapping that was created for
     that folder should now also be gone — this auto-cleanup is the
     cache-specific assertion that the spec cannot make.

> Scope note: this file does **not** re-test file/folder create / rename /
> delete — that lifecycle is automated in `files-cache-spec.ts`. Only the
> connection-cache configuration and its mapping auto-cleanup (no JS API) live
> here, so there is no overlap with the autotest.

---
{
  "order": 4
}
