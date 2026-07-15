---
feature: projects
target_layer: apitest
coverage_type: regression
priority: p0
realizes_atlas: [search-list-recent]
realizes: []
pyramid_layer: integration
ui_coverage_responsibility: []
ui_coverage_delegated_to: null
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/lifecycle-api.md
migration_date: 2026-05-04
related_bugs: []
---

# Projects — REST lifecycle contract

Exercises the project REST/JS-API contract directly — save, fetch by
id, list with filters and pagination, count, and delete — independent
of any UI. Confirms that a deleted project is truly gone (not
soft-deleted-but-still-returned), that list and count results agree,
and that deleting a non-existent or already-deleted project is
rejected rather than silently succeeding.

## Setup

1. Authenticate to a Datagrok server using a service-level token; the
   contract test must run against a real server, not a mock.
2. Pick a unique project name like `api-lifecycle-${Date.now()}` to
   avoid collisions across parallel test runs.

## Scenarios

### Scenario 1: Save → fetch → delete round-trip

Single-project lifecycle: save a project, find it by id, count it,
delete it, and verify it is gone.

1. Construct an empty `Project` with `friendlyName` set to the unique
   name and call `grok.dapi.projects.save(project)`. Verify the
   response resolves and the returned `Project` carries a non-empty
   `id`.
2. Call `grok.dapi.projects.find(<id>)` (`GET /projects/<id>`) and
   verify the returned project's `id` and `name` match what was
   saved.
3. Call `grok.dapi.projects.list().filter('name = "<name>"').count()`
   (`GET /projects/count?name=...`) and verify it returns `1`.
4. Delete the project via `grok.dapi.projects.delete(project)`
   (`DELETE /projects/<id>`).
5. Repeat the count query from step 3 and verify it returns `0`.
6. Call `grok.dapi.projects.find(<id>)` again and verify the call
   rejects (the entity is gone — the contract is "deleted means
   not-findable", not "soft-deleted but still returned").

### Scenario 2: List with filter and pagination

Batch lifecycle covering list filtering, pagination, and count
agreement.

1. Save three projects with names sharing a unique prefix
   `api-lifecycle-batch-<n>` for n in {1, 2, 3}.
2. Call
   `grok.dapi.projects.list().filter('name like "api-lifecycle-batch-%"').by(2).first()`
   and verify the first page returns 2 entries.
3. Call
   `grok.dapi.projects.list().filter('name like "api-lifecycle-batch-%"').page(1).by(2)`
   and verify the second page contains the third entry.
4. Call `.count()` over the same filter and verify it returns `3` —
   the count must agree with the sum of the paginated entries.
5. Cleanup: delete all three projects and verify the count returns to
   `0`.

### Scenario 3: Delete of non-existent id

Negative-path contract for the delete endpoint.

1. Call `grok.dapi.projects.delete({id: '<random uuid>'})` for an id
   that is known not to exist on the server. Verify the call
   rejects.
2. Save a fresh project, delete it, then immediately re-call delete
   on the same id. Verify the second call also rejects.
3. The contract under test: a successfully-deleted project is gone,
   and the next delete on the same id is treated as missing — not
   as a duplicate-success.

## Notes

- **No UI surface.** Every step is a direct JS API call against
  `grok.dapi.projects`; there's nothing to delegate to a UI-driven
  scenario.
- **Self-cleaning.** Every sub-scenario deletes the projects it
  creates; its ephemeral `api-lifecycle-${Date.now()}` and
  `api-lifecycle-batch-<n>` projects don't survive past the contract
  test, so this scenario doesn't need to be cleaned up by
  `projects-ui-smoke.md`.
- **Deferrals.**
  - File-storage and relation sub-APIs (`projects.api.files.*`,
    `projects.api.relations.*`) are not covered here — a file-upload
    contract test needs fixture binaries on the server, which needs
    a curated test data set that doesn't exist yet.
  - Namespace endpoints (`projects.api.namespaces.*`) are not
    covered here either — that needs a pre-existing namespace tree
    with known names, which has no fixture builder yet.
- **Complements the per-source-class lifecycle scenarios.** This
  scenario covers the source-agnostic CRUD contract; the other
  `projects-lifecycle-*.md` scenarios cover per-source-class
  behavior (files / query / script / spaces / db_table / derived).
