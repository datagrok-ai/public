---
feature: projects
sub_features_covered:
  - projects.api.save
  - projects.api.get-by-id
  - projects.api.list
  - projects.api.delete
  - projects.api.count
target_layer: api-contract
priority: regression
produced_from: atlas-driven
related_bugs: []
---

# Projects — REST lifecycle contract

## Setup

1. Authenticate to a Datagrok server using a service-level token; the
   contract test must run against a real server, not a mock.
2. Pick a unique project name like `api-lifecycle-${Date.now()}`.

## Scenarios

### Scenario 1: Save → fetch → delete round-trip (regression)

Priority: regression

Steps:
1. Construct an empty `Project` with `friendlyName = <unique name>` and
   call `grok.dapi.projects.save(project)` — verify the response
   resolves and the returned `Project` carries a non-empty `id`.
2. Call `grok.dapi.projects.find(<id>)` (`GET /projects/<id>`) and
   verify the returned project's `id` and `name` match.
3. Pre-delete count: call `grok.dapi.projects.list().filter('name = "<name>"').count()`
   (`GET /projects/count?name=...`) and verify `1`.
4. Delete the project (`DELETE /projects/<id>`).
5. Repeat the count query and verify `0`.

Expected:
- Each call returns the documented HTTP status (200 / 200 / 200 / 200 / 200)
  and the entity is consistent across save, find, and count.
- After delete, `find` of the same id rejects (404) and the count goes
  to 0.

### Scenario 2: List with filter and pagination (regression)

Priority: regression

Steps:
1. Save three projects with names sharing a unique prefix
   `api-lifecycle-batch-<n>`.
2. Call `grok.dapi.projects.list().filter('name like "api-lifecycle-batch-%"').by(2).first()`
   to verify pagination first page returns 2 entries.
3. Call `grok.dapi.projects.list().filter('name like "api-lifecycle-batch-%"').page(1).by(2)`
   and verify the second page has the third entry.
4. `.count()` over the same filter returns `3`.
5. Cleanup: delete all three.

Expected:
- The list endpoint honors `name`, paging (`page`, `limit`) and the
  count endpoint returns the same total as the filter result.

### Scenario 3: Delete of non-existent id (edge)

Priority: edge

Steps:
1. Call `grok.dapi.projects.delete({id: '<random uuid>'})` for an id
   that is known not to exist on the server.
2. Capture the rejection.
3. Save then delete a real project; immediately re-delete it.
4. Capture the second rejection on the same id.

Expected:
- Both delete-of-missing calls reject with the documented status
  (404). The contract: a successfully-deleted project is gone, and
  the next delete on the same id is treated as missing — not as a
  duplicate-success.

## Notes

- **target_layer rationale:** every step is a JS API call against
  `grok.dapi.projects`; no UI surface is exercised. This is the
  `api-contract` shape per STEP D.
- **STEP B coverage map status:** at `version_of_scan: 2026-04-27`,
  the coverage map shows `projects.upload: covered`, with every other
  sub_feature labelled `partial` (incidental coverage only — the bulk
  of the catalog has no dedicated contract test). This file targets
  the partial-coverage subset most central to the lifecycle.
- **STEP C fallback path taken:** no `gap_description` was provided.
  The `projects` atlas (revision 1, `atlas_level: lightweight`) has
  `critical_paths: []` — it has not been curated for critical_paths
  yet. Additionally, no sub_feature is `uncovered` in the coverage
  map. Per STEP C-2, generator emitted the warning and fell back
  to picking the highest-interaction `partial` sub_features grouped
  into a coherent CRUD lifecycle.
- **Atlas level note:** `lightweight` — `interactions`, `edge_cases`,
  `manual_only`, `coverage_exceptions`, `help_docs`, `impl_docs` are
  all empty. Scenario quality is bounded by source-derived behavior
  only; the curated path through this surface should be added when
  Olesia / Olena enrich `projects.yaml`.
- **Deferrals:**
  - File-storage and relation sub-APIs (`projects.api.files.*`,
    `projects.api.relations.*`) — deferred. Real dependency
    cited (per A-MERIT-02 / Lattice Rule 13): file-upload contract
    requires fixture binaries already on the server, which depend
    on a curated test data set that is not yet in
    `helpers-registry.yaml`. To be picked up after Olesia adds a
    fixture-builder helper for project files.
  - Namespace endpoints (`projects.api.namespaces.*`) — deferred.
    Requires a pre-existing namespace tree with known names; the
    `helpers-registry.yaml` does not yet expose a namespace-fixture
    builder. Tracked alongside the file-storage gap.
- **No gap_description was passed**; this is the `critical_paths`-empty
  fallback path explicitly noted in STEP G.
