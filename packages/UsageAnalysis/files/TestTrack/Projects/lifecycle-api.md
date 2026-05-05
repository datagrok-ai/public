---
feature: projects
sub_features_covered:
  - projects.api.save
  - projects.api.get-by-id
  - projects.api.list
  - projects.api.delete
  - projects.api.count
target_layer: api-contract
coverage_type: regression
pyramid_layer: integration
ui_coverage_responsibility: []
ui_coverage_delegated_to: null
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/lifecycle-api.md
migration_date: 2026-05-04
migration_report: lifecycle-api-migration-report.md
related_bugs: []
---

# Projects — REST lifecycle contract

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

- **Origin (atlas-driven):** this scenario is `produced_from:
  atlas-driven`. Its content is derived from the `projects` feature
  atlas (rev 11) `sub_features` for the REST surface — `projects.api.save`,
  `projects.api.get-by-id`, `projects.api.list`, `projects.api.delete`,
  `projects.api.count` — and from the `search-list-recent` critical
  path. The scenario is NOT a migration of a TestTrack-original `.md`;
  it was generated to fill the api-contract layer that the chain's
  TestTrack scenarios do not cover (every other scenario in the chain
  is `target_layer: playwright`).
- **target_layer rationale:** every step is a JS API call against
  `grok.dapi.projects`; no UI surface is exercised. This is the
  `api-contract` shape per chain YAML rev 3 `output_plan` (target_layer:
  api-contract, strategy: simple).
- **pyramid_layer: integration (Rule 4 fallthrough):** Rule 1
  (ui-smoke) inapplicable — no UI surface to drive. Rule 2
  inapplicable — no matrix. Rule 3 inapplicable — no
  `related_bugs` frontmatter. Rule 5 inapplicable — not
  `manual_only`. Falls through to Rule 4 (source-agnostic op): the
  REST contract holds for any project regardless of its underlying
  source class.
- **No UI delegation:** `ui_coverage_responsibility: []` and
  `ui_coverage_delegated_to: null` — there is no UI by design, so
  there is nothing to delegate. Per chain rev 3
  `ui_coverage_plan.delegated_scenarios` entry for lifecycle-api.md,
  rationale recorded: "UI coverage N/A by design — api-contract
  target_layer".
- **Self-cleaning:** every sub-scenario deletes the projects it
  creates. lifecycle-api.md is intentionally NOT in
  `deleting.md.depends_on` (chain rev 3) — its ephemeral
  `api-lifecycle-${Date.now()}` and `api-lifecycle-batch-<n>`
  projects do not survive past the contract test.
- **Deferrals (cite real prerequisites):**
  - File-storage and relation sub-APIs (`projects.api.files.*`,
    `projects.api.relations.*`) — deferred. Real dependency
    cited: file-upload contract requires fixture binaries already
    on the server, which depend on a curated test data set that
    is not yet in `helpers-registry.yaml`. To be picked up after
    a fixture-builder helper for project files is added.
  - Namespace endpoints (`projects.api.namespaces.*`) — deferred.
    Requires a pre-existing namespace tree with known names; the
    `helpers-registry.yaml` does not yet expose a namespace-fixture
    builder. Tracked alongside the file-storage gap.
- **Coverage gap candidates for atlas-driven follow-up:** chain rev
  3 `proactive_lifecycle_specs` already enumerates source-class ×
  dep-op cells (files / query / script / spaces / db_table /
  derived). The REST surface in this scenario complements those
  per-source-class specs — it covers the source-agnostic CRUD
  contract, not the per-source lifecycle.
