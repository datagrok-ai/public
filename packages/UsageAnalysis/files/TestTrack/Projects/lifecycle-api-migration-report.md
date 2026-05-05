# Migration Report — lifecycle-api.md

First migration report for this scenario. Prior cycles produced the
`lifecycle-api.md` body but no co-located migration report (the
scenario was authored as `produced_from: atlas-driven`, not as a
TestTrack-original migration). This report is generated against
`scenario-chains/projects.yaml` revision 3 (2026-05-04) and
validates that the Stage-1 reconciliation correctly handles
atlas-driven scenarios alongside TestTrack-sourced ones in the
chain analyzer's two-mode output.

## Step mapping

This scenario has no TestTrack-original numbered steps to map (it
is `produced_from: atlas-driven`, not migrated from a manual `.md`).
The table below maps the atlas inputs that drove each migrated
section, per chain rev 3 `dependency_graph.lifecycle-api.md` and
atlas rev 11 `sub_features`.

| Original input (atlas) | Migrated location | Decision |
|------------------------|-------------------|----------|
| atlas sub_feature `projects.api.save` (`POST /projects`) | Setup steps 1–2 (precondition) + Scenario 1 step 1 + Scenario 2 step 1 + Scenario 3 step 2 | preserved (covered by every save call across the three sub-scenarios) |
| atlas sub_feature `projects.api.get-by-id` (`GET /projects/<id>`) | Scenario 1 step 2 + Scenario 1 step 6 (negative — verify gone after delete) | preserved (positive + negative contract verified) |
| atlas sub_feature `projects.api.list` (`GET /projects?<filters>`) | Scenario 2 steps 2–3 (filter + pagination) | preserved |
| atlas sub_feature `projects.api.count` (`GET /projects/count?<filters>`) | Scenario 1 steps 3, 5 + Scenario 2 step 4 + Scenario 2 step 5 | preserved (count-agreement invariant exercised on every save/delete edge) |
| atlas sub_feature `projects.api.delete` (`DELETE /projects/<id>`) | Scenario 1 step 4 + Scenario 2 step 5 + Scenario 3 steps 1–2 | preserved (positive + idempotency-failure contract verified) |
| atlas critical_path `search-list-recent` (priority p3) | Scenario 2 (filter + count + pagination) | preserved as a sub-scenario |
| atlas critical_path `upload-save-reopen-golden` REST projection | Scenario 1 (save + find + delete round-trip) | preserved as the primary sub-scenario; "reopen" maps to `find` at the api-contract layer (no shell.open at this layer) |
| atlas sub_feature `projects.api.files.*` (file storage) | (deferred) | deferred — fixture-builder helper for project files not yet in `helpers-registry.yaml` |
| atlas sub_feature `projects.api.relations.*` | (deferred) | deferred — requires curated test data set not yet in `helpers-registry.yaml` |
| atlas sub_feature `projects.api.namespaces.*` | (deferred) | deferred — requires namespace-fixture builder not yet in `helpers-registry.yaml` |
| atlas sub_feature `projects.api.search`, `projects.api.search-namespace`, `projects.api.recent`, `projects.api.tags`, `projects.api.guest-session-token` | (not in this scenario) | (moved to atlas — left for follow-up atlas-driven cycles; flagged in `## Edge cases` below) |

## Decisions

- **Why this `target_layer`:** chose `api-contract` per chain rev 3
  `output_plan` entry for lifecycle-api.md
  (`target_layer: api-contract`, `strategy: simple`). Per
  migration-prompt.md "Migrated scenario shape" decision tree, rule
  #3 applies (scenario is purely REST-call shape — every step is a
  `grok.dapi.projects.*` call with no UI surface). No persistence
  across navigations, no viewers, no widgets.
- **Why this `strategy`:** `simple` per chain rev 3
  `output_plan.lifecycle-api.md.strategy: simple`. Three stand-alone
  sub-scenarios; no state shared across them; no chained_tests, no
  data_driven (the three sub-scenarios are not parameterised
  variants of one shape — they exercise distinct contracts).
- **Schema fields populated from chain rev 3 (Edit 2 + atlas-driven
  classification):**
  - `pyramid_layer: integration` — sourced from chain YAML
    `dependency_graph.lifecycle-api.md.pyramid_layer`. Rationale
    cited in chain notes (Rule 4 fallthrough — Rule 1 inapplicable
    because no UI surface; Rules 2/3/5 inapplicable; Rule 4 fits
    because the REST contract is source-agnostic by construction).
  - `ui_coverage_responsibility: []` — sourced from chain YAML
    same node; correctly empty because no UI flows are driven.
  - `ui_coverage_delegated_to: null` — sourced from chain YAML
    same node; by-design no UI to delegate (per chain
    `ui_coverage_plan.delegated_scenarios` rationale: "UI coverage
    N/A by design — api-contract target_layer").
  - `produced_from: atlas-driven` — RETAINED from prior cycle
    (this is the key validation property — the chain analyzer
    rev 3 must continue to recognize and preserve atlas-driven
    origin alongside `migrated` originating from TestTrack `.md`
    files).
- **UI coverage ownership:** N/A by design. Per migration-prompt.md
  "UI delegation in SCOPE_REDUCTION proposals", absence of any UI
  flow means no UI-coverage-gap citation is required and no
  delegation language applies. Chain rev 3
  `ui_coverage_plan.delegated_scenarios` carries the explicit
  `delegated_to: null` entry for this scenario with the
  by-design rationale.
- **Sibling tests consulted (per migration-prompt.md "Sibling
  test consultation"):**
  - `existing-test-index.yaml` summary: `api-contract: 0` tests
    in the index as of 2026-05-04. lifecycle-api.md is the
    pioneer scenario for this layer; no sibling api-contract spec
    exists yet to align house style with.
  - `public/packages/UsageAnalysis/files/TestTrack/Projects/upload-project-spec.ts`
    (target_layer: playwright, style: c-playwright) — consulted
    for `grok.dapi.projects.*` usage patterns. Confirms the
    canonical `grok.dapi.projects.save(...)` shape; downstream
    Automator should align the api-contract spec at this call
    surface even though the broader test architecture differs
    (no Playwright `page` object).
  - `public/packages/ApiSamples/scripts/dapi/` (referenced by
    project CLAUDE.md as the canonical samples for `grok.dapi`
    usage) — consulted as the conceptual sibling for
    api-contract calls.
- **Helpers reused (per migration-prompt.md "Helpers
  discipline"):**
  - None directly. `helpers-registry.yaml` does not currently
    register any api-contract / Dart-style project helpers
    (every existing project helper is `grok_test_layer` or
    `playwright_layer` — see registry entry for
    `uploadProject` at line 4, `style_compatibility: [a, b]`).
  - **Candidate helper (not invented; flagged for registry):**
    `helpers.api.projects.saveAndCleanup({prefix}, fn)` — a
    fixture-style wrapper that allocates a unique-name project,
    runs the test body, and guarantees deletion in a `finally`
    block. Used implicitly across all three sub-scenarios.
    Propose adding to `helpers-registry.yaml :: api_contract_layer`
    if Block 2 of the orchestrator chain produces a second
    api-contract scenario that would reuse it.
  - **Candidate helper (not invented; flagged for registry):**
    `helpers.api.projects.assertGoneById(id)` — wraps the
    "find rejects + count is 0" double-check used in
    Scenario 1 step 5/6 and Scenario 3.
- **Bug library consulted:** yes —
  `bug-library/projects.yaml` rev 2 (8 curated bugs as of
  2026-04-28). `related_bugs: []` retained per prior cycle
  (consistent with the existing frontmatter on disk). Rationale:
  the curated bugs all have UI/multi-user repro paths
  (Save Copy with Link, rename + sync, share + un-shared deps,
  Spaces + datasync) that lie outside the source-agnostic REST
  contract this scenario exercises. Specifically:
  - GROK-19750 (`projects.api.save`, `projects.add-relation`):
    repro requires Save Copy dialog; out of scope.
  - GROK-19212 (`projects.api.save`, `projects.api.files.sync`,
    `projects.api.relations.list`): repro requires data sync +
    rename; `projects.api.files.sync` is in the deferred set.
  - GROK-19103 (`projects.api.save`, `projects.add-relation`):
    repro requires UI Join → Save flow; out of scope.
  - GROK-19403 (`projects.api.get-by-id`, `projects.add-relation`):
    repro requires share + recipient open; multi-user path out
    of scope.
  - GROK-18345 (`projects.api.get-by-id`,
    `projects.api.files.sync`, `projects.api.namespaces`): repro
    requires Spaces + datasync + share; out of scope.
  - GROK-19728 (`projects.api.get-by-id`): repro requires share
    + view-and-use access edit; multi-user path out of scope.
  - github-3550 (`projects.api.relations.list`,
    `projects.add-relation`): repro requires Query rename;
    `projects.api.relations.list` deferred.
- **Cross-cutting bug citations from chain YAML rev 3
  `bug_focused_candidates` (per migration-prompt.md
  "Cross-cutting bug citations from chain YAML" —
  RECOMMENDED):**
  - Bug GROK-19403 has cross-cutting candidate spec — see chain
    `bug_focused_candidates` entry, proposed_spec:
    `projects-grok-19403-spec.ts` (spans
    `share-project.md:Step 4`, `complex.md:Step 12`,
    `projects-copy-clone.md:Step 5`). Chain rationale notes
    intersection with `lifecycle-api.md` via
    `projects.api.get-by-id` sub_feature; this scenario
    exercises the positive-path get-by-id contract (Scenario 1
    step 2) and the negative-path (Scenario 1 step 6 +
    Scenario 3) but NOT the share + un-shared-deps repro path.
  - Bug GROK-18345 has cross-cutting candidate spec — see chain
    `bug_focused_candidates` entry, proposed_spec:
    `projects-grok-18345-spec.ts` (spans
    `uploading.md:Step 7`, `complex.md:Step 12`). Chain
    rationale notes intersection with `lifecycle-api.md` via
    `projects.api.get-by-id`; this scenario does not cover
    Spaces + datasync (deferred items).
  - Bug GROK-19728 has cross-cutting candidate spec — see chain
    `bug_focused_candidates` entry, proposed_spec:
    `projects-grok-19728-spec.ts` (spans `complex.md:Step 12`).
    Chain rationale notes the stricter intersection yields
    exactly 1 affecting scenario (lifecycle-api.md via
    `projects.api.get-by-id`); this scenario exercises only the
    happy-path get-by-id, not the failure-state repro.
  - F-BUG-COVERAGE-01 at section-complete is the authoritative
    gate; the citations above are early-visibility only.
- **Decision log queried (D10 active query):** yes — entries
  filtered by `feature: projects` AND `scenario: lifecycle-api`
  OR `lifecycle`:
  - `mig-2026-04-30-priority-enum-drift-discovery` (resolution
    applied 2026-04-30): canonical scenario test-kind field is
    `coverage_type` (values `smoke|regression|edge|perf`) per
    A-STRUCT-MECH-06 (Phase 0 Block 11). Honored in this
    re-migration: frontmatter uses `coverage_type: regression`,
    NOT `priority: <severity>`. The on-disk frontmatter from
    the prior cycle already used `coverage_type: regression`
    correctly; preserved.
  - `mig-2026-04-30-custom-creation-scripts-migration` (line 646
    of decision-log): identifies the slot in the chain that
    lifecycle-api.md was added to fill after the renamed
    custom-creation-scripts.md → -ui.md left a vacancy. Chain
    rev 3 explicitly notes lifecycle-api.md is "NOT a 1:1
    replacement" — it has different `produced_from`,
    different `target_layer`, different `pyramid_layer`.
  - `wave-3-d3-classification-d4-verification` (referenced by
    chain rev 3 commentary): the 2026-05-01 file rename of
    custom-creation-scripts.md → custom-creation-scripts-ui.md
    that vacated the slot. Documented for audit traceability.
  - No `failed_attempts` matching this scenario —
    `approach_tried` table empty for lifecycle-api.md.

## Opt-outs (SCOPE_REDUCTION proposals)

(none) — no UI step was reduced to a JS API equivalent because
this scenario has no UI to begin with. atlas-driven origin means
the scenario was generated against the REST surface directly; no
SCOPE_REDUCTION semantics apply. Per migration-prompt.md
"UI delegation in SCOPE_REDUCTION proposals", absence of
SCOPE_REDUCTION means no UI-coverage-gap citation is required.

## Deferred items (NOT opt-outs)

- **`projects.api.files.*` sub-APIs (file-storage contract):**
  deferred. Real prerequisite cited: file-upload contract
  requires fixture binaries already on the server, which depend
  on a curated test data set that is not yet in
  `helpers-registry.yaml :: api_contract_layer`. To be picked
  up after a fixture-builder helper for project files is added.
- **`projects.api.relations.*` sub-APIs (relations contract):**
  deferred. Real prerequisite cited: relation save/list/delete
  requires entity-id fixtures (table info, query, script) that
  are not yet exposed via api-contract helpers. Tracked
  alongside the file-storage gap.
- **`projects.api.namespaces.*` sub-APIs (namespace contract):**
  deferred. Real prerequisite cited: requires a pre-existing
  namespace tree with known names; the `helpers-registry.yaml`
  does not yet expose a namespace-fixture builder.

## Edge cases

- **Concurrent name collision:** the
  `api-lifecycle-${Date.now()}` naming pattern in Setup step 2
  guards against two parallel test runs colliding on the same
  name (millisecond-resolution timestamp). Implicit edge case;
  preserved in Setup. If parallelism > N tests/ms is observed,
  add a per-process suffix — flag for atlas curator to capture
  in `edge_cases` array.
- **Server returns 404 vs rejects:** Scenario 3 expects the
  delete-of-missing call to *reject*, deferring HTTP-status-code
  specifics to the JS API client implementation (it may surface
  as a thrown error, a rejected Promise, or an HTTP 404 status
  depending on the dapi client layer). The contract under test
  is "next delete on the same id is treated as missing" — not
  "exact HTTP status code". Flagged for atlas curator: should
  the api-contract atlas note record the canonical error shape
  for `projects.api.delete` of missing id? Proposed atlas
  addition.
- **List ordering vs pagination consistency:** Scenario 2 step 2
  retrieves 2 of 3 entries via `.by(2).first()` and step 3
  retrieves the third via `.page(1).by(2)`. Contract assumes
  the underlying ordering is stable across the two calls
  (otherwise the same entry could appear on both pages). Implicit
  edge case; preserved as a verification of the count-agreement
  invariant in step 4 (count = 3 must equal 2-on-first-page +
  1-on-second-page). Proposed atlas addition: explicit
  ordering-stability invariant on `projects.api.list` + paging
  endpoints.
- **Read after delete (entity gone vs soft-deleted):** Scenario 1
  step 6 verifies `find(<deleted_id>)` rejects, asserting the
  contract is "deleted means not-findable". This is implicit in
  the atlas — `projects.api.delete` description says "delete a
  project and reply 'ok'" but does not state the post-condition
  on `get-by-id`. Proposed atlas addition: explicit
  not-findable-after-delete invariant on `projects.api.delete`.
- **Atlas sub_features outside this scenario:**
  `projects.api.search`, `projects.api.search-namespace`,
  `projects.api.recent`, `projects.api.tags`,
  `projects.api.guest-session-token` are not exercised by this
  scenario. They are NOT atlas `manual_only` (atlas rev 11
  `manual_only` block lists only drag-drop / move sub_features),
  so they are valid api-contract candidates. Flagged for atlas
  curator as candidates for follow-up atlas-driven cycles —
  proactive_lifecycle_specs candidates per
  Block 7 / Block 11 retro:
  - `projects-api-search-spec.ts` (search + search-namespace
    cross-namespace scoping; chain
    `unresolved_ambiguities[6]` partially relates).
  - `projects-api-recent-spec.ts` (recent + page/limit
    contract).
  - `projects-api-tags-spec.ts` (tags filter + tags arg).
  - `projects-api-guest-session-token-spec.ts` (guest token
    contract; may have multi-user prerequisite).

## Unresolved ambiguities

- **`coverage_type` choice (smoke vs regression):** prior cycle
  on-disk frontmatter declared `coverage_type: regression`.
  Re-migration preserves the value. The scenario could
  arguably be `smoke` (golden-path REST round-trip) but
  `regression` better reflects its role as the catch-all
  api-contract scenario for the source-agnostic CRUD surface,
  not a primary smoke (the chain's smoke is upload-project.md
  per `ui_coverage_plan.smoke_scenario`). Flag for retro if
  Gate D challenges.
- **`related_bugs: []` vs explicit list:** the prior on-disk
  frontmatter carried `related_bugs: []`. Chain rev 3
  `bug_focused_candidates` lists multiple bugs whose
  `bug.affects` intersects this scenario's
  `sub_features_covered` (GROK-19403, GROK-18345, GROK-19728
  via `projects.api.get-by-id`). The `bug_focused_candidates`
  rationale acknowledges these are cross-cutting bugs whose
  repro paths lie OUTSIDE this scenario (UI / multi-user /
  Spaces flows). Decision: preserve `related_bugs: []` —
  per migration-prompt.md "related_bugs", the field captures
  bugs whose `affects` intersects sub_features AND whose
  repro path is exercised in this scenario; merely sharing a
  sub_feature is insufficient. Flagged for retro: should the
  semantic of `related_bugs` be "affects intersects" (broader)
  or "repro covered" (stricter)?
- **Atlas-driven origin preserved across re-migration:** this
  re-migration validates that chain analyzer rev 3 correctly
  preserves `produced_from: atlas-driven` (the file is NOT
  treated as a TestTrack-original migration — its origin
  marker is retained). If a future migrate cycle silently
  rewrites this to `produced_from: migrated`, that is a
  regression in the chain analyzer's two-mode logic. No
  rewrite occurred this cycle — atlas-driven retained.
- **Atlas REST sub_features not yet covered by any scenario:**
  `projects.api.search`, `projects.api.search-namespace`,
  `projects.api.recent`, `projects.api.tags`,
  `projects.api.entity-paths`, `projects.api.get-entity-project`,
  `projects.api.guest-session-token`, plus the deferred
  `projects.api.files.*`, `projects.api.relations.*`,
  `projects.api.namespaces.*` families. Atlas rev 11 lists
  these as sub_features but no scenario in the Projects chain
  currently exercises them at the api-contract layer. Flagged
  for `proactive_lifecycle_specs` follow-up additions per
  source-class-orthogonal API surface (these are
  source-agnostic REST endpoints, so the existing
  `proactive_lifecycle_specs[]` source-class × dep-op cells
  do not naturally cover them).
