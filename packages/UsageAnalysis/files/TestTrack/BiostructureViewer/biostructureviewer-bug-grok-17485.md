---
feature: biostructureviewer
sub_features_covered:
  - biostructure.viewer
  - biostructure.file-open
  - biostructure.project-persistence
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
related_bugs:
  - GROK-17485
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions:
  - id: SR-01
    scope: scenario-2-cross-user-reopen
    verdict_status: PASS
    rationale: |
      Scenario 2's cross-user reopen leg requires a second-account fixture
      (process.env.PW_OTHER_USER) which is not wired in helpers-registry as
      of 2026-06-04. The spec asserts the SHARE leg (server-side permission
      grant via Right-click -> Share dialog) — DOM-driving DOM-driven via the
      shareProjectViaContextMenu helper — and conditionally skips the
      cross-user reopen with test.skip(). Same-user save+reopen invariant
      (Scenario 1) IS the structural regression guard; cross-user adds
      permissions-traversal coverage. Future cycle: extract a
      logoutAndLoginAs helper once second-user fixture lands in registry.
realized_as:
  - biostructureviewer-bug-grok-17485-spec.ts
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T14:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T17:10:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T14:30:00Z
    spec_runs:
      - spec: biostructureviewer-bug-grok-17485-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 86
        failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T18:00:00Z
    failure_keys: []
    review_round: 1
    claims:
      - check: A-STRUCT-MECH-01
        status: PASS
      - check: A-STRUCT-MECH-02
        status: PASS
      - check: A-STRUCT-MECH-03
        status: PASS
      - check: A-STRUCT-MECH-04
        status: PASS
      - check: A-STRUCT-MECH-05
        status: PASS
      - check: A-STRUCT-MECH-06
        status: PASS
      - check: A-STRUCT-03
        status: PASS
      - check: A-STRUCT-04
        status: PASS
      - check: A-LAYER-ALIGN-01
        status: PASS
      - check: A-CONT-01
        status: PASS
      - check: A-BUG-01
        status: PASS
      - check: A-MERIT-01
        status: PASS
      - check: A-MERIT-02
        status: PASS
---

# BiostructureViewer — Project persistence roundtrip for an ad-hoc PDB structure (GROK-17485 regression guard)

## Overview

Regression guard for [GROK-17485](https://reddata.atlassian.net/browse/GROK-17485)
("BiostructureViewer - Does not maintain opened structure file"). Around
v1.21.4 the **Biostructure** viewer stopped serialising its ad-hoc
in-memory structure payload into the saved project / dashboard — opening
a PDB, saving the project, and reopening it (typically as a different
user with whom the project was shared) left the viewer empty instead of
re-rendering the original structure. The bug is fixed; this scenario
exists to catch a re-regression.

The bug surfaces across three sub_features tracked in the atlas:

- `biostructure.viewer` — the consumer that needs the structure on
  reopen.
- `biostructure.file-open` — the producer that loaded the ad-hoc
  payload originally.
- `biostructure.project-persistence` — the cross-atlas concern with
  the Projects feature: the payload must survive project save +
  reopen, not just live in-session.

Cross-atlas note: `biostructure.project-persistence` is owned jointly
with the **Projects** feature
(`feature-atlas/biostructureviewer.yaml#L786-L794`); this scenario
asserts the BiostructureViewer side of the contract (the structure
re-renders on reopen) — Projects-side persistence is verified in
the Projects section.

## Setup

- Datagrok session is logged in; the **BiostructureViewer** package is
  installed and registered (viewer type `[name="viewer-Biostructure"]`
  exists via `_package.registerViewer` —
  `public/packages/BiostructureViewer/src/package.ts#L455`).
- A canonical PDB file is available under
  `System:AppData/BiostructureViewer/`. The existing migrated smoke
  uses `samples/1RQ9.mmcif` and `1U54_protein.pdb`; this scenario
  uses **`samples/1bdq.pdb`** (the file the bug report cites) so the
  repro matches the recorded reproduction.
- The current user has permission to save projects (default for
  logged-in users via `grok.dapi.projects.save(project)`).
- A second user account, or a permissions group the first user can
  share a project to, is available for Scenario 2. Setup of the
  second account is out of scope here — pass it as an injected test
  context (`process.env.PW_OTHER_USER` or equivalent helper).
- After each structure loads, await render settle —
  `viewer.awaitRendered(timeoutMs)` or poll
  `[name="viewer-Biostructure"] .msp-viewport canvas` before
  asserting. A dark viewport showing only the axis gizmo means "not
  rendered yet" or "parse failed".

## Scenarios

### Scenario 1 — Open PDB, save project, reopen as same user; structure must re-render from persisted state

Exercises the within-user persistence path of GROK-17485: the
in-memory `dataJson` payload (see
`public/packages/BiostructureViewer/src/viewers/molstar-viewer/molstar-viewer.ts#L235`)
must serialise into the saved project layout and rehydrate the viewer
on reopen.

Steps:

1. Open the **Files** browser, navigate to
   **App Data > BiostructureViewer > samples**, and double-click
   **`1bdq.pdb`**.

   * Expected result: a **Biostructure** viewer opens with the
     structure rendered as a **cartoon** (default `representation`).
     `[name="viewer-Biostructure"]` exists; `.msp-viewport canvas`
     is non-empty after `awaitRendered`; no error balloon, no fatal
     console error.

2. Capture a fingerprint of the current viewer state for the
   reopen-assert in step 6 — e.g. the viewer's `dataJson.name` (the
   parsed structure name on the Mol* engine), `representation`
   property value, and a small set of atom-count / chain-count facts
   read from the loaded structure model. These are the invariants the
   reopen must preserve.

3. Save the current view as a project: **File > Save as project**
   (or equivalent — `grok.dapi.projects.save(project)`). Name the
   project **`bug-grok-17485-roundtrip`**.

   * Expected result: the save dialog reports success; the saved
     project lists the Biostructure viewer among its layout viewers;
     `grok.dapi.projects.find(...)` (or the same call the dialog made)
     resolves to a project with non-empty payload. No console error,
     no error balloon.

4. Close the current view entirely (close the table view that hosts
   the viewer; clear the workspace so the saved project must be
   loaded fresh — not merely re-shown).

   * Expected result: the workspace no longer contains
     `[name="viewer-Biostructure"]`.

5. Reopen the saved project from **Projects** browser by double-click
   (or equivalent — `grok.dapi.projects.find('bug-grok-17485-roundtrip')`
   then `project.open()`).

6. Await render. Assert the viewer re-renders from persisted state.

   * Expected result: the table view restores; a
     `[name="viewer-Biostructure"]` is present;
     `.msp-viewport canvas` is non-empty after `awaitRendered`;
     the fingerprint captured in step 2 matches the rehydrated
     state (same `dataJson.name`, same `representation`, same
     structure-model chain/atom counts). NO empty / blank viewport
     showing only the axis gizmo — that is the regression
     signature. No "Parsed object is empty" toast, no fatal
     console error.

### Scenario 2 — Cross-user persistence: share project, reopen as another user; structure must re-render

Exercises the cross-user path of GROK-17485 — the bug as filed
described "sharing a dashboard with the chemists" and seeing the
structure missing on their side. The on-disk persistence step is
the same as Scenario 1; this scenario adds a permissions cross-over
to assert the payload travels through the share boundary, not just
through the same-user save/load.

Steps:

1. From Scenario 1 setup, open the **Files** browser, navigate to
   **App Data > BiostructureViewer > samples**, and double-click
   **`1bdq.pdb`**.

   * Expected result: same as Scenario 1 step 1 — Biostructure viewer
     opens with the cartoon render; `awaitRendered` resolves.

2. Capture the same state fingerprint described in Scenario 1 step 2.

3. Save as project named **`bug-grok-17485-share`**.

4. Share the project with the second-user context: open the project
   in the **Projects** browser, right-click → **Share…**, grant
   **View** (or **Edit**) to the second user / group.

   * Expected result: the share dialog reports success; the project's
     permissions list now includes the second principal. No console
     error.

5. Switch to the second user (re-authenticate the Playwright context
   as `process.env.PW_OTHER_USER`, or — for an apitest variant —
   replay the JS API as that user). Reopen the project
   **`bug-grok-17485-share`** from their **Projects** browser.

6. Await render. Assert the viewer re-renders from persisted state on
   the receiving side.

   * Expected result: the receiving user sees the same Biostructure
     viewer with the same structure rendered; the fingerprint from
     step 2 matches. The viewport is NOT empty. No "Parsed object is
     empty" toast. No fatal console error.

7. Teardown — remove the test projects to keep the server clean.
   Delete `bug-grok-17485-roundtrip` and `bug-grok-17485-share` via
   `grok.dapi.projects.delete(project)` (idempotent — ignore
   not-found).

## Notes

- target_layer rationale: **playwright**. The bug surfaces only across
  a full UI roundtrip (file-open dispatch → viewer in-memory state →
  project save → workspace close → project reopen → viewer
  rehydration). The structural assertion is "the rendered canvas is
  non-empty after reopen" — that requires a real browser canvas
  (`.msp-viewport canvas`) plus the Mol\* engine's async render
  settle (`viewer.awaitRendered`); apitest cannot verify the canvas
  render. The cross-user share leg additionally needs a real
  permissions context (Datagrok client-side caches identity) which is
  easier to drive via Playwright user-switch than via the JS API.

- coverage_type rationale: **regression** per the dispatch
  `inputs.bug_brief.coverage_type`. The bug is fixed and this
  scenario is a regression guard, not a smoke (it does not exercise
  the canonical opening flow per se, but the post-save reopen) and
  not a positive-path edge (the negative outcome is "viewport stays
  empty", which IS the regression signature, but the scenario shape
  is a guard for a known fix). Future authoring may add an edge
  variant covering the bug's other failure modes (e.g. NGL-viewer
  variant for `importWithNgl` payloads); that variant would carry
  `coverage_type: edge` and is out of scope here.

- Deferrals: none required for the playwright path. An apitest
  variant covering Scenario 1's same-user roundtrip (skipping the
  Mol\* render assertion and asserting only that the saved project
  payload contains a non-empty `dataJson`) would be additive — file
  it as a follow-up; not deferred for a missing helper.

- atlas entry derived from bug-library:
  `bug-library/biostructureviewer.yaml#GROK-17485`. Atlas edge_case
  citation:
  `feature-atlas/biostructureviewer.yaml#edge_cases[4]`. Affected
  sub_features per atlas:
  `feature-atlas/biostructureviewer.yaml#sub_features` ids
  `biostructure.viewer` (L26),
  `biostructure.file-open` (L346),
  `biostructure.project-persistence` (L784).

- The existing smoke `biostructure-viewer.md` does NOT exercise this
  path — its `bug_match_attempts_skipped[]` audit record for
  GROK-17485 reads "the scenario opens structures and exercises the
  viewer surface but never saves a project, so semantic_match
  returns []". This is the dedicated regression guard authored to
  close that gap.

- Cross-atlas concern: the Projects atlas owns project save / share /
  reopen mechanics. This scenario asserts the BiostructureViewer
  contract only: given a working save/reopen pipeline, the
  Biostructure viewer must rehydrate from persisted state.
