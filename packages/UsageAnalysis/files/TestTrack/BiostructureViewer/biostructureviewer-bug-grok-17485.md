---
feature: biostructureviewer
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [GROK-17485]
realizes: []
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

Regression guard for [GROK-17485](https://reddata.atlassian.net/browse/GROK-17485)
("BiostructureViewer - Does not maintain opened structure file"). Around
v1.21.4 the **Biostructure** viewer stopped saving its in-memory structure
data into the saved project / dashboard — opening a PDB, saving the
project, and reopening it (typically as a different user with whom the
project was shared) left the viewer empty instead of re-rendering the
original structure. The bug is fixed; this scenario exists to catch a
re-regression. This scenario asserts the BiostructureViewer side of the
contract (the structure re-renders on reopen) — project save/share
mechanics themselves are verified in the Projects section's own tests.

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

- Deferrals: Scenario 2's cross-user reopen leg is skipped when a
  second test-account fixture isn't wired up in the environment (see
  the frontmatter for the exact condition). The same-user save +
  reopen invariant in Scenario 1 is the structural regression guard;
  the cross-user leg adds permissions-traversal coverage on top of
  it and is not required for the core regression check.

- An apitest variant covering Scenario 1's same-user roundtrip
  (skipping the rendered-canvas assertion and asserting only that
  the saved project payload contains a non-empty `dataJson`) would
  be a useful additional check — filed as a follow-up, not a gap in
  this scenario.

- The existing smoke scenario (`biostructure-viewer.md`) opens
  structures but never saves a project, so it doesn't exercise this
  path. This scenario is the dedicated regression guard for that gap.

- Project save / share / reopen mechanics themselves belong to the
  Projects feature's own tests. This scenario only asserts the
  BiostructureViewer side of the contract: given a working
  save/reopen pipeline, the viewer must rehydrate from persisted
  state.
