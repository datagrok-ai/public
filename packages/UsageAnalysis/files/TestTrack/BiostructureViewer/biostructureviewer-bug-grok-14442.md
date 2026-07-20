---
feature: biostructureviewer
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [GROK-14442]
realizes: []
produced_from: atlas-driven
related_bugs:
  - GROK-14442
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - biostructureviewer-bug-grok-14442-spec.ts
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T18:15:00Z
    failure_keys: []
    review_round: 1
  f:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T14:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-automate-01
    timestamp: 2026-06-04T19:05:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-automate-01
    timestamp: 2026-06-04T21:33:00Z
    spec_runs:
      - spec: biostructureviewer-bug-grok-14442-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 75
        failure_keys: []
---

# BiostructureViewer — File-handler search disambiguation between `.pdb` and `.pdbqt` (GROK-14442 regression guard)

Regression guard for [GROK-14442](https://reddata.atlassian.net/browse/GROK-14442)
("Fix file-handler search"). The bug surfaced as a file-handler-search
collision: opening a plain `.pdb` file from the Files browser routed to
the AutoDock-pose importer (`importPdbqt`) instead of the regular
Biostructure importer (`importPdb`). The two handlers are declared side
by side in the package (`importPdb` covers `mmcif, cifCore, pdb, gro`;
`importPdbqt` covers `pdbqt`); the search-by-extension logic was matching
the longer prefix (`pdbqt`) even when the user opened a file with the
shorter extension (`pdb`). Because `importPdbqt` expects AutoDock pose
records, a plain `.pdb` routed to it either failed to parse or displayed
under the wrong UI, masking the normal "open a PDB and see it in the
Biostructure viewer" happy path entirely.

## Setup

- Datagrok session is logged in; the **BiostructureViewer** package is
  installed and registered. The two file handlers exist as registered
  functions:
  `BiostructureViewer:importPdb`
  (`public/packages/BiostructureViewer/src/package.ts#L142`,
  extensions `mmcif, cifCore, pdb, gro`) and
  `BiostructureViewer:importPdbqt`
  (`public/packages/BiostructureViewer/src/package.ts#L177`,
  extension `pdbqt`).
- Two structure fixtures are available under
  `System:AppData/BiostructureViewer/samples/`:
  - `1bdq.pdb` — a plain `.pdb` file (the file used in the bug
    report's reproduction).
  - A `.pdbqt` AutoDock-pose fixture under the same `samples/`
    directory (any package-shipped `.pdbqt`, e.g. a docking output).
    If no in-package `.pdbqt` is shipped, fall back to the AutoDock
    sample bundled with the Docking package (search for `*.pdbqt`
    under `System:AppData/Docking/`).
- The Mol* (Biostructure) viewer is registered under
  `[name="viewer-Biostructure"]`
  (`public/packages/BiostructureViewer/src/package.ts#L455`); after a
  successful importPdb route, the viewer mounts and the structure
  renders within `viewer.awaitRendered(timeoutMs)`. After a
  successful importPdbqt route, the AutoDock-pose UI assembled by
  `importPdbqtUI` opens (not the Mol* Biostructure viewer).
- For the routing assertions, the load-bearing question is "which
  registered function fired on file double-click?". Two assertion
  paths are equivalent and both are acceptable:
  - **Function-call spy.** Wrap or monitor `grok.functions.call`
    (or read the platform's function-call audit log) and capture the
    fully-qualified name of the function invoked by the double-click.
    The expected exact-string match is
    `BiostructureViewer:importPdb` for `.pdb` and
    `BiostructureViewer:importPdbqt` for `.pdbqt`. NEVER the inverse.
  - **Resulting-UI inference.** Where the function-call spy is not
    available, infer from the UI that appears: a Mol* viewport
    (`[name="viewer-Biostructure"]` + `.msp-viewport canvas`
    non-empty) implies `importPdb` fired; the AutoDock-pose UI
    (assembled by `importPdbqtUI`) implies `importPdbqt` fired.
    A `.pdb` that opens into the AutoDock-pose UI is the bug
    regression signature; a `.pdbqt` that opens into a Mol* viewport
    is the inverse regression.
- Await render / route settle before asserting: for the Mol* side,
  `viewer.awaitRendered(timeoutMs)`; for the importPdbqt side, await
  the AutoDock-pose UI mount.

## Scenarios

### Scenario 1 — Double-click `.pdb` MUST route to `BiostructureViewer:importPdb` (NOT `importPdbqt`)

Exercises the regressed direction of the bug. The bug report's exact
reproduction was: install BiostructureViewer, open
`samples/1bdq.pdb` by double-click, observe `importPdbqt()` firing
instead of `importPdb()`. This scenario asserts the fix: the dispatcher
MUST route by exact extension, so `.pdb` MUST fire `importPdb`.

Steps:

1. Open the **Files** browser via the left sidebar (Browse → Files).
   Navigate to
   **App Data > BiostructureViewer > samples**.

   * Expected result: the Files browser lists at least one `.pdb`
     fixture (`1bdq.pdb`); no error balloon, no fatal console error.

2. Begin function-call monitoring (capture the fully-qualified name
   of any `BiostructureViewer:*` function invoked over the next
   action). The assertion path documented in Setup applies — either
   wrap `grok.functions.call` to record the function name, or be
   ready to infer from the resulting UI.

   * Expected result: monitoring is active; no captured calls yet.

3. **Double-click `1bdq.pdb`.** This is the bug's reproduction
   action verbatim.

   * Expected result: exactly one
     `BiostructureViewer:importPdb` call is captured. NO
     `BiostructureViewer:importPdbqt` call is captured.
     **Regression signature**: if
     `BiostructureViewer:importPdbqt` fires instead, the test
     FAILS with diagnostic "GROK-14442 file-handler search
     regressed: .pdb routed to importPdbqt".

4. Await the resulting UI mount via
   `viewer.awaitRendered(timeoutMs)`.

   * Expected result: a **Biostructure** viewer is mounted —
     `[name="viewer-Biostructure"]` exists; `.msp-viewport canvas`
     is non-empty; default `representation: cartoon`. The
     AutoDock-pose UI (the destination of `importPdbqt` via
     `importPdbqtUI`) is NOT mounted. No "Parsed object is empty"
     toast.

5. Close the Biostructure viewer to leave a clean state for
   Scenario 2 (`viewer.close()` or `grok.shell.closeAll()` on the
   table view).

   * Expected result: viewer closed; no fatal console error during
     teardown.

### Scenario 2 — Double-click `.pdbqt` MUST route to `BiostructureViewer:importPdbqt` (NOT `importPdb`)

Exercises the inverse direction of the disambiguation invariant. A
correct fix MUST keep the `.pdbqt` → `importPdbqt` route intact while
fixing the `.pdb` → `importPdb` route. A naive fix that simply makes
extension matching more aggressive could break the
`.pdbqt` route in the opposite direction; this scenario guards
against that regression.

Steps:

1. From the same **Files** browser (Browse → Files), navigate to a
   `.pdbqt` fixture. Preferred location is under
   **App Data > BiostructureViewer > samples** if a `.pdbqt` is
   shipped there; otherwise fall back to a `.pdbqt` under
   **App Data > Docking** as documented in Setup.

   * Expected result: a `.pdbqt` file is visible in the Files
     browser; no error balloon, no fatal console error.

2. Begin function-call monitoring again (capture the
   fully-qualified name of any `BiostructureViewer:*` function
   invoked over the next action).

   * Expected result: monitoring is active; no captured calls yet.

3. **Double-click the `.pdbqt` file.**

   * Expected result: exactly one
     `BiostructureViewer:importPdbqt` call is captured. NO
     `BiostructureViewer:importPdb` call is captured.
     **Inverse-regression signature**: if
     `BiostructureViewer:importPdb` fires instead, the test
     FAILS with diagnostic "GROK-14442 file-handler search inverse
     regression: .pdbqt routed to importPdb".

4. Await the AutoDock-pose UI mount (the destination of
   `importPdbqt` via `importPdbqtUI`).

   * Expected result: the AutoDock-pose UI is assembled and
     visible. A Mol* (Biostructure) viewer is NOT auto-mounted as
     the file-open destination — the AutoDock-pose UI owns the
     `.pdbqt` opening flow.

5. **Disambiguation cross-check** — at this point both directions
   of the invariant have been asserted independently. Re-state
   the joint invariant explicitly: the test passes iff
   (a) `.pdb` fired `importPdb` AND did NOT fire `importPdbqt`,
   AND
   (b) `.pdbqt` fired `importPdbqt` AND did NOT fire `importPdb`.
   This is the GROK-14442 disambiguation invariant — the
   file-handler search MUST resolve by exact extension equality,
   not by prefix containment.

   * Expected result: both clauses of the joint invariant hold.
     Either-side mis-routing is a regression; only the matched
     handler-fires-per-extension pattern passes.

6. Teardown — close the AutoDock-pose UI and any incidental
   viewers created during the scenario. `grok.shell.closeAll()`
   or per-viewer `viewer.close()` is sufficient; no server-side
   state is created by this scenario.

## Notes

- Deferrals: a fallback resulting-UI-inference assertion path is
  documented in Setup for environments where wrapping
  `grok.functions.call` to spy on the dispatched function name is
  not feasible — a documented fallback, not a coverage gap.

- Fixture availability note: the `.pdbqt` fixture in step 1 of
  Scenario 2 may not ship inside the BiostructureViewer package's
  `samples/` directory. The documented fallback (a `.pdbqt` under
  `System:AppData/Docking/`) keeps the scenario runnable without
  authoring net-new fixtures. If neither location yields a
  `.pdbqt`, stage a minimal AutoDock-pose `.pdbqt` under
  `System:AppData/BiostructureViewer/samples/` as a one-time setup.

- The existing smoke scenario (`biostructure-viewer.md`) opens `.pdb`
  files but never asserts which import handler actually fires, so it
  does not exercise this discrimination invariant. This scenario is
  the dedicated regression guard for that gap.

- This scenario is the `.pdb` vs `.pdbqt` slice of the broader
  extension-routing invariant across `importPdb` / `importPdbqt` /
  `importXYZ` / `importWithNgl`. A future breadth-coverage scenario
  exercising `importXYZ` / `importWithNgl` can extend coverage to the
  remaining slices without touching this scenario.
