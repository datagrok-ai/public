---
feature: biostructureviewer
sub_features_covered:
  - biostructure.file-open
  - biostructure.file-open.importPdb
  - biostructure.file-open.importPdbqt
target_layer: playwright
coverage_type: regression
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

## Overview

Regression guard for [GROK-14442](https://reddata.atlassian.net/browse/GROK-14442)
("Fix file-handler search"). The bug surfaced as a file-handler-search
collision in the platform's extension-to-importer dispatcher: opening a
plain `.pdb` file from the Files browser routed to
`BiostructureViewer:importPdbqt()` instead of
`BiostructureViewer:importPdb()`. The two handlers exist side by side in
`public/packages/BiostructureViewer/src/package.ts` (importPdb at L142
declares extensions `mmcif, cifCore, pdb, gro`; importPdbqt at L177
declares extension `pdbqt`); the search-by-extension logic was selecting
the longer-prefix match (`pdbqt`) when the user opened a file with the
shorter prefix (`pdb`). Because importPdbqt expects AutoDock pose
records (PDBQT parser, AutoDock-pose UI route via `importPdbqtUI`), a
plain `.pdb` routed to it either fails the parse or displays the file
under the wrong pipeline, masking the canonical
`viewBiostructure(content, 'pdb')` happy path entirely.

The bug surfaces across three sub_features tracked in the atlas:

- `biostructure.file-open` — the file-handler family registered by the
  package (`package.ts#L142`). The dispatcher invariant lives at this
  level: each declared extension MUST resolve to exactly one handler,
  and the choice MUST be by exact extension equality, not prefix
  containment.
- `biostructure.file-open.importPdb` — the handler that MUST fire when
  the file extension is `.pdb` (or `.mmcif`, `.cif`, `.gro`). Opens
  content via `viewBiostructure(content, 'pdb')` in the Mol*
  Biostructure viewer.
- `biostructure.file-open.importPdbqt` — the handler that MUST fire
  when (and only when) the file extension is `.pdbqt`. Routes to the
  AutoDock-pose UI via `importPdbqtUI`.

Atlas cross-references:

- `feature-atlas/biostructureviewer.yaml#edge_cases[2]`
  (file-handler extension disambiguation,
  `source_bug: GROK-14442`,
  `derived_from: bug-library:biostructureviewer.yaml#GROK-14442`).
- `feature-atlas/biostructureviewer.yaml#interactions[biostructure-file-open-routing]`
  (cross-feature interaction entry — extension routing across
  importPdb / importPdbqt / importXYZ / importWithNgl — that already
  carries `related_bugs: [GROK-14442]`).
- `feature-atlas/biostructureviewer.yaml#critical_paths[biostructure-file-open-pdb-routes-to-molstar]`
  (p0 smoke path, `derived_from` cites this bug-library entry directly).

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
    report's reproduction). It is also the fixture cited by atlas
    `critical_paths[biostructure-file-open-pdb-routes-to-molstar]`.
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
  `viewer.awaitRendered(timeoutMs)` (see atlas edge_cases[6] —
  async-render await pattern); for the importPdbqt side, await the
  AutoDock-pose UI mount.

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
     Atlas reference: `biostructure.file-open` family
     (`package.ts#L142`).

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
     `BiostructureViewer:importPdbqt` call is captured. Atlas
     reference: `biostructure.file-open.importPdb` →
     `viewBiostructure(content, 'pdb')` (`package.ts#L142`).
     **Regression signature**: if
     `BiostructureViewer:importPdbqt` fires instead, the test
     FAILS with diagnostic "GROK-14442 file-handler search
     regressed: .pdb routed to importPdbqt".

4. Await the resulting UI mount via
   `viewer.awaitRendered(timeoutMs)` (see atlas edge_cases[6] —
   async-render await pattern).

   * Expected result: a **Biostructure** viewer is mounted —
     `[name="viewer-Biostructure"]` exists; `.msp-viewport canvas`
     is non-empty; default `representation: cartoon`. The
     AutoDock-pose UI (the destination of `importPdbqt` via
     `importPdbqtUI`) is NOT mounted. No "Parsed object is empty"
     toast (atlas edge_cases[5] pitfall — the canonical
     `viewBiostructure(content, 'pdb')` entry point supplies the
     name internally, so this guard holds).

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
     `BiostructureViewer:importPdb` call is captured. Atlas
     reference: `biostructure.file-open.importPdbqt` (`package.ts#L177`).
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

- target_layer rationale: **playwright**. The bug manifests at the
  intersection of the platform's file-handler dispatcher (a Datagrok
  shell behaviour) and the Files-browser double-click UI affordance.
  The canonical reproduction action is "double-click a file in the
  Files browser"; the assertion is "which registered function fired".
  A pure-apitest variant could call `grok.functions.call('BiostructureViewer:importPdb', {fileContent})`
  and `grok.functions.call('BiostructureViewer:importPdbqt', {fileContent})`
  directly, but that BYPASSES the file-handler search (the layer
  where the bug lives) — it asserts only that each importer parses
  its own content, not that the dispatcher routes by exact extension.
  The bug is the dispatcher choosing the wrong handler, which only
  manifests through the Files-browser double-click path; that path is
  playwright-only.

- coverage_type rationale: **regression** per the dispatch
  `inputs.bug_brief.coverage_type`. The bug is fixed and this
  scenario is a regression guard against re-emergence of the
  file-handler search collision. Note that the chain-YAML SR
  proposal recommended `coverage_type: edge` for the three residual
  bug-focused scenarios so each authoring would double up against
  F-STRUCT-NEGATIVE-01 (atlas edge_cases[2] carries
  `coverage_type: edge`). The dispatch explicitly resolved to
  `regression` — operator may re-classify this scenario to `edge`
  post-author to satisfy F-STRUCT-NEGATIVE-01 in addition to
  F-BUG-COVERAGE-01.

- Deferrals: a fallback resulting-UI-inference assertion path is
  documented in Setup for environments where wrapping
  `grok.functions.call` to spy on the dispatched function name is
  not feasible. This is NOT a deferral against a missing helper —
  the playwright primitives (UI selectors, DOM queries) are
  first-party — but a documented fallback for environments without
  function-call instrumentation. No sub_features are deferred to
  another layer.

- Fixture availability note: the `.pdbqt` fixture in step 1 of
  Scenario 2 may not ship inside the BiostructureViewer package's
  `samples/` directory. The documented fallback (a `.pdbqt` under
  `System:AppData/Docking/`) keeps the scenario runnable without
  authoring net-new fixtures. If neither location yields a
  `.pdbqt`, the scenario operator may stage a minimal AutoDock-pose
  `.pdbqt` under `System:AppData/BiostructureViewer/samples/` as a
  one-time setup; this is fixture-staging, not test logic, and
  does not affect the assertions.

- Cross-reference with the existing migrated smoke
  `biostructure-viewer.md`: its `bug_match_attempts_skipped[]`
  audit record for GROK-14442 reads "the scenario opens .pdb files
  (Blocks A/B/E) but does not assert which import handler fires,
  so semantic_match_all_steps cannot deterministically extract step
  references that exercise the discrimination invariant". This is
  the dedicated regression guard authored to close that gap
  (F-BUG-COVERAGE-01 branch (ii) anchor via
  `related_bugs: [GROK-14442]`).

- Cross-reference with atlas `interactions[biostructure-file-open-routing]`
  (`feature-atlas/biostructureviewer.yaml#L857-L872`): that
  cross-feature interaction declares the broader extension-routing
  invariant across importPdb / importPdbqt / importXYZ /
  importWithNgl, and already carries `related_bugs: [GROK-14442]`.
  This scenario is the concrete realization of the `.pdb` vs
  `.pdbqt` slice of that broader invariant — the slice where the
  bug actually manifested. A future breadth-coverage scenario
  exercising importXYZ / importWithNgl can extend coverage to the
  remaining slices without touching this scenario.
