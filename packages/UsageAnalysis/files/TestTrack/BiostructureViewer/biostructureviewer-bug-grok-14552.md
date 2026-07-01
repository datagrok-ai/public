---
feature: biostructureviewer
sub_features_covered:
  - biostructure.grid-context-menu.copy-raw
  - biostructure.grid-context-menu.download-raw
  - biostructure.grid-context-menu.show-biostructure-viewer
  - biostructure.cell-renderer.molecule3d
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
related_bugs:
  - GROK-14552
realized_as:
  - biostructureviewer-bug-grok-14552-spec.ts
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
gate_verdicts:
  f:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T14:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T18:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T15:55:00Z
    spec_runs:
      - spec: biostructureviewer-bug-grok-14552-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 87
        failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T19:30:00Z
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
        status: NA
      - check: A-CONT-01
        status: PASS
      - check: A-BUG-01
        status: PASS
      - check: A-MERIT-01
        status: PASS
      - check: A-MERIT-02
        status: PASS
---

# BiostructureViewer — Grid null-cell right-click safety (GROK-14552 regression guard)

## Overview

Regression guard for [GROK-14552](https://reddata.atlassian.net/browse/GROK-14552)
("Grid: right-clicking the white space of the row triggers an error"). The bug
surfaced as a `TypeError: Cannot read properties of null (reading 'semType')`
raised from the BiostructureViewer package's grid-cell context-menu hook
(`addContextMenuForCell` in
`public/packages/BiostructureViewer/src/utils/context-menu.ts`) when the user
right-clicks past the populated cells of a grid row (the empty whitespace
beyond the last column). The platform delivers a `null` cell argument to the
hook for that whitespace zone, and the unguarded `cell.semType` access threw.
The fix (shipped in BiostructureViewer 1.19.0) added a null-cell guard to
the hook so the right-click silently no-ops on whitespace rather than
crashing.

This scenario is the dedicated regression guard authored to close the
F-BUG-COVERAGE-01 gap for this bug (no prior scenario's `related_bugs[]`
listed `GROK-14552`, and the affected sub_features are outside atlas
`manual_only[]`).

The bug surfaces across four sub_features tracked in the atlas:

- `biostructure.cell-renderer.molecule3d` — the registered cell renderer
  for `Molecule3D` semType
  (`public/packages/BiostructureViewer/src/package.g.ts#L14`). Its mere
  presence on a grid column activates the package's grid-cell hooks for
  that grid view; this is the activation surface for the bug.
- `biostructure.grid-context-menu.copy-raw` — "Copy Biostructure raw
  value" item registered for Molecule3D cells
  (`utils/context-menu.ts#L18`). Its `addContextMenuForCell` predicate
  reads `cell.semType` — the exact null-deref site.
- `biostructure.grid-context-menu.download-raw` — "Download Biostructure
  raw value" item (`utils/context-menu.ts#L29`). Same registration code
  path; same null-cell exposure.
- `biostructure.grid-context-menu.show-biostructure-viewer` — "Show
  Biostructure Viewer" item (`utils/context-menu.ts#L41`). Same code
  path; same null-cell exposure.

Atlas cross-references:

- `feature-atlas/biostructureviewer.yaml#edge_cases[0]`
  (null-cell right-click; `source_bug: GROK-14552`;
  `derived_from: bug-library:biostructureviewer.yaml#GROK-14552`).
- `feature-atlas/biostructureviewer.yaml#known_issues[GROK-14552]`
  (`affects_sub_features: [copy-raw, download-raw,
  show-biostructure-viewer]`; `test_coverage.exists: false` — the
  hole this scenario fills).
- Chain YAML
  `scenario-chains/biostructureviewer.yaml#bug_focused_candidates[GROK-14552]`
  proposed_spec verbatim:
  "Open dataframe with Molecule3D column, right-click in grid
  whitespace past populated cells → assert no console error from
  addContextMenuForCell (null-cell safety)."

## Setup

- Datagrok session is logged in; the **BiostructureViewer** package is
  installed and registered. The grid-cell context-menu hook
  (`addContextMenuForCell`) is wired by the package on init
  (`public/packages/BiostructureViewer/src/utils/context-menu.ts`).
  The hook receives a `cell` argument from the platform on every grid
  right-click event and reads `cell.semType` to decide which menu items
  to inject. After the GROK-14552 fix it must treat a `null` cell as a
  no-op (skip menu injection entirely, raise no error).
- A table with at least one `Molecule3D` column is open in the table
  view, so the grid context-menu hook is live for that grid. Two
  equivalent fixtures are acceptable:
  - `1bdq.pdb` opened via the Files browser
    (`System:AppData/BiostructureViewer/samples/1bdq.pdb`); the
    package's `importPdb` handler routes it into a Molecule3D column.
    Atlas reference: `biostructure.file-open.importPdb`
    (`package.ts#L142`).
  - A table whose Molecule3D column is staged programmatically via
    `grok.functions.call('BiostructureViewer:viewBiostructure', {...})`
    or via a `.sdf` whose Molecule3D column is detected by the
    semantic-type detector. Atlas reference:
    `biostructure.api.viewBiostructure` (`package.ts#L130`).
- For the regression assertion ("no error from the hook"), the
  load-bearing question is "did `addContextMenuForCell` raise
  `Cannot read properties of null (reading 'semType')` during the
  right-click?". Two assertion paths are equivalent:
  - **Console-error capture.** Attach a Playwright `page.on('pageerror')`
    listener (or read the platform's console buffer) and assert that
    no error matching `TypeError: Cannot read properties of null` (or
    the broader `TypeError.*semType`) was raised between right-click
    and assertion. This is the canonical signature of the bug.
  - **No-toast / no-balloon-error inference.** The platform surfaces
    uncaught errors via balloon notifications; assert that no
    error-balloon containing the substring `semType` appears after
    the right-click action.
- Cleanup: no server-side state is created by this scenario (no
  saved project, no shared connection); a `grok.shell.closeAll()`
  in teardown is sufficient.

## Scenarios

### Scenario 1 — Right-click grid whitespace (past the last column) MUST NOT raise `Cannot read properties of null (reading 'semType')`

Exercises the exact reproduction path from the bug report. The platform
delivers a `null` cell to `addContextMenuForCell` for the empty
whitespace beyond the populated row; the hook MUST treat it as a no-op.

Steps:

1. Open the **Files** browser via the left sidebar (Browse → Files).
   Navigate to **App Data > BiostructureViewer > samples** and
   double-click `1bdq.pdb`. The Mol* (Biostructure) viewer mounts and
   the table view holds at least one `Molecule3D` column.

   * Expected result: the Biostructure viewer renders the structure
     (atlas edge_cases[6] async-render await applies —
     `viewer.awaitRendered(timeoutMs)`); the table view shows a grid
     with a `Molecule3D` column. No error balloon, no console error.

2. Locate the grid for the open table view. Identify a row whose
   populated cells stop before the grid's right edge — i.e. there is
   visible whitespace to the right of the last column on that row.
   Width-wise, the grid container is wider than the sum of its
   populated columns (true by default for tables imported from
   `.pdb` / `.sdf`).

   * Expected result: a target row is identifiable; the whitespace
     zone past the last column is non-empty (a clickable region in
     the grid container, not in any cell).

3. Begin console-error capture (or balloon-error capture per the
   alternative assertion path documented in Setup). The capture
   buffer must be clean before the next step.

   * Expected result: capture is armed; the buffer is empty.

4. **Right-click in the empty whitespace of the target row, past
   the last populated cell.** This is the bug's reproduction action
   verbatim — the right-click target is the row's whitespace, NOT a
   populated cell.

   * Expected result: NO `TypeError: Cannot read properties of null
     (reading 'semType')` is raised. NO error matching
     `TypeError.*semType` appears in the capture buffer. NO error
     balloon containing `semType` is surfaced. The platform's
     default context menu MAY appear (or no menu at all — both are
     acceptable; the load-bearing assertion is "no crash from the
     BiostructureViewer hook"). **Regression signature**: if the
     `TypeError: Cannot read properties of null (reading 'semType')`
     fires, the test FAILS with diagnostic "GROK-14552 grid null-cell
     right-click crash regressed: addContextMenuForCell did not
     null-guard the cell argument".

5. Dismiss any context menu that may have appeared (Escape, or
   click elsewhere). Assert the capture buffer remains free of
   the bug's signature error string.

   * Expected result: menu dismissed; capture buffer still clean.

### Scenario 2 — Right-click a populated `Molecule3D` cell MUST still inject the package's three context-menu items

Exercises the positive-direction invariant — a correct fix MUST guard
against null cells WITHOUT regressing the documented behaviour on
populated Molecule3D cells. The three menu items registered by
`addContextMenuForCell` ("Copy Biostructure raw value", "Download
Biostructure raw value", "Show Biostructure Viewer") MUST still appear
for a populated Molecule3D cell. A naive fix that early-returns the
hook for every cell would silence the bug but break this invariant; this
scenario guards against that direction.

Steps:

1. With the table view from Scenario 1 still open (or re-open via
   `1bdq.pdb` per Setup), locate a populated `Molecule3D` cell in
   the grid (a cell rendered by the BiostructureViewer cell renderer
   — atlas reference: `biostructure.cell-renderer.molecule3d` at
   `package.g.ts#L14`).

   * Expected result: at least one populated Molecule3D cell is
     identifiable (it renders an inline 3D preview per the cell
     renderer). No error balloon, no console error.

2. Begin console-error capture again (clean buffer).

   * Expected result: capture is armed; the buffer is empty.

3. **Right-click the populated Molecule3D cell.**

   * Expected result: a context menu appears containing the three
     items registered by `addContextMenuForCell`:
     - "Copy Biostructure raw value" — atlas reference
       `biostructure.grid-context-menu.copy-raw`
       (`utils/context-menu.ts#L18`).
     - "Download Biostructure raw value" — atlas reference
       `biostructure.grid-context-menu.download-raw`
       (`utils/context-menu.ts#L29`).
     - "Show Biostructure Viewer" — atlas reference
       `biostructure.grid-context-menu.show-biostructure-viewer`
       (`utils/context-menu.ts#L41`).
     NO `TypeError.*semType` is raised. Capture buffer remains
     clean. **Inverse-regression signature**: if the three menu
     items do NOT appear (the null-cell guard accidentally
     over-fires on populated cells), the test FAILS with
     diagnostic "GROK-14552 fix over-applied: positive-cell menu
     items missing".

4. **Joint invariant cross-check.** Both directions of the
   null-cell-safety invariant have now been asserted independently.
   Re-state explicitly: the test passes iff
   (a) right-click on row whitespace raises NO `semType`-null error
   (Scenario 1), AND
   (b) right-click on a populated Molecule3D cell still injects the
   three BiostructureViewer menu items (Scenario 2, this step).
   This is the GROK-14552 null-cell-safety invariant — the hook
   MUST null-guard the cell argument WITHOUT regressing the
   populated-cell behaviour.

   * Expected result: both clauses of the joint invariant hold.

5. Dismiss the context menu (Escape) and tear down — close the
   Biostructure viewer (`viewer.close()` or `grok.shell.closeAll()`).
   No server-side state to clean up.

   * Expected result: viewer closed; no fatal console error during
     teardown; capture buffer still clean.

## Notes

- target_layer rationale: **playwright**. The bug manifests at the
  intersection of the Datagrok grid's right-click event-dispatch and
  the BiostructureViewer package's `addContextMenuForCell` hook. The
  canonical reproduction action is "right-click on grid whitespace"
  — a pointer event delivered to the grid container, with the cell
  argument the platform computes from the click coordinates. The
  assertion is "no `semType`-null TypeError fires AND the populated
  cell still gets the three menu items". A pure-apitest variant
  could invoke `addContextMenuForCell` directly with a `null` cell
  argument and assert no throw, but that BYPASSES the grid
  right-click event-dispatch (the layer through which the bug
  actually manifests) — the platform's coordinate-to-cell mapping
  for whitespace clicks is the load-bearing piece, and it is
  playwright-only.

- coverage_type rationale: **regression** per the dispatch
  `inputs.bug_brief.coverage_type`. The bug is fixed (shipped in
  BiostructureViewer 1.19.0 per the bug-library) and this scenario
  is a regression guard against re-emergence of the null-cell crash.
  Note that the chain-YAML SR proposal for this bug recommended
  `coverage_type: edge` so each authoring would double up against
  F-STRUCT-NEGATIVE-01 (atlas `edge_cases[0]` carries
  `coverage_type: edge`). The dispatch explicitly resolved to
  `regression` — the operator may re-classify this scenario to
  `edge` post-author to satisfy F-STRUCT-NEGATIVE-01 in addition to
  F-BUG-COVERAGE-01.

- Deferrals: a fallback no-toast / no-balloon assertion path is
  documented in Setup for environments where `page.on('pageerror')`
  console-error capture is not feasible. This is NOT a deferral
  against a missing helper — the playwright primitives (UI
  selectors, DOM queries, error-listener wiring) are first-party —
  but a documented fallback for environments without direct console
  instrumentation. No sub_features are deferred to another layer.

- Cross-reference with the existing migrated smoke
  `biostructure-viewer.md`: its `bug_match_attempts_skipped[]` audit
  record for GROK-14552 reads "bug.affects intersects no scenario
  sub_features_covered — the scenario does not exercise grid-cell
  right-click flows (Molecule3D / PDB_ID grid cells are present in
  Block D and Block G but the scenario does not right-click in empty
  grid whitespace)". This is the dedicated regression guard authored
  to close that gap (F-BUG-COVERAGE-01 branch (ii) anchor via
  `related_bugs: [GROK-14552]`).

- Cross-reference with atlas `edge_cases[0]`
  (`feature-atlas/biostructureviewer.yaml`): that edge case is
  exactly this bug's null-cell scenario, derived from
  `bug-library:biostructureviewer.yaml#GROK-14552`. This scenario is
  its concrete realization — the regression-guard test for the
  documented edge case.

- Cross-reference with the SR breadth proposal in chain YAML
  `gate_f_verdict.scope_reduction_proposal` (Grid context-menu
  scenario proposal): that proposal would cover seven grid
  context-menu / cell-renderer sub_features in a breadth scenario.
  This bug-focused scenario partially overlaps it
  (`biostructure.cell-renderer.molecule3d` plus three of the four
  grid-context-menu items it would propose) and is value-bearing
  under the bug-repro justification independently of the broader
  breadth scenario — the net-new refusal explicitly exempts
  bug-repro scenarios from the breadth-loop overlap check.
