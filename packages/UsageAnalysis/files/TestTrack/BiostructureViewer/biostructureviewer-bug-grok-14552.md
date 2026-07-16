---
feature: biostructureviewer
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [GROK-14552]
realizes: []
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

Regression guard for [GROK-14552](https://reddata.atlassian.net/browse/GROK-14552)
("Grid: right-clicking the white space of the row triggers an error"). The
bug surfaced as a `TypeError: Cannot read properties of null (reading
'semType')` when the user right-clicked past the populated cells of a grid
row — the empty whitespace beyond the last column — on any grid holding a
`Molecule3D` column. The platform delivers a `null` cell argument to the
package's grid-cell context-menu hook for that whitespace zone, and the
hook's unguarded access to the cell's semantic type threw. The fix (shipped
in BiostructureViewer 1.19.0) added a null-cell guard so the right-click
silently no-ops on whitespace instead of crashing. The same hook also
registers three context-menu items on populated `Molecule3D` cells — "Copy
Biostructure raw value", "Download Biostructure raw value", and "Show
Biostructure Viewer" — which this scenario also confirms still appear.

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
  - A table whose Molecule3D column is staged programmatically via
    `grok.functions.call('BiostructureViewer:viewBiostructure', {...})`
    or via a `.sdf` whose Molecule3D column is detected by the
    semantic-type detector.
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
     (await `viewer.awaitRendered(timeoutMs)`); the table view shows
     a grid with a `Molecule3D` column. No error balloon, no console
     error.

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
   the grid (a cell rendered by the BiostructureViewer cell
   renderer).

   * Expected result: at least one populated Molecule3D cell is
     identifiable (it renders an inline 3D preview per the cell
     renderer). No error balloon, no console error.

2. Begin console-error capture again (clean buffer).

   * Expected result: capture is armed; the buffer is empty.

3. **Right-click the populated Molecule3D cell.**

   * Expected result: a context menu appears containing the three
     items registered by `addContextMenuForCell`:
     - "Copy Biostructure raw value"
     - "Download Biostructure raw value"
     - "Show Biostructure Viewer"
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

- Deferrals: a fallback no-toast / no-balloon assertion path is
  documented in Setup for environments where `page.on('pageerror')`
  console-error capture is not feasible — a documented fallback, not
  a coverage gap.

- The existing smoke scenario (`biostructure-viewer.md`) does not
  exercise grid-cell right-click flows at all (it never right-clicks
  in empty grid whitespace), so it doesn't cover this bug. This
  scenario is the dedicated regression guard for that gap.

- A separate, broader breadth scenario proposal exists for covering
  all seven grid context-menu / cell-renderer sub-features together;
  this bug-focused scenario overlaps three of those four
  grid-context-menu items but is justified independently as the
  regression guard for GROK-14552.
