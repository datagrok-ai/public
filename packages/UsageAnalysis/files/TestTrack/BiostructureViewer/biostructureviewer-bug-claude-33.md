---
feature: biostructureviewer
target_layer: playwright
coverage_type: regression
priority: p1
realizes: [CLAUDE-33]
produced_from: atlas-driven
related_bugs:
  - CLAUDE-33
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T17:00:00Z
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
  f:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T14:00:00Z
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T16:45:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-04-biostructureviewer-migrate-02
    timestamp: 2026-06-04T16:30:00Z
    spec_runs:
      - spec: biostructureviewer-bug-claude-33-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 28
        failure_keys: []
realized_as:
  - biostructureviewer-bug-claude-33-spec.ts
---

# BiostructureViewer — Molstar viewer survives unrelated-view close (CLAUDE-33 regression guard)

Regression guard for [CLAUDE-33](https://reddata.atlassian.net/browse/CLAUDE-33)
("BiostructureViewer: Molstar viewer crashes on unrelated view close"). The
bug surfaced as a `TypeError: Cannot read properties of undefined (reading
'children')` when ANY view — not just the view that hosts the Molstar
viewer — is closed while a Molstar (Biostructure) viewer is docked.

The root cause is an over-broad event subscription with an unguarded DOM
access: the package's `onViewRemoved` handler read
`evtView.root.children[0].children[0].classList.contains('msp-plugin')`
before checking whether the closed view was even the one hosting the
viewer. Because it fired for every closed view regardless of identity, it
crashed whenever the DOM of the closed (unrelated) view didn't have that
nested shape — a leaf preview view, a detached help pane, a view whose DOM
was already torn down. The fix null-guards that access and gates it on the
view-ID equality check, so closing an unrelated view is a no-op.

## Setup

- Datagrok session is logged in; the **BiostructureViewer** package is
  installed and registered. The viewer's per-instance lifecycle hook
  in
  `public/packages/BiostructureViewer/src/viewers/molstar-viewer/utils.ts`
  installs the `grok.events.onViewRemoved.subscribe(...)` handler whose
  `evtView.root.children[0].children[0]` access is the bug's defect
  site. After the CLAUDE-33 fix it MUST null-guard the nested DOM
  access (or gate it on the view-ID equality check) so that closing an
  unrelated view raises no error.
- A Biostructure (Molstar) viewer is docked over an open table view.
  Two equivalent fixtures are acceptable:
  - Open `1bdq.pdb` via the Files browser
    (`System:AppData/BiostructureViewer/samples/1bdq.pdb`); the
    package's `importPdb` handler routes it into the Biostructure
    viewer.
  - Add the viewer programmatically via
    `tv.addViewer('Biostructure')` or
    `grok.functions.call('BiostructureViewer:viewBiostructure',
    {content, format, name})`.
  - In either path, `await viewer.awaitRendered(timeoutMs)` before
    closing any sibling view, so the Molstar viewport is fully
    mounted and the `onViewRemoved` subscription is active.
- For the regression assertion ("no `TypeError: Cannot read properties
  of undefined (reading 'children')`"), the load-bearing question is
  "did the `onViewRemoved` handler raise that signature error when an
  unrelated view was closed?". Two assertion paths are equivalent:
  - **Console-error capture.** Attach a Playwright `page.on('pageerror')`
    listener (or read the platform's console buffer) and assert that
    no error matching `TypeError: Cannot read properties of undefined
    \(reading 'children'\)` (or the broader
    `TypeError.*reading 'children'`) was raised between the
    close-unrelated-view action and assertion. This is the canonical
    signature of the bug.
  - **No-toast / no-balloon-error inference.** The platform surfaces
    uncaught errors via balloon notifications; assert that no
    error-balloon containing the substring `children` (or the
    BiostructureViewer source-stack frame
    `molstar-viewer/utils`) appears after the close action.
- Cleanup: no server-side state is created by this scenario (no
  saved project, no shared connection); a `grok.shell.closeAll()`
  in teardown is sufficient.

## Scenarios

### Scenario 1 — Closing an unrelated view while Molstar is docked MUST NOT raise `Cannot read properties of undefined (reading 'children')`

Exercises the exact reproduction path from the bug report. With a
Molstar (Biostructure) viewer docked, the user opens any other,
unrelated view and then closes it. The package's `onViewRemoved`
subscription fires for the unrelated view's close event; the handler
MUST treat it as a no-op.

Steps:

1. Open the **Files** browser via the left sidebar (Browse → Files).
   Navigate to **App Data > BiostructureViewer > samples** and
   double-click `1bdq.pdb`. The Mol* (Biostructure) viewer mounts
   over the resulting table view.

   * Expected result: the Biostructure viewer renders the structure
     (await `viewer.awaitRendered(timeoutMs)`); the table view holds
     at least one `Molecule3D` column; the `.msp-viewport canvas` is
     non-empty. No error balloon, no console error.

2. Confirm the Molstar viewer is active and that its
   `onViewRemoved` subscription is registered. Inspecting the
   subscription list is implementation-internal; the load-bearing
   precondition is "the viewer's `init()` has completed", which is
   guaranteed by step 1's `awaitRendered`.

   * Expected result: the Molstar viewport is docked, rendered, and
     responsive. No console error.

3. Begin console-error capture (or balloon-error capture per the
   alternative assertion path documented in Setup). The capture
   buffer must be clean before the next step.

   * Expected result: capture is armed; the buffer is empty.

4. **Open an unrelated view.** Any view whose root DOM does NOT
   nest a `.msp-plugin`-classed element at
   `root.children[0].children[0]` is a sufficient reproduction
   target. Cheapest paths (in order):
   - Browse → Functions → search any registered function and open
     its parameter editor (creates a function view).
   - Browse → Help → open the platform-help view.
   - File → New → Table → create a small empty table view
     (`grok.shell.addTableView(DG.DataFrame.fromColumns([...]))`).

   * Expected result: the second (unrelated) view opens and becomes
     active; the Molstar viewer remains docked underneath. No
     console error.

5. **Close the unrelated view.** Click the view's close button
   (`×` on the view tab), or invoke `view.close()` programmatically
   on the unrelated view's handle. This fires
   `grok.events.onViewRemoved` with the unrelated view as
   `evtView`. The BiostructureViewer subscription handler
   evaluates
   `evtView.root.children[0].children[0].classList.contains('msp-plugin')`
   against the unrelated view's root DOM. This is the bug's
   reproduction action verbatim.

   * Expected result: NO `TypeError: Cannot read properties of
     undefined (reading 'children')` is raised. NO error matching
     `TypeError.*reading 'children'` appears in the capture buffer.
     NO error balloon containing the BiostructureViewer
     `molstar-viewer/utils` source-stack frame is surfaced. The
     unrelated view tears down cleanly; the Molstar viewer remains
     docked and responsive. **Regression signature**: if
     `TypeError: Cannot read properties of undefined (reading
     'children')` fires from the `onViewRemoved` handler, the test
     FAILS with diagnostic "CLAUDE-33 Molstar onViewRemoved
     unrelated-view-close crash regressed: handler dereferenced
     `evtView.root.children[0].children[0]` without null-guard
     before the view-ID equality check".

6. **Repeat the open / close cycle for two more unrelated views**
   (different shapes — e.g. a function-parameter view AND an empty
   table view) to exercise the failure mode against a range of
   `evtView.root` shapes. Re-assert after each close that the
   capture buffer remains free of the bug's signature error string.

   * Expected result: each close raises no `TypeError`; capture
     buffer stays clean; the Molstar viewer remains docked and
     responsive across all three open/close cycles.

### Scenario 2 — Closing the Molstar view itself MUST still tear down cleanly

Exercises the positive-direction invariant — a correct fix MUST
null-guard the unrelated-view branch WITHOUT regressing the
documented behaviour when the Molstar viewer's OWN host view is
closed. The `onViewRemoved` handler is supposed to perform
cleanup when `evtView.id === view.id` (or when the fallback
`msp-plugin`-class probe legitimately matches); a naive fix that
early-returns the handler for every event would silence the crash
but break the Molstar viewer's own teardown. This scenario guards
against that direction.

Steps:

1. With the table view from Scenario 1 still open (or re-open via
   `1bdq.pdb` per Setup), confirm the Molstar viewer is docked and
   rendered.

   * Expected result: `.msp-viewport canvas` is non-empty; the
     viewer's settings panel is reachable. No error balloon, no
     console error.

2. Begin console-error capture again (clean buffer).

   * Expected result: capture is armed; the buffer is empty.

3. **Close the Molstar host view** by clicking the view-tab close
   button, or invoke `view.close()` programmatically on the view
   that hosts the Molstar viewer. This fires
   `grok.events.onViewRemoved` where
   `evtView.id === view.id` — the legitimate match path that the
   handler is supposed to act on.

   * Expected result: the Molstar viewer tears down cleanly (the
     viewport is removed from the DOM; the viewer instance's
     subscriptions are disposed). NO `TypeError.*reading
     'children'` is raised. Capture buffer remains clean.
     **Inverse-regression signature**: if the Molstar viewer
     fails to tear down (the viewport DOM lingers, or subsequent
     viewer-add operations regress because the prior instance's
     subscriptions were never disposed), the test FAILS with
     diagnostic "CLAUDE-33 fix over-applied: onViewRemoved
     no-ops for the Molstar host view too — legitimate teardown
     branch missing".

4. **Joint invariant cross-check.** Both directions of the
   onViewRemoved-safety invariant have now been asserted
   independently. Re-state explicitly: the test passes iff
   (a) closing ANY unrelated view raises NO
   `Cannot read properties of undefined (reading 'children')`
   error from the BiostructureViewer subscription (Scenario 1),
   AND
   (b) closing the Molstar host view itself still tears the
   viewer down cleanly via the legitimate
   `evtView.id === view.id` branch (Scenario 2, step 3).
   This is the CLAUDE-33 onViewRemoved-safety invariant —
   the handler MUST null-guard the unrelated-view DOM access
   WITHOUT regressing its own teardown path.

   * Expected result: both clauses of the joint invariant hold.

5. Tear down — `grok.shell.closeAll()`. No server-side state to
   clean up.

   * Expected result: all views closed; no fatal console error
     during teardown; capture buffer still clean.

## Notes

- Deferrals: a fallback no-toast / no-balloon assertion path is
  documented in Setup for environments where `page.on('pageerror')`
  console-error capture is not feasible. This is a documented
  fallback for environments without direct console instrumentation,
  not a gap in coverage.

- The existing smoke scenario (`biostructure-viewer.md`) does not
  cover this bug — it opens and interacts with the Molstar viewer
  but never opens and closes an unrelated view while it is docked,
  which is the reproduction path. This scenario is the dedicated
  regression guard for that gap.

- Source citation for the defect site:
  `public/packages/BiostructureViewer/src/viewers/molstar-viewer/utils.ts#L155-L157`
  — the `grok.events.onViewRemoved.subscribe(...)` handler whose
  unguarded `evtView.root.children[0].children[0]` access was the
  reproduction target.
