---
feature: powerpack
target_layer: apitest
coverage_type: regression
priority: p1
realizes: [package-cleanup-keeps-widgets]
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Powerpack/widgets-after-debug-delete.md
migration_date: 2026-05-23
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
related_bugs:
  - GROK-16915
realized_as:
  - widgets-after-debug-delete-api-spec.ts
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: 2026-05-23-powerpack-migrate-02
    timestamp: 2026-05-23T12:00:00Z
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-05-23-powerpack-migrate-02
    timestamp: 2026-05-23T13:00:00Z
    review_round: 1
    failure_keys: []
  e:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T00:00:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T12:48:42Z
    spec_runs:
      - spec: widgets-after-debug-delete-api-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 19
        failure_keys: []
---

# PowerPack — Widget registration after debug-version package delete (GROK-16915 regression)

Bug-focused regression scenario for `GROK-16915`: deleting the
debug version of a package (PowerPack, UsageAnalysis, or any
widget-registering package) MUST NOT leave the home view with
"phantom uninstalled" widgets. After deletion, the home view's
Welcome / Dashboards / Recent Projects / KPI / Spotlight widgets
registered by the remaining bleeding-edge version of the package
must continue to surface — without the manual version-switching
workaround the bug reporter discovered. Widget registration state
must synchronize with package state changes.

## Setup

1. Connect to a sandbox / dev server with `grok.dapi.packages` write
   access (the test will create and delete debug-version package
   records — destructive on prod data; never run against a shared
   environment).
2. Capture the initial widget-registry baseline: list registered
   `@dashboard` functions and `@func` widgets currently surfaced by
   PowerPack + UsageAnalysis on the home view
   (`DG.Func.find({tags: [DG.FUNC_TYPES.DASHBOARD]})` and the
   `@func` widget registry). Record the count and the per-widget
   identifiers (Spotlight, Community, Recent Projects, KPI, Web,
   HTML, About, Learning, etc.). This baseline is the post-deletion
   expected state.
3. Ensure PowerPack (and UsageAnalysis, if available) is installed
   in **both** a bleeding-edge version and a debug-version on the
   target server — the bug requires the two versions to coexist
   pre-deletion. If a debug version is not currently installed, the
   helper `installDebugPackageVersion(pkg)` (or the manual
   equivalent: `grok s packages publish <pkg> --debug` from the
   plugin source dir) provisions one.

## Scenarios

### Scenario 1: Delete debug-version PowerPack — bleeding-edge widgets remain registered

This is the main reproduction path for GROK-16915. The deletion is
performed via the JS API so the check is deterministic and doesn't
depend on UI interactions.

1. **Snapshot pre-deletion widget registry.** Call
   `grok.dapi.packages.list()` and identify the debug-version
   PowerPack record (the entry whose `version` matches the
   provisioned debug version from Setup Step 3). Also enumerate
   `DG.Func.find({tags: [DG.FUNC_TYPES.DASHBOARD], package: 'PowerPack'})`
   and the `@func` widget registry; record the widget identifiers
   and count.
2. **Delete the debug-version PowerPack package.** Invoke
   `grok.dapi.packages.delete(debugPackage)` on the debug-version
   record from Step 1.
3. **Trigger a home-view refresh equivalent.** Re-fetch the platform
   widget registry — call `DG.Func.find(...)` again with the same
   filter as Step 1, restricted to the bleeding-edge PowerPack
   version. The platform's normal "home view refresh" path drives
   re-enumeration of `@dashboard` / `@func` widgets; the apitest
   exercises that path by re-querying the function registry.
4. **Verify widget registry post-deletion.** Assert: every widget
   identifier captured in the Setup Step 2 baseline is still
   present in the registry after Step 3. The remaining bleeding-edge
   PowerPack version must continue to surface all its widgets
   (Spotlight, Community, Recent Projects, KPI, Web, HTML, About,
   Learning) without manual intervention.
5. **Sanity-check the welcome view binding.** Resolve
   `DG.Func.find({name: 'welcomeView', package: 'PowerPack'})` and
   assert exactly one match — the bleeding-edge version's
   autostart-immediate home view (`powerpack.welcome.view`) is
   still registered.

Expected:
- After deleting the debug-version package, all `@dashboard` and
  `@func` widget registrations belonging to the bleeding-edge
  version remain present in the registry.
- No "phantom uninstalled" state — the registry never falsely
  reports the bleeding-edge widgets as absent.
- The `welcomeView` autostart binding resolves to the bleeding-
  edge package exactly once.

### Scenario 2: Edge — delete debug-version UsageAnalysis under same conditions

This scenario covers the bug's cross-package dimension: the
reporter observed the regression after deleting BOTH debug-PowerPack
AND debug-UsageAnalysis. Verifying the invariant holds for
UsageAnalysis independently establishes that the registration-state
synchronization issue is platform-level, not PowerPack-specific.

1. **Re-establish the debug-version of UsageAnalysis.** Provision a
   debug-version UsageAnalysis package if not present (Setup Step 3
   helper, this time for UsageAnalysis).
2. **Snapshot pre-deletion widget registry for UsageAnalysis.**
   Enumerate `DG.Func.find({tags: [DG.FUNC_TYPES.DASHBOARD], package: 'UsageAnalysis'})`
   and any `@func` widget registrations from UsageAnalysis; record
   the count.
3. **Delete the debug-version UsageAnalysis package.** Invoke
   `grok.dapi.packages.delete(debugUsageAnalysisPackage)`.
4. **Re-query the widget registry.** Same pattern as Scenario 1
   Step 3.
5. **Verify widget registry post-deletion.** Assert: every
   UsageAnalysis dashboard / widget identifier captured in Step 2
   is still present in the registry post-deletion. The bleeding-edge
   UsageAnalysis version's widgets must still be registered.

Expected:
- Same invariant as Scenario 1 — deleting the debug-version of
  ANY package that registers home-view widgets must NOT remove
  bleeding-edge widget registrations of the same package.
- The bug's pre-fix workaround (switch package version → switch
  back → widgets reappear) must NOT be required after the
  GROK-16915 fix.

### Scenario 3: Edge — back-to-back delete of both debug packages (combined regression surface)

The bug's original reproduction exercised PowerPack + UsageAnalysis
debug deletions together. Combining both deletions in a single
test sequence catches regressions where the registration-state
synchronization works for one package in isolation but breaks
under concurrent / sequential cross-package deletion.

1. **Re-provision both debug-version packages** (PowerPack +
   UsageAnalysis) per Setup Step 3.
2. **Snapshot the unified widget registry** for both packages
   (union of Scenario 1 Step 1 + Scenario 2 Step 2).
3. **Delete the debug-version PowerPack package, then the
   debug-version UsageAnalysis package** — back-to-back, in the
   same test sequence.
4. **Re-query the widget registry.**
5. **Verify all bleeding-edge widgets remain registered** —
   union assertion: every widget identifier from Step 2 is still
   present.

Expected:
- The bug's original cross-package reproduction surface is
  guarded: deleting both debug versions back-to-back does not
  leave the home view widget-less.

## Notes

- **Related bug.** GROK-16915 (fixed): deleting the debug versions of
  PowerPack and UsageAnalysis caused their bleeding-edge widgets to
  disappear from the home view, even after a cache clear and refresh.
  The only workaround was switching the package to another version and
  back. Widget registration must now stay in sync with package state
  without that workaround.
- **Deferrals.** A UI-driven manual variant (deleting the package via
  **Platform | Plugins | Delete** instead of the JS API) was
  considered but deferred — destructive package-management actions
  can't safely be automated against a shared server, since it risks
  deleting packages other developers depend on.
