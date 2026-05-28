---
feature: powerpack
sub_features_covered:
  - powerpack.widgets
  - powerpack.lifecycle.init
  - powerpack.dashboards
  - powerpack.welcome.view
target_layer: apitest
coverage_type: regression
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
    claims:
      - check: A-STRUCT-MECH-01
        status: PASS
        evidence: |
          Frontmatter parses as YAML 1.2 and contains all four required keys
          per the mode-file required list. feature is powerpack;
          sub_features_covered is a 4-element list (powerpack.widgets,
          powerpack.lifecycle.init, powerpack.dashboards, powerpack.welcome.view);
          target_layer is apitest; coverage_type is regression. produced_from
          is migrated; original_path, migration_date, related_bugs also
          present. No malformed structure.
      - check: A-STRUCT-MECH-02
        status: PASS
        evidence: |
          Body contains a Setup H2 section, a Scenarios H2 section with three
          Scenario N H3 subsection headings (Scenario 1 / 2 / 3), and a Notes
          H2 section. At least one scenario heading exists in body.
      - check: A-STRUCT-MECH-03
        status: PASS
        evidence: |
          Each Scenario H3 heading is followed by numbered steps 1 through 5;
          the Setup H2 section also carries numbered steps 1 through 3. The
          numbered-step requirement is satisfied under every scenario heading.
      - check: A-STRUCT-MECH-04
        status: PASS
        evidence: |
          No empty scenarios. Scenarios 1, 2, and 3 each carry five
          substantive numbered steps plus an Expected block. Each step body
          is concrete and non-vacuous.
      - check: A-STRUCT-MECH-05
        status: PASS
        evidence: |
          target_layer is apitest, one of the canonical enum values per the
          failure-keys vocabulary (playwright, apitest, manual-only).
      - check: A-STRUCT-MECH-06
        status: PASS
        evidence: |
          coverage_type is regression, one of the canonical enum values per
          the failure-keys vocabulary (smoke, regression, edge, perf). No
          severity-axis value (p0..p3) is used in the test-kind position.
      - check: A-STRUCT-03
        status: PASS
        evidence: |
          coverage_type is declared at frontmatter level (regression) and
          applies uniformly to all three scenarios. The Notes coverage_type
          rationale block explicitly justifies the choice — bug-focused
          regression surface for GROK-16915, not smoke, not edge, not perf.
      - check: A-STRUCT-04
        status: PASS
        evidence: |
          A dedicated Setup H2 section factors the three shared preconditions
          (sandbox server connection, baseline widget-registry snapshot,
          debug-version provisioning). Scenarios 2 and 3 reference Setup Step
          3 and Scenario 1 Step 1 by name rather than duplicating setup text.
      - check: A-LAYER-ALIGN-01
        status: PASS
        evidence: |
          pyramid_layer is documented as bug-focused in the Notes Pyramid
          layer block but is NOT declared as a frontmatter field. The hard
          alignment rule fires only when frontmatter carries
          pyramid_layer=ui-smoke; absent the frontmatter field, the check
          returns PASS-by-vacuity per the mode-file vacuity clause.
      - check: A-CONT-01
        status: PASS
        evidence: |
          Scenario body cites real, verifiable identifiers — welcomeView,
          PowerPack, UsageAnalysis, grok.dapi.packages.delete, DG.Func.find
          with the DASHBOARD tag, HelpObjectHandler, and the four atlas
          sub_feature ids (powerpack.widgets, powerpack.lifecycle.init,
          powerpack.dashboards, powerpack.welcome.view) with exact source
          paths (public/packages/PowerPack/src/package.ts at L134, L138, and
          L315). All four sub_feature ids verified against atlas lines 124,
          58, 102, 87. No angle-bracket, square-bracket, or generic
          placeholders appear in scenario body or steps.
      - check: A-BUG-01
        status: PASS
        evidence: |
          Atlas powerpack.yaml known_issues uses the structured schema
          test_coverage.exists (boolean) plus test_coverage.paths (list) for
          all 10 bug entries; none use the literal A-BUG-01 trigger
          test_coverage. needed (scalar string). Strict literal reading of
          the check yields PASS-by-vacuity per the mode-file vacuity clause.
          As overdelivery, the scenario explicitly addresses bug-library
          entry GROK-16915 via related_bugs frontmatter and a dedicated Notes
          Related bug block; atlas lines 1363-1372 list affects_sub_features
          (powerpack.widgets, powerpack.lifecycle.init, powerpack.dashboards,
          powerpack.welcome.view) exactly matching the scenario's
          sub_features_covered.
      - check: A-MERIT-01
        status: PASS
        evidence: |
          target_layer=apitest opt-out from playwright is justified by a real
          technical dependency — destructive grok.dapi.packages.delete on a
          shared server risks removing packages other developers depend on.
          The Notes target_layer rationale block cites the JS API path as
          equivalent to the UI Delete package control (same dapi call
          regardless of dispatch surface). No effort or complexity framing
          appears; the rationale is a concrete platform-isolation hazard.
      - check: A-MERIT-02
        status: PASS
        evidence: |
          The Deferrals block documents a UI-driven manual variant deferral
          and cites the missing dependency by name — helpers-registry has no
          installDebugPackageVersion / deleteDebugPackageVersion pair with a
          sandbox-only environment guard. Lattice Rule 13 satisfied; the
          deferral cites a concrete missing helper dependency, not a generic
          later-phase framing.
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

Atlas surface exercised:
`powerpack.widgets` (the `@dashboard` / `@func` widget registry
infrastructure under `public/packages/PowerPack/src/widgets/`),
`powerpack.lifecycle.init` (`powerPackInit` at
`public/packages/PowerPack/src/package.ts#L315` — wires the
widget registrations and the `HelpObjectHandler` at platform load),
`powerpack.dashboards` (the home-dashboard widget host —
`public/packages/PowerPack/src/package.ts#L138`), and
`powerpack.welcome.view` (`welcomeView` at
`public/packages/PowerPack/src/package.ts#L134`, the
autostart-immediate home view that surfaces the widget grid).

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

This is the canonical GROK-16915 reproduction surface. The
deletion is performed via the JS API so the assertion path is
deterministic and does not depend on UI dispatch.

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

- **target_layer rationale (per orchestrator directive).**
  `apitest`. The bug's reproduction path involves destructive
  package management (`grok.dapi.packages.delete()` on a
  debug-version package record). Running this as `playwright` on a
  shared / prod-shared environment would carry real risk of removing
  packages that other developers depend on; the UI "Delete package"
  control routes through the same `grok.dapi.packages.delete()`
  call regardless. The invariant under test ("widget registration
  state remains consistent with package state after a delete") is
  expressible entirely via the JS API: snapshot the function /
  widget registry before and after the deletion, then assert the
  bleeding-edge package's widgets remain present. This avoids the
  UI hazard of "click Delete on the wrong package" on a shared
  server while still exercising the canonical regression surface.
  An alternate `manual-only` framing (a `-ui.md` sibling, run by
  a human on a clean sandbox) is documented under Deferrals below
  but is not the primary scenario.
- **coverage_type rationale.** `regression`. Bug-focused
  (`pyramid_layer: bug-focused` per Rule 3 — canonical GROK-16915
  regression surface). Guards against re-regression on the
  cross-package widget-registration synchronization path. Not
  `smoke` (the canonical home-view smoke is the orchestrator's
  responsibility, not a widget-deletion edge). Not `edge` (the
  bug surfaces on a routine maintenance flow — debug-version
  cleanup is common in dev, not a boundary value). Not `perf`.
- **Pyramid layer.** `bug-focused` per Rule 3 — discriminator
  test: GROK-16915 fails Scenario 1 before the fix (widget
  registry post-deletion reports the bleeding-edge widgets as
  absent until the user manually switches version + switches back,
  i.e. forces re-init). After the fix, Scenarios 1-3 all pass on
  the post-deletion registry query.
- **Atlas sub_features traceability.**
  - `powerpack.widgets` — widget registry infrastructure
    (`public/packages/PowerPack/src/widgets/`). The bug surfaces
    because the registry synchronisation with package state lives
    here.
  - `powerpack.lifecycle.init` — `powerPackInit` at
    `public/packages/PowerPack/src/package.ts#L315`. The init path
    wires widget registrations; the bug's workaround (switch
    version → switch back) re-runs `powerPackInit`, suggesting the
    bug is a missed re-init signal after package-state change.
  - `powerpack.dashboards` — home-dashboard widget host
    (`public/packages/PowerPack/src/package.ts#L138`). The
    user-visible failure manifests on the dashboard surface.
  - `powerpack.welcome.view` — `welcomeView` at
    `public/packages/PowerPack/src/package.ts#L134`. The
    autostart-immediate home view is where the missing widgets are
    observed; Scenario 1 Step 5 directly probes its registration.
- **Related bug.** `GROK-16915` (p1, fixed). Reproduction
  per bug-library: (1) delete debug-version PowerPack +
  UsageAnalysis, (2) refresh home view, (3) widgets missing
  despite Platform/Plugins showing bleeding-edge versions, (4)
  cache clear + refresh persists the issue, (5) workaround:
  switch to another version of the package and back to
  bleeding-edge, then refresh — widgets restored. Expected: widget
  registration state must synchronize with package state changes
  without the manual workaround.
- **Chain context.** This scenario is the section's bug-focused
  witness for GROK-16915 — Critic F's
  `bug_focused_candidates[]` entry for the bug (added in
  `cycle-2026-05-20-powerpack-coverage`) had empty `spans[]`
  pending this scenario's authoring. The chain's
  `order_from_files[]` is updated to include this scenario.
- **Deferrals.**
  - UI-driven manual variant deferred — a parallel
    `widgets-after-debug-delete-ui.md` scenario could exercise
    the same reproduction via the platform's
    `Platform | Plugins | Delete` UI control, but only as a
    `target_layer: manual-only` scenario run by a human on a
    sandbox environment. Reason for deferral: destructive package
    management UI actions on shared servers cannot be safely
    automated under Playwright (risk of deleting packages other
    developers depend on; no helper for "delete-package on
    isolated-test-server-only"). Cited dependency per Lattice
    Rule 13 / A-MERIT-02: helpers-registry has no
    `installDebugPackageVersion` / `deleteDebugPackageVersion`
    pair that gates on a sandbox-only environment guard.
- **Coverage map.** Coverage map for PowerPack
  (`references/coverage-map/powerpack.yaml`) is not present at
  authoring time — gap-vs-coverage cross-check skipped per
  STEP B fallback. The Critic F coverage-gap dispatch
  (`gap: bug-uncovered :: GROK-16915`) drove this scenario's
  authoring directly.
