---
feature: projects
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [derive-then-save-inside-project, share_with_recipient_open, rename_project]
realizes: [views.projects]
realized_as:
  - projects-lifecycle-derived-spec.ts
pyramid_layer: proactive-lifecycle
ui_coverage_responsibility: []
ui_coverage_delegated_to: projects-ui-smoke.md
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Projects/projects-lifecycle-derived.md
migration_date: 2026-05-04
related_bugs:
  - GROK-19103
---

# Projects — Derived-source lifecycle (Pivot / Aggregate / Join)

Covers the lifecycle of a project containing derived tables — Pivot
Table, Aggregate Rows, and Join Tables results built on top of a
file-share source. Verifies save/reopen, sharing with a second user,
and project rename all correctly preserve the derivation chain, and
reinforces GROK-19103: a Join result must land in the active
workspace, not get silently saved as a separate, broken project.

UI coverage delegated to `projects-ui-smoke.md`. The Pivot / Aggregate
/ Join dialogs themselves are exercised by `uploading.md`'s
source-matrix scenarios, not here.

## Setup

1. Authenticate as test user.
2. Project name: `lifecycle-derived-${Date.now()}`.
3. Recipient placeholder: `<RECIPIENT_USERNAME_TBD>`.
4. Helper 3 dependency: `helpers.playwright.session.logoutAndLoginAs`
   (NOT YET REGISTERED). Same deferral pattern as
   `projects-lifecycle-files.md`.
5. Cleanup: delete project; revoke permissions.

## Scenarios

### Main flow — derived-tables lifecycle

1. **Open parent table from File share.** Open
   `System:DemoFiles/demog.csv`. Verify table loaded.
2. **Build derived tables in the workspace.**
   - Pivot: configure rows / columns / values, click **Add**.
     Resulting Pivot table appears as a new tab.
   - Aggregate Rows: group + aggregate on the demog table, click
     **Add**. Resulting aggregate appears as a new tab.
   - Join Tables: join the Pivot result with the original demog on
     a key column. Resulting Join appears as a new tab.
   - **GROK-19103 invariant assertion:** the Join result lands in
     the **active workspace** (visible in `grok.shell.tables`),
     NOT as a stray separate project. Verify
     `grok.shell.tables.length` increased by 1 (Join), not by 2
     (Join + stray project).
3. **Save project with Data Sync ON.** Trigger Save Project (ribbon
   button). Project name from Setup, Data Sync **ON** for parent
   AND each derived table, click **OK**. Cancel auto-share.
4. **Share project with second user (View-and-Use + Full).** Use
   `grok.dapi.permissions.grant`. Verify Sharing tab lists the
   recipient.
   - **Original-user assertion:** project reopens; all 4 tables
     (parent + 3 derived) load; derivation links preserved.
   - **Recipient-side assertion (Helper 3 — deferred):** second
     user opens; same 4 tables; same derivation links.
5. **Rename external dependency — N/A for derived.** Derived
   tables inherit their parent (here: `demog.csv`); the parent
   File source has no rename surface (path is fixed). Per chain
   rev 3 `proactive_lifecycle_specs[5].dep_lifecycle_ops_covered:
   [share_with_recipient_open]`, `rename_external_dep` is NOT
   in scope.
6. **Rename project itself.** Via JS API. Verify rename persists.
   - Original-user assertion: project still opens under new name;
     all derivations preserved.
   - Recipient-side assertion (Helper 3 — deferred).
7. **Cleanup.** Delete project; revoke permissions.

### Expected results

- Derived tables persist correctly across save → reopen.
- Join result does NOT leak into a stray separate project (GROK-19103
  invariant).
- Project rename does not break the derivation chain.
- Share + recipient-open works for derived-source projects.

## Notes

- **Reinforces GROK-19103.** Step 2 explicitly asserts the
  GROK-19103 invariant (a Join result lands in the active project,
  not a stray one). This is a reinforcement — the primary bug-focused
  spec for that bug is `complex-derived-tables-spec.ts`.
- **No external rename for derived sources.** Derived tables inherit
  their parent's source (here, a fixed file path with no rename
  surface), so external-dependency rename isn't applicable to this
  entry.
- **UI coverage delegated.** The Pivot / Aggregate / Join UI surfaces
  are owned by `uploading.md`'s source-matrix cases; this scenario
  uses the JS API and assumes that matrix coverage.
- **Deferred.** Recipient-side assertions are blocked on a
  not-yet-registered login-as-another-user test helper, same as
  `projects-lifecycle-files.md`.
- **Self-cleaning.** Step 7 deletes the project.
