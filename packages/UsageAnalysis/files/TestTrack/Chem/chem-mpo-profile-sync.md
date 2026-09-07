---
feature: chem
realizes_atlas:
  - chem.int.mpo-profile-sync
  - chem.cp.mpo-profile-crud
realizes:
  - chem.calculate.mpo-score
priority: p1
target_layer: playwright
pyramid_layer: integration
coverage_type: regression
produced_from: atlas-driven
related_bugs:
  - id: GROK-19624
    status: fixed
source_text_fixes:
  - >-
    Setup Step 1, Scenario 1 Step 3 and Scenario 3 Steps 11 and 14 named MW as
    the property rule, and Setup Step 1 further claimed smiles.csv carries an MW
    column from a prior Calculate > Chemical Properties run. No copy of
    smiles.csv in the repo has an MW column — data/demo/chem/smiles.csv,
    Chem/files/smiles.csv and Chem/files/demo_files/smiles.csv all carry
    HeavyAtomCount and no MW — so the parenthetical was a false claim about the
    fixture. All four places now name HeavyAtomCount, which is the column the
    paired spec has always typed into the profile editor.
realized_as:
  - chem-mpo-profile-sync-spec.ts
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions:
  - anchor: "Scenario 1 Step 1"
    reduction: >-
      The MPO Profiles view is reached with
      grok.functions.call('Chem:mpoProfilesApp'), not through either route the
      step names (Chem > Calculate > MPO Score... > Manage Profiles, or Browse >
      Chem > MPO profiles).
    rationale: >-
      The subject of this scenario is cross-surface sync after create and
      delete, not app entry. Everything from "Create profile" onward — the
      editor, the property row, the name field, the row actions menu, the delete
      confirmation — is real UI actuation. Reaching the view programmatically
      also keeps Scenario 1 independent of an open molecule table, which the
      Chem top-menu route would otherwise require.
    uncovered: >-
      Neither entry point is exercised: the MPO Score dialog's
      [name="button-Manage..."] button, and the Browse tree route into the MPO
      Profiles view. The scoring dialog is already listed as uncovered below,
      and these are the same dialog's other affordance.
    verdict_status: SCOPE_REDUCTION
  - anchor: "Scenario 3 Steps 11-13"
    reduction: >-
      Profile creation and MPO scoring are driven through the JS API
      (grok.dapi.files.writeAsText of the profile JSON, then DG.Func
      mpoScoreByProfile.prepare().call()) instead of through Chem > Calculate >
      MPO Score... and its Manage Profiles > Create profile route.
    rationale: >-
      The profile-editor UI path is actuated end-to-end in Scenario 1, so the
      editor itself is not left untested. Scenario 3 exists for the scoring
      pipeline's product-state signal — the added MPO column and the spread of
      its values — and the direct call exercises that identically. The
      HeavyAtomCount-for-MW swap that used to be recorded here is not a scope
      reduction at all: smiles.csv never had an MW column, so the scenario text
      was simply wrong and is corrected under source_text_fixes, which applies
      to Scenario 1 and the Setup as well as to this anchor.
    uncovered: >-
      The MPO scoring dialog itself — profile selection, molecule-column
      mapping, OK — is not exercised by any spec in this section.
    verdict_status: SCOPE_REDUCTION
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-08-19-chem-new-05
    timestamp: "2026-08-19T23:00:00Z"
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
  e:
    verdict: PASS
    cycle_id: direct-gate-e-2026-08-22-chem-mpo-profile-sync-r2
    timestamp: 2026-08-22T13:56:40Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: direct-gate-b-2026-08-22-chem-mpo-profile-sync-r4
    timestamp: "2026-08-22T12:48:00Z"
    spec_runs:
      - spec: chem-mpo-profile-sync-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 35
        run_mode: headless-cold
        failure_keys: []
expected_results:
  - anchor: "Step 4"
    expectation: >-
      Saving announces itself through the platform's own
      chem-mpo-profile-changed event within a bounded 6 s — not merely "before
      the wait expired" — and the still-open MPO Profiles list, returned to
      through the editor's breadcrumb rather than rebuilt or reloaded, then
      shows the new profile name. The profile is also readable back from disk
      under exactly one file whose internal name is the name that was typed.
  - anchor: "Step 5"
    expectation: >-
      The Browse tree sidebar under Chem shows the newly created profile name
      without a manual refresh — the entry count increases by exactly one
      compared to before the create action. A tree read that failed is rejected
      before presence or count is decided.
  - anchor: "Step 9"
    expectation: >-
      Deleting through the row actions menu announces itself through
      chem-mpo-profile-deleted within a bounded 6 s, and the profile is then
      absent both from disk — confirmed through files.exists, which unlike
      files.list is not stale straight after a delete — and from the rendered
      MPO Profiles list, which must still be showing its other profiles for that
      absence to mean anything. No page reload.
  - anchor: "Step 10"
    expectation: >-
      The deleted profile name is absent from the Browse tree sidebar and the
      entry count returns to the count observed before Step 4 (GROK-19624
      regression guard: the tree must not retain the stale entry until a manual
      refresh). The read assigns nothing — no `expanded`, no Browse-panel
      rebuild — so the entry can only have gone because the product's own event
      drove the tree's refresh. A tree read that failed, or came back empty, is
      rejected before absence or count is decided, at both ends.
  - anchor: "Step 14"
    expectation: >-
      An MPO score column named "MPO <profile name>" is added to the open table,
      holding more than one distinct value, at least one value strictly below 1,
      and every value within [0, 1]. A uniformly zero or uniformly one column
      fails: it means the scoring never read the property.
---

# Chem — MPO profile CRUD and Browse tree sync

## Setup

1. Open a dataset that has a numeric property column suitable for MPO scoring —
   for example, open Demo Files > chem > smiles.csv, which ships with a
   `HeavyAtomCount` column. That is the column the property rules below use; it
   is present in the file itself, so no prior Calculate run is needed.
   Actuation channel: a manual run opens the file through Browse; the paired
   spec opens it programmatically (`grok.dapi.files.readCsv`) because file-open
   actuation is not the subject of this scenario — it is covered by
   `chem-import-export-formats`. Everything else the MPO surfaces own is driven
   through the UI, with the two exceptions recorded under `scope_reductions`:
   how the MPO Profiles app is reached (Scenario 1 Step 1) and how Scenario 3
   creates and applies its scoring profile.
2. Navigate to the Browse panel (left sidebar) and expand Apps > Chem >
   "MPO profiles" so the profile entries are visible. This expansion happens
   once, in Setup. Every later check reads the already-expanded tree without
   expanding anything again — re-expanding rebuilds the entry list from Chem's
   in-memory profile cache and would mask the very staleness Scenario 2 guards.
3. Note the current number of profile entries shown under that node (the
   baseline count) and the number of profile files on disk. Both must be
   non-zero: on a stand with no profiles every later absence check would pass
   while observing nothing. Dev carries 15 stored profiles.

## Scenarios

### Scenario 1: create a profile and verify immediate sync to Browse tree

Steps:
1. In the top menu choose Chem > Calculate > MPO Score... to open the MPO
   dialog, then click "Manage Profiles" (or navigate directly to the MPO Profiles
   app via Browse > Chem > MPO profiles) to reach the MPO Profiles grok-view.
2. In the MPO Profiles view, click the "Create profile" button.
3. In the profile creation form, enter a distinct profile name such as
   `SyncTest-Profile`, add one property rule for `HeavyAtomCount` with a Linear
   desirability ramp, and click Save.
4. Verify: saving fires the platform's `chem-mpo-profile-changed` event within a
   bounded 6 s, and the still-open MPO Profiles list — returned to through the
   editor's breadcrumb, not rebuilt and not reloaded — then shows the new
   profile name. The profile is also on disk, in exactly one file whose internal
   name is the name that was typed.
5. Verify: the Browse tree sidebar under Chem shows the newly created profile
   name without a manual refresh — the entry count increases by exactly one
   compared to the baseline count recorded in Setup Step 3.

Expected:
- The grok-view list gains one row naming `SyncTest-Profile` immediately after
  Save (Step 4).
- The Browse tree entry count under Chem increases by one immediately after Save
  (Step 5).
- Both surfaces are consistent: the profile visible in the list is the same one
  visible in the Browse tree.

### Scenario 2: delete a profile and verify immediate sync to Browse tree (GROK-19624 regression)

Steps:
6. In the MPO Profiles grok-view list, locate the `SyncTest-Profile` row created
   in Scenario 1.
7. Click the actions menu button (the three-dot "⋮" button) on the
   `SyncTest-Profile` row.
8. Choose "Delete" from the popup menu and confirm if a confirmation prompt appears.
9. Verify: the delete fires the platform's `chem-mpo-profile-deleted` event
   within a bounded 6 s, and the profile row is then gone both from disk and
   from the rendered MPO Profiles list — which must still be showing its other
   profiles, or its being empty would look like a successful delete. Confirm
   the disk side with `files.exists`, not by re-listing the directory: a
   listing keeps serving a just-deleted file for a while in the same session,
   which would report a delete that worked as a failure. No reload.
10. Verify: the deleted profile name is absent from the Browse tree sidebar
    immediately after deletion — the Browse tree entry count returns to the
    baseline count recorded in Setup Step 3 (GROK-19624 regression guard: the
    tree must not retain the stale entry until a manual refresh).

Expected:
- The MPO Profiles list no longer shows `SyncTest-Profile` immediately after
  deletion (Step 9).
- The Browse tree sidebar under Chem no longer shows `SyncTest-Profile` without
  a manual refresh — entry count matches the baseline (Step 10).
- Neither surface requires a page reload or a manual tree-refresh to reflect
  the deletion.

### Scenario 3: MPO scoring against an open dataset writes desirability columns

Steps:
11. Create a new minimal MPO profile (or use an existing profile if one is
    present) via Chem > Calculate > MPO Score... > Manage Profiles > Create
    profile. Name it `ScoreTest-Profile`, add one property rule for
    `HeavyAtomCount` with a Linear desirability ramp spanning that column's own
    min–max, and save.
12. Return to the open table from Setup Step 1 (smiles.csv or equivalent).
13. In the top menu choose Chem > Calculate > MPO Score... to open the MPO
    scoring dialog. In the dialog select `ScoreTest-Profile` as the active
    profile and confirm the molecule column is mapped correctly, then click OK.
14. Verify: a column named `MPO ScoreTest-Profile` is added to the open table,
    holding more than one distinct value, at least one value strictly below 1,
    and every value within [0, 1] — the ramp is sloped across the column's own
    range, so a uniformly zero or uniformly one result means the scoring never
    read `HeavyAtomCount`.
15. Delete `ScoreTest-Profile` via the MPO Profiles view to restore the
    baseline state (no cleanup prompt should block the next test run).

Expected:
- An MPO score column is present in the table after Step 13 (Step 14).
- Column values span a range within [0, 1] reflecting real `HeavyAtomCount`
  desirability scores — a uniform-zero or uniform-one result indicates the
  scoring did not apply to the actual property values (Step 14).

## Automation notes

- target_layer rationale: the scenario exercises cross-view state sync between
  three surfaces that only coexist in a live browser session — the MPO Profiles
  grok-view, the Browse tree sidebar, and the open table's column list. The
  gestures under test (Create profile, the row actions menu, Delete + confirm)
  are DOM actuation, and the grok-view list is asserted through its rendered
  rows. The Browse tree is read through the JS model (`mainTree.children`)
  rather than through the DOM, because the tree virtualizes its rows and a DOM
  count would measure what is scrolled into view; that read still needs a
  browser, since `grok.shell.browsePanel` exists only in a running client.
  Playwright is the only layer that gives both at once.
- The Browse-tree read cannot be performing the refresh it asserts is absent —
  closed by construction, not by argument. Expanding happens once in Setup;
  every later read locates nodes and reads `children` and assigns nothing. This
  matters because `mpoProfilesAppTreeBrowser`'s `refresh()` drops all items and
  repopulates from `MpoProfileManager.ensureLoaded()`, the in-memory cache
  (`Chem/src/package.ts:3170-3175`), and the UI delete empties that cache at
  `mpo-profile-manager.ts:69`, one line before firing its event at `:70`/`:72`.
  A read that re-assigned `expanded` could therefore rebuild the list from an
  already-emptied cache and satisfy the guard with the event never having fired.
  With the read-only walk there is no such branch. Verified end-to-end on dev
  2026-08-22: after a real UI create the read-only walk went 15 → 16 and held
  the new name; after a real UI delete it returned to 15 and did not, with the
  changed and deleted events observed firing. Do not re-add an `expanded`
  assignment and do not "simplify" the read into a fresh Browse-panel build —
  either one refreshes the surface under test and makes the guard a tautology.
- The `grok-browser` chem.md note that used to say
  `grok.functions.call('Chem:mpoProfilesApp', {})` does not build the UI (zero
  `.chem-mpo-action-button` after 12 s, `[probe32]`) was refuted by this pair and
  has been corrected in place at `chem.md:1409-1435`: the call does build the UI,
  and the old probe simply never added the returned view to the shell. On dev
  2026-08-22 the call plus `grok.shell.addView` produced two such buttons 104 ms
  later, twice. Do not re-derive the old claim from an older copy of the doc.
- One chem.md note still contradicts this pair and is wrong; see the report
  accompanying this change. `chem.md:1498-1507`
  says a newly created profile does not appear in the Browse tree and calls
  `Browse > Apps > Chem > "MPO profiles"` a wrong path — that node exists,
  carries the 15 stored profiles, and gained the new profile on create. The
  operator observation it records is most likely about the Files subtree
  (`Browse > Files > AppData > Chem > mpo`), which is a different surface and is
  not what this scenario reads.
- GROK-19624 (status: fixed) is the regression surface for Scenario 2. The
  fixed behaviour is that deleting a profile propagates to the Browse tree
  immediately via custom events (`MPO_PROFILE_DELETED_EVENT`) rather than
  requiring a manual refresh. The scenario guards against re-regression.
- The scenario realizes `chem.int.mpo-profile-sync` (interaction: cross-surface
  sync contract between MPO grok-view and Browse tree) and `chem.cp.mpo-profile-crud`
  (p2 critical path: MPO profile CRUD round-trip including scoring). Both are on
  `realizes_atlas`, which is the channel coverage reads; a claim made only in
  this prose would never be credited. No other scenario in the section claims
  the CRUD path, and this one does create, delete and score end to end.
- MPO Profiles app is accessible via Browse > Chem > MPO profiles OR via
  Chem | Calculate | MPO Score... > Manage Profiles. Both entry points reach the
  same `MpoProfilesView`.
- The scoring step (Scenario 3) asserts that the two-surface sync contract is
  not the only tested behavior — the scoring pipeline that consumes a profile
  also has an observable product-state signal (the added column).
