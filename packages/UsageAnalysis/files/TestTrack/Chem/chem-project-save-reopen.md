---
feature: chem
realizes_atlas:
  - chem.int.project-save-reopen-chem-state
realizes: []
priority: p1
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: smoke
produced_from: atlas-driven
related_bugs:
  - id: GROK-17595
    status: regression-risk
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities:
  - "scaffold-tree-for-nx-testing.tree is absent from DemoFiles/chem and
    AppData/Chem (grok.dapi.files.exists=false both); per scenario Notes the
    spec substitutes auto-generation, driven from the viewer's magic-wand
    Generate control, for the custom-tree upload in Scenario 1 Step 3."
  - "RETRACTED 2026-08-21 (operator): the claim that Scaffold Tree never
    populates nodes on dev was wrong, and it was wrong because of the dataset.
    Generation is bounded to small tables — on the 1000-row smiles.csv the
    Generate action is disabled and the tree stays empty, which reads as a
    silent no-op. On a 100-row table it generates: chem-filters-spec.ts renders
    and clicks scaffold nodes on spgi-100 and is Gate B PASS 3/3. This scenario
    now runs on spgi-100 and asserts node presence rather than waiving it."
  - "Whether the datasync leg of GROK-17595 can be reached from a different save
    route is an OPEN OPERATOR QUESTION. Source evidence (see scope_reductions
    'DATASYNC NOT EXERCISED') says the toggle is presented only for tables that
    carry a creation script, i.e. tables opened through the OpenFile function
    rather than through grok.dapi.files.readCsv. Changing the Setup actuation
    channel is not a decision this spec took on its own."
scope_reductions:
  - "WITHDRAWN 2026-08-21: this reduction rested on the retracted claim that
    generation never produces nodes. It does, on a small enough table. The spec
    now checks the first generated scaffold node — which is what turns scaffold
    filtering ON — and the round-trip is asserted against a table that is
    actually filtered."
  - "CUSTOM .tree BLOB NOT EXERCISED: scaffold-tree-for-nx-testing.tree does not
    exist on the server, so the spec builds the tree in place through the
    viewer's own magic-wand Generate control instead of uploading the fixture.
    The generated tree is exercised rather than waived — Scenario 1 Step 3
    checks the first node and asserts the row count it narrows, on both sides of
    the save — but what the saved project carries is a locally generated tree,
    not an uploaded blob. The specific surface GROK-17595 is about, a custom
    .tree file deserializing out of a saved project, is therefore NOT covered.
    Re-strengthen when the fixture exists."
  - "NON-UI ACTUATION CHANNELS: three gestures this scenario's prose describes
    as UI actions are driven through the JS API instead, and each is listed here
    rather than covered by a general claim. (1) REOPEN: the saved project is
    reopened through `grok.dapi.projects.filter(...).first()` then
    `project.open()` (spec:329-332), not by double-clicking its tile in Browse >
    Projects — a deserialization defect that only surfaces on the Browse open
    path is invisible here. (2) CLOSE ALL: both scenarios say `Select File >
    Close All`, but the spec calls `grok.shell.closeAll()` (spec:326 in the
    reopen path, :86 in Setup, :605 in cleanup). The POST-state is asserted
    (`expect(reopened.tvAfterClose).toBe(0)`, spec:505) but the menu path itself
    is never exercised, so a File > Close All that fails to dismiss views is
    invisible here. (3) TABLE OPEN: Setup says `via File > Browse`, but the spec
    calls `grok.dapi.files.readCsv` (spec:90); that channel also has a datasync
    consequence — see the 'DATASYNC NOT EXERCISED' entry. On a `pyramid_layer:
    ui-smoke` scenario all three gaps are load-bearing. Re-strengthen via the
    Browse-tile open route, the File menu for Close All, and File > Browse for
    the table open."
  - "UI CONTROL SUBSTITUTIONS: gestures this scenario's prose describes as one
    UI control but the spec drives through a different UI control. These do not
    belong in the NON-UI ACTUATION CHANNELS entry above — the channel is still
    the UI — but they are divergences between the prose and the driven control
    and are enumerated here for the same reason. (1) SAVE PROJECT: Scenario 1
    Step 4 and Scenario 2 Step 3 both say `Select File > Save Project`, but the
    spec clicks the table-view toolbar Save button (`[name=\"button-Save\"]`,
    spec:272) and drives the `[name=\"dialog-Save-project\"]` dialog it opens.
    Both routes end in the same dialog, so the save itself is exercised, but the
    File menu is not: a File > Save Project item that is missing, disabled, or
    wired to the wrong command is invisible here, with nothing in the spec
    recording that. Re-strengthen by opening the dialog from the File menu."
  - "DATASYNC NOT EXERCISED: the GROK-17595 triple is exercised as a PAIR
    (active Substructure Search + Scaffold Tree), without datasync. The
    Data-sync switch in the Save-project dialog is rendered display:none for
    this table
    (core/client/xamgle/lib/src/features/project_entity_move.dart:513-514 hides
    the switch when the TableInfo carries no Tags.CreationScript; the creation
    script is stamped by DataHistory.logDataFrameCreationScript only for tables
    produced by a FuncCall such as OpenFile, and this scenario opens the table
    through grok.dapi.files.readCsv). The spec reports the observed datasync
    state per save instead of silently skipping it. Re-strengthen only via an
    operator decision on the Setup actuation channel."
gate_verdicts:
  b:
    verdict: PASS
    cycle_id: direct-gate-b-2026-08-22-chem-project-save-reopen-r6
    timestamp: "2026-08-22T14:53:10Z"
    spec_runs:
      - spec: chem-project-save-reopen-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 97
        run_mode: headless-cold
        failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-08-19-chem-new-05
    timestamp: "2026-08-19T12:00:00Z"
    failure_keys: []
    review_round: 2
    claims:
      - check_id: A-STRUCT-MECH-01
        status: PASS
      - check_id: A-STRUCT-MECH-02
        status: PASS
      - check_id: A-STRUCT-MECH-03
        status: PASS
      - check_id: A-STRUCT-MECH-04
        status: PASS
      - check_id: A-STRUCT-MECH-05
        status: PASS
      - check_id: A-STRUCT-MECH-06
        status: PASS
      - check_id: A-STRUCT-03
        status: PASS
      - check_id: A-STRUCT-04
        status: PASS
      - check_id: A-LAYER-ALIGN-01
        status: PASS
      - check_id: A-CONT-01
        status: PASS
      - check_id: A-BUG-01
        status: PASS
      - check_id: A-MERIT-01
        status: PASS
      - check_id: A-MERIT-02
        status: PASS
  e:
    verdict: PASS
    cycle_id: direct-gate-e-2026-08-22-chem-project-save-reopen-r6
    timestamp: 2026-08-22T14:59:41Z
    failure_keys: []
expected_results:
  - anchor: "Step 5"
    expectation: >-
      After reopening the project the Scaffold Tree viewer is restored in the
      layout — present both among the table view's viewers and as a DOM element.
      The scaffold leg's CRITERION is restored too, and it is separately
      observable: checking a scaffold node ANDs a scaffold BitSet into the
      dataframe filter, so the row count drops below what the Substructure
      Search passes on its own. Scenario 1 Step 3 records both numbers before
      the save — the substructure-only count and the narrower count with the
      node checked — and Step 5 requires the reopened table to settle on the
      narrower one. A restored project that came back with the substructure
      filter alone lands on the wider count and fails. That equality is
      one-sided on its own — see the next entry for the release step that makes
      the substructure leg observable too. The tree is generated in place rather
      than uploaded from the custom fixture — see the scope_reductions entry
      "CUSTOM .tree BLOB NOT EXERCISED".
  - anchor: "Step 5"
    expectation: >-
      After reopening the project the Substructure Search filter is restored,
      and it is observable independently of the scaffold criterion. The equality
      alone cannot show this: the root scaffold contains a benzene ring, so
      every row it matches also matches the benzene query and the scaffold-only
      row set is a SUBSET of the substructure-only one — a reopen that restored
      the scaffold criterion and dropped the Substructure Search lands on the
      same number. So Step 5 also RELEASES the scaffold criterion after the
      count has settled: unticking the restored node must widen the row set back
      to exactly the substructure-only count noted in Step 3. Discrimination,
      both directions. With the node still ticked: the substructure-only count
      means the SCAFFOLD criterion was lost, an empty row set means a BitSet
      deserialized empty, the full row count means nothing deserialized at all.
      With the node unticked: the full row count means the SUBSTRUCTURE
      criterion was lost, and a count still at the narrower ticked value means
      unticking released nothing — so the equality before it was not measuring a
      live scaffold filter either.
  - anchor: "Step 5"
    expectation: >-
      No console errors of the form "error: null" or "null" appear after the
      project is reopened, and no error balloon is raised during the reopen (the
      two channels a GROK-17595 deserialization error can reach — see the "Error
      channels" note).
  - anchor: "Scenario 2 (step 3)"
    expectation: >-
      Reopening the project that carries only the Scaffold Tree (no Substructure
      Search active) completes without console errors and the Scaffold Tree
      viewer is restored in the layout, with no error balloon raised. This
      baseline deliberately leaves the tree unchecked — no scaffold criterion is
      applied, so there is no filtered row set to compare and restoration is
      judged by viewer presence alone, unlike Scenario 1 where the checked node
      is what makes the round-trip observable. Datasync is not enabled: see the
      scope_reductions entry "DATASYNC NOT EXERCISED".
realized_as:
  - chem-project-save-reopen-spec.ts
---

# Chem — project save and reopen with Chem state

## Setup

1. Open `AppData/Chem/tests/spgi-100.csv` via File > Browse. A 100-row table is required:
   Scaffold Tree generation is bounded to small tables and does nothing on a 1000-row one.
   Actuation channel: the automated spec opens the table programmatically
   (`grok.dapi.files.readCsv`) rather than through the Browse UI. This has a
   consequence specific to this scenario — see the datasync note under Notes.
2. Confirm the molecule column is recognised (cells render structures, not raw text).

## Scenarios

### Scenario 1: Chem state round-trip (Substructure Search + Scaffold Tree; datasync not exercised)

GROK-17595 describes a fragile triple: active Substructure Search + Scaffold Tree with a
custom uploaded tree blob + datasync ON. This scenario exercises the pair — Substructure
Search plus Scaffold Tree — with an auto-generated tree and without datasync. The two
missing legs are declared in `scope_reductions`.

Steps:
1. With `spgi-100.csv` open, select top menu Chem > Search > Substructure Search.
   In the sketcher that opens, draw or type the SMILES `c1ccccc1` (benzene ring) and
   confirm the query. Verify the table filters: the filtered row count displayed in
   the status bar is greater than zero and less than the total row count of the table.
2. Select top menu Chem > Analyze > Scaffold Tree. Verify the Scaffold Tree viewer
   panel appears attached to the table view.
3. Generate the tree with the magic-wand control in the Scaffold Tree viewer toolbar.
   The control is only enabled below 500 distinct molecules, so confirm the table is
   inside that bound first; above it the click is swallowed and the tree stays empty.
   Verify the generation completes without raising an error and that scaffold nodes
   appear. Then tick the checkbox on the first scaffold node: generation alone filters
   nothing, and checking a node is what turns scaffold filtering on. Verify the effect
   on the table, not the checkbox: the row count must drop below what the Substructure
   Search alone was passing. Note both counts — they are the reference the reopen is
   judged against. The custom `.tree` fixture the bug report uses does not exist on
   the server, so the tree is generated rather than uploaded — declared in
   `scope_reductions`.
4. Select File > Save Project. In the save dialog set the project name to
   `chem-state-roundtrip-test` and confirm save. The Data-sync switch is not presented
   for this table, so datasync stays off; record which state was observed. Verify the
   project is queryable on the server after the dialog closes.
5. Select File > Close All to dismiss all open views. Then reopen the saved project and
   wait for the views to finish loading. Once the row count has stopped moving, untick
   the checkbox on the restored scaffold node and watch the row count again: releasing
   the scaffold criterion must leave the Substructure Search filtering on its own.

Expected:
- The Scaffold Tree viewer is present in the restored layout (viewer list and DOM).
- The restored filters reproduce the exact row set that was saved: the row count settles
  on the number noted in Step 3 with the scaffold node checked — not on the wider
  substructure-only number, not on zero, not on the full table. Filtering is not applied
  at the moment the layout returns — the criteria come back first and the row set is
  recomputed after — so this is read once the count settles and holds, not the moment it
  first touches the number.
- With the restored node unticked the row count widens back to exactly the
  substructure-only number noted in Step 3. That is what makes the two criteria
  separately observable: the root scaffold carries a benzene ring, so the rows it matches
  are a subset of the rows the benzene query matches, and the settled count above would
  be reproduced by a reopen that restored the scaffold criterion and lost the Substructure
  Search. Untick and the full table means the substructure criterion was lost; a count
  still at the narrower value means unticking released nothing.
- No console errors of the form `"error: null"` or containing the string `"null"` fire
  during or after reopen, and no error balloon is raised (GROK-17595 deserialization
  surface; see the "Error channels" note for why both are watched).

### Scenario 2: Scaffold Tree alone — minimal round-trip guard

This narrower path isolates whether the Scaffold Tree state itself (without a concurrent
Substructure Search) round-trips cleanly, providing a baseline that lets Scenario 1
failures be attributed to the interaction rather than the Scaffold Tree alone.

Steps:
1. Open `AppData/Chem/tests/spgi-100.csv`. Select top menu Chem > Analyze > Scaffold Tree.
2. Generate the tree with the magic-wand control in the Scaffold Tree viewer toolbar.
   Verify the generation completes without raising an error, that scaffold nodes
   appear, and that the viewer is attached to the table view. No
   scaffold node is checked here: this baseline saves an unfiltered table on purpose,
   so that a Scenario 1 failure can be attributed to the Substructure Search
   interaction rather than to the Scaffold Tree on its own.
3. Select File > Save Project with name `chem-scaffold-only-test`. Datasync is not
   presented and stays off. Confirm save. Select File > Close All, then reopen the
   saved project.

Expected:
- The Scaffold Tree viewer is present in the restored layout (viewer list and DOM).
- No console errors of the form `"error: null"` appear during or after reopen, and no
  error balloon is raised.

## Notes

- **Correction 2026-08-22 — the scaffold criterion IS observable, and the earlier
  claim that it was not was wrong.** The first `expected_results` entry used to say
  that "the scaffold leg's round-trip is NOT separately observable in this
  configuration", on the reasoning that after reopen the only signal is the row set
  and any row count below the full table is already explained by the Substructure
  Search filter. The section's own UI reference refutes that:
  `.claude/skills/grok-browser/references/chem.md:441-443` records a measured,
  falsifiable signal — "Checking a node's checkbox filters the table ... 50 → 47 rows
  after one click, read via `viewer.dataFrame.filter.trueCount`". Checking a scaffold
  node ANDs a scaffold BitSet into the filter, so the row set with the node checked is
  strictly narrower than the substructure-only row set, and both numbers exist on both
  sides of the save. The round-trip is therefore asserted as an exact equality against
  the pre-save number, with the substructure-only number named in the failure message
  as the value a lost scaffold criterion would land on. This is a factual correction,
  not a scope reduction: nothing was narrowed, a claim that had been given up was
  recovered.
- **Error channels.** A GROK-17595-style `"error: null"` can surface through two
  independent channels, and only one of them was watched. `console.error` is real —
  the reopen path calls `log.error(x, s)`
  (`core/client/xamgle/lib/src/shell/shell_project.dart:412-414`) and the client
  logger writes every `LogLevel.ERROR` to `window.console.error`
  (`core/client/xamgle/lib/src/features/logging/logger.dart:35`). But it is not the
  only one: `Balloon.error(...)` logs to the console only when called with
  `logError: true` (`core/client/d4/lib/src/widgets/balloon/balloon.dart:131-134`),
  so a user-visible error balloon can sit on screen with a clean console. The spec
  now watches `.d4-balloon.error` with a MutationObserver installed before the reopen
  (collecting as balloons appear, rather than snapshotting at the end, so an early
  error that has since faded is still caught) and, as a positive control on the
  console arm, round-trips a sentinel through `console.error` and requires the trap to
  have captured it — without that, an empty error list is equally consistent with "no
  errors" and "the trap was never installed".

- target_layer rationale: the scenario drives multi-step UI state across three
  Chem subsystems (Substructure Search modal, Scaffold Tree viewer, project
  save/reopen dialog), requires a real project persistence round-trip (save →
  close all → reopen), and the critical assertion (no `"error: null"` on reopen)
  can only be observed in a live browser session. playwright is the only layer
  capable of this.
- realizes_atlas: `chem.int.project-save-reopen-chem-state` is the interaction id
  that addresses the cross-package contract: Chem + Projects + datasync.
- GROK-17595 repro surface: the fragile triple from the bug library is
  "substructure search + scaffold tree (custom uploaded tree blob) + datasync =
  fragile interaction during serialization". This scenario reproduces two of the
  three legs. What is NOT reproduced, and why:
  - **The custom tree blob.** `scaffold-tree-for-nx-testing.tree` does not exist on
    the server, so the spec builds the tree in place through the viewer's own
    magic-wand Generate control and then checks its first node, which is what turns
    scaffold filtering on. The
    tree driven into the save is therefore a locally generated one rather than an
    uploaded blob. Whether the scaffold criterion survives the round-trip IS
    asserted — through the row count it narrows, see the correction note above.
    What stays uncovered is deserializing a custom `.tree` file that was uploaded
    into the project, which is the specific surface GROK-17595 is about.
  - **Datasync.** The Data-sync switch in the Save-project dialog is only rendered
    for a table that carries a creation script, and a table opened through
    `grok.dapi.files.readCsv` has none, so the switch is `display:none` and the
    save is a plain, non-synced save. The spec records the observed datasync state
    at each save rather than skipping the question silently.
- Actuation channels: the automated spec reopens the project through the API rather
  than double-clicking a tile in Browse > Projects, calls `grok.shell.closeAll()`
  rather than File > Close All, and opens the table with `grok.dapi.files.readCsv`
  rather than File > Browse — all three declared in `scope_reductions` as
  "NON-UI ACTUATION CHANNELS". It also opens the Save-project dialog from the
  table-view toolbar Save button rather than from File > Save Project; that one
  stays inside the UI, so it is declared separately as "UI CONTROL
  SUBSTITUTIONS".
- Cleanup: after both scenarios the test projects `chem-state-roundtrip-test` and
  `chem-scaffold-only-test` are deleted, in a `finally` block so that server-side
  state is removed even when an assertion fails.
- See: public/help/datagrok/solutions/domains/chem/chem.md (project save/reopen
  with Chem features).
- See: public/help/visualize/viewers/scaffold-tree.md (Scaffold Tree upload /
  generate controls).
