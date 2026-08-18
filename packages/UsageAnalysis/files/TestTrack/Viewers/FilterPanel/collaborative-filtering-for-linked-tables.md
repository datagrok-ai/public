---
feature: filters
realizes_atlas: [ filters.cp.linked-tables-collaborative ]
realizes: [ viewers.filters ]
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-19137
    status: fixed
realized_as:
  - collaborative-filtering-for-linked-tables-spec.ts
precondition_guards:
  - description: >-
      Three tables are open as named table views — spgi-100, SPGI-linked1 and
      SPGI-linked2 — each with its expected row count, and the two links have
      been requested (df3↔df2 FILTER_TO_FILTER on four key columns, df1↔df2
      SELECTION_TO_FILTER on Id/Concept Id). The link registry itself is not
      exposed in the JS API, so link correctness is proved by its effect at
      Steps 2, 6 and 8 rather than by inspecting the registry.
    check: >-
      grok.shell.tableViews contains views named spgi-100, SPGI-linked1 and
      SPGI-linked2, with rowCount 100 / 3624 / 224 respectively
scope_reductions:
  - id: SR-01
    check: E-TRACE-02
    rationale: |
      "Switch to the <name> view" is written as a click on the view's tab. The
      spec substitutes the JS API instead: it assigns grok.shell.v to the table
      view whose dataFrame carries that name. View switching is not the
      behaviour under test here — the subject is how a filter or a selection
      propagates across a link — and the tab strip carries no stable per-view
      handle to address by name. The substitution is made non-silent: the
      helper asserts that a view with that name was actually found (a
      no-op switch used to leave a later gesture on whatever view happened to
      be active), and then polls until grok.shell.tv reports that name, so a
      switch that did not take is a failure rather than a wrong-frame pass.
      Every step that then acts on the active frame also asserts the frame's
      own identity at the point of the gesture.
    verdict_status: SCOPE_REDUCTION
  - id: SR-02
    check: E-TRACE-02
    rationale: |
      Steps 5 and 8 name the CARD gesture: "in the link column 3 filter card,
      select the category v ii", and the same shape for PAMPA Classification on
      SPGI-linked1. The spec does not click those category rows — it sets the
      criterion through the filter group's updateOrAdd. A categorical card's
      category list is a canvas grid: a plain click on the checkbox column of an
      all-selected card trips the product's own inversion rule and ends up
      selecting only the clicked category, and Ctrl/Shift-click suppresses the
      inversion but also moves the row SELECTION — which in this scenario is the
      very channel one of the two links carries, so a selection-moving gesture
      would contaminate the propagation under test. The subject here is
      propagation ACROSS a link, not the click that sets a criterion.
      The substitution is adequate for a propagation claim here, and that is
      established rather than assumed. A link subscribes to
      `dataFrame1.onFilterChanged` / `onSelectionChanged` and contributes its
      contribution from `dataFrame2.onRowsFiltering`
      (`ddt/lib/src/data_frame/tools/link.dart:86-106,159-186`); nothing in the
      link machinery listens to `FILTER_CRITERIA_CHANGED`, the one notification
      the card click raises and `updateOrAdd` does not. Both paths therefore
      reach the link through the same door, and both were measured on dev
      2026-08-18: driving Step 5's criterion by clicking the "v ii" category-name
      row on the SPGI-linked2 card left SPGI-linked2 at 148 rows and moved
      SPGI-linked1 from 9 to 5, and driving Step 8's criterion by clicking the
      "inconclusive" row on SPGI-linked1 left it at 2 rows while SPGI-linked2
      held 148 and spgi-100 held 100 filtered / 5 selected at every sample of a
      4 s window — the same numbers, including the same bystander zeroes, that
      the API path produces.
      What compensates: the Filter Panel is still opened through the ribbon on
      the frame the criterion lands on, the criterion is read back off the
      surviving rows on both steps (every surviving row must carry the chosen
      category), the frame's own identity is asserted at the point of the
      gesture, and every count the steps assert is derived from the raw columns.
      What is NOT covered is the card click itself, on either step.
    verdict_status: SCOPE_REDUCTION
expected_results:
  - anchor: Step 1
    expectation: >-
      Exactly 5 rows are selected in spgi-100 (df.selection.trueCount == 5),
      the frame the selection landed on IS spgi-100 (its name and its 100 rows
      are read at the point of the gesture), and spgi-100's OWN filtered row
      count is unchanged at its full 100 rows — selecting is not filtering, so
      the master table must not narrow itself.
  - anchor: Step 2
    expectation: >-
      After selecting 5 rows in spgi-100 under a SELECTION_TO_FILTER link, the
      SPGI-linked1 view shows exactly the DERIVED key-matched row count: the
      number of SPGI-linked1 rows whose Concept Id is among the first five
      spgi-100 Id values, computed from the raw columns rather than pinned as a
      literal. In the same breath the two frames that must NOT have moved are
      read twice, a settle apart: SPGI-linked2 still shows all 224 rows
      unfiltered and unselected, and spgi-100 still shows 100 filtered rows
      with its 5 selected.
  - anchor: Step 5
    expectation: >-
      The "v ii" criterion narrows SPGI-linked2 itself: its filtered row count
      is strictly between 0 and its full 224 rows, and every surviving row
      carries "v ii" in link column 3 — the criterion is read back, not only
      the count. spgi-100, which is two links away and on neither end of this
      one, is read twice a settle apart and must still show 100 filtered rows
      with its 5 selected.
  - anchor: Step 6
    expectation: >-
      After applying the "v ii" category filter on link column 3 in the
      SPGI-linked2 Filter Panel, the SPGI-linked1 view shows exactly the
      DERIVED count of rows that are both key-matched to the selected spgi-100
      rows and carried by a "v ii" SPGI-linked2 row, because the
      FILTER_TO_FILTER link propagates the narrowed set from df3 back to df2.
      The frames that must NOT have moved are read twice, a settle apart:
      SPGI-linked2 still carries exactly the state Step 5 left it in (the link
      did not feed back into its own source) and spgi-100 is untouched.
  - anchor: Step 8
    expectation: >-
      After applying the "inconclusive" category filter on PAMPA Classification
      in the SPGI-linked1 Filter Panel, the SPGI-linked1 view shows exactly the
      DERIVED count of rows satisfying selection + "v ii" + "inconclusive" at
      once, and every surviving row carries "inconclusive". SPGI-linked2 is
      read twice, a settle apart, and must still carry exactly the state Step 5
      left it in — a FILTER_TO_FILTER link that also ran df2 -> df3 would
      contaminate it, which is the whole point of the link having a direction.
      spgi-100 is likewise unmoved.
  - anchor: Step 9
    expectation: >-
      Clearing the selection in spgi-100 releases only the selection-driven
      narrowing: SPGI-linked1's filtered row count rises strictly above the
      Step 8 count (so the selection link demonstrably contributed) and settles
      at exactly the DERIVED no-selection-link count — the number of
      SPGI-linked1 rows carried by a "v ii" SPGI-linked2 row AND classified
      "inconclusive", computed from the raw columns. That derived number, not a
      value read back out of the product, is what Step 11 is judged against.
      SPGI-linked2 is unmoved, and spgi-100 settles at 0 selected over its full
      100 rows.
  - anchor: Step 10
    expectation: >-
      Re-selecting the same 5 rows in spgi-100 — with the frame identity
      asserted at the gesture, so a switch that silently did not take cannot
      pass as "5 rows selected" — brings SPGI-linked1 back to exactly the Step
      8 derived count, proving the selection link is still live and that the
      Step 9 movement was the link releasing its grip rather than unrelated
      drift. SPGI-linked2 and spgi-100's own filtered count are read twice, a
      settle apart, and are unmoved.
  - anchor: Step 11
    expectation: >-
      GROK-19137. After the link type between spgi-100 and SPGI-linked1 is
      changed from "selection to filter" to "selection to selection" in the
      Data > Link Tables dialog, the dialog shows the new type and no editor
      still shows the old one; SPGI-linked1's SELECTION carries exactly the
      derived key-matched row count (the new link type is in effect); and
      SPGI-linked1's filtered row count equals the DERIVED no-selection-link
      count and is NOT the Step 8 count the old link type produced — the
      previous link type's narrowing was cleared rather than left stale.
      SPGI-linked2 is unmoved, and spgi-100 still holds its own 5-row selection
      over its full 100 rows.
---

# Collaborative Filtering for Linked Tables

Verifies that Filter Panel criteria propagate correctly across tables that are
linked with FILTER_TO_FILTER and SELECTION_TO_FILTER link types — both directions
of propagation and the composition of a link-driven filter with a panel-local one.

## Setup

1. Open the table **System:AppData/Chem/tests/spgi-100.csv** and note its name as
   **spgi-100**.
2. Open the table **System:AppData/ApiTests/datasets/SPGI-linked1.csv** and note its
   name as **SPGI-linked1**.
3. Open the table **System:AppData/ApiTests/datasets/SPGI-linked2.csv** and note its
   name as **SPGI-linked2**.
4. Link **SPGI-linked2** to **SPGI-linked1** using a FILTER_TO_FILTER link on all
   four key columns: Sample Name, link column 1, link column 2, link column 3.
5. Link **spgi-100** to **SPGI-linked1** using a SELECTION_TO_FILTER link on Id
   (spgi-100) and Concept Id (SPGI-linked1).

## Scenarios

### 1. Selection in one table propagates as a filter to a linked table

1. Switch to the **spgi-100** view and select the top 5 rows of the table.

   **Expected result (Step 1)**: exactly 5 rows are selected, and spgi-100's own
   filtered row count is still its full 100 rows. Selecting is not filtering — the
   master table of a selection-to-filter link must not narrow itself.

2. Switch to the **SPGI-linked1** view.

   **Expected result (Step 2)**: Exactly 9 rows remain visible in the
   SPGI-linked1 view. The selection-to-filter link propagated the selection from
   spgi-100 as a filter on SPGI-linked1.

3. Switch to the **SPGI-linked2** view.
4. Open the **Filter Panel** for this view.
5. In the **link column 3** filter card, select the category **v ii**.

   **Expected result (Step 5)**: SPGI-linked2 itself narrows — its filtered row
   count is above 0 and below its full 224 rows — and every row that survives the
   filter carries **v ii** in link column 3. The criterion is read back from the
   surviving rows, not inferred from the count alone.

6. Switch to the **SPGI-linked1** view.

   **Expected result (Step 6)**: Exactly 5 rows remain visible in the
   SPGI-linked1 view. The filter-to-filter link propagated the narrowed set from
   SPGI-linked2 down to SPGI-linked1, combining with the selection-driven filter
   from spgi-100.

7. Open the **Filter Panel** for the **SPGI-linked1** view.
8. In the **PAMPA Classification** filter card, select the category **inconclusive**.

   **Expected result (Step 8)**: Exactly 2 rows remain visible in the
   SPGI-linked1 view. The criterion set in this panel combined with the
   link-driven narrowing from SPGI-linked2, so both conditions apply at once.

### 2. Changing the link type clears the previous link type's effect (GROK-19137)

This scenario continues from the state Scenario 1 leaves behind: spgi-100 has 5 rows
selected, SPGI-linked2 carries the **v ii** criterion, SPGI-linked1 carries the
**inconclusive** criterion, and SPGI-linked1 shows 2 rows.

9. Switch to the **spgi-100** view and clear its selection.

   **Expected result (Step 9)**: SPGI-linked1's filtered row count rises strictly
   above the 2 rows of Step 8 — so the selection link was demonstrably contributing
   narrowing — while staying strictly below its full 3624 rows, because the panel
   criterion and the filter-to-filter link from SPGI-linked2 still apply. Record this
   value as the **no-selection-link baseline**: it is what SPGI-linked1 must show
   whenever nothing is arriving over the spgi-100 link.

10. Select the same top 5 rows in **spgi-100** again.

    **Expected result (Step 10)**: SPGI-linked1 returns to exactly 2 filtered rows.
    This proves the selection link is still live and that the Step 9 movement was the
    link releasing its grip, not unrelated drift — without it, "the narrowing is gone"
    in the next step could not be told apart from "the link was never working".

11. Open **Data > Link Tables...**, select the tab for the **spgi-100 -> SPGI-linked1**
    link, and change its **Link Type** from `selection to filter` to
    `selection to selection`. Close the dialog.

    **Expected result (Step 11)**: this is the GROK-19137 surface.

    - The dialog reflects the change: a Link Type editor now reads
      `selection to selection`, and no editor is left reading `selection to filter`.
    - The new link type is in effect: SPGI-linked1's **selection** holds exactly the 9
      key-matched rows — the same 9 rows the old link type used to push onto the filter.
    - The old link type's effect is CLEARED, not left stale: SPGI-linked1's filtered row
      count is back at the no-selection-link baseline recorded in Step 9, and is NOT the
      2 rows the selection-to-filter link produced. A stale link leaves the target still
      carrying the previous type's narrowing — that is the bug.

## Teardown

In a finally block that runs even on failure: close the Link Tables dialog if it is
still open, then close all probe views.

**Known limitation, declared:** "closing the views drops the probe tables and, with them,
their links" is the intent, not something the teardown verifies. There is no JS API handle
for a link — no getter, no remover — so the only thing the teardown can assert is what it
does assert: the Link Tables dialog is gone and no table view is left open. Whether
the link registry is actually emptied by that is unobservable from the spec, and a leaked
link would surface only as cross-talk in a later spec that happens to reopen one of the
three fixtures under the same name.

## Automation notes

The Link Tables dialog anatomy, the absence of a link getter, the fact that `grok.data.linkTables`
only ever creates, the `attach()`/`detach()` mechanism behind a link-type change and the
close-the-views teardown are in `.claude/skills/grok-browser/references/viewers/filters.md`.

- **Setup builds the three-table linked graph programmatically** — `grok.data.files.openTable`
  plus `grok.data.linkTables` with `DG.SYNC_TYPE.FILTER_TO_FILTER` and
  `DG.SYNC_TYPE.SELECTION_TO_FILTER` — before the UI steps begin.
- **Every expected row count is DERIVED from the raw columns, never pinned as a literal.** Setup
  reads spgi-100's first five `Id` values, the key tuples of the SPGI-linked2 rows carrying `v ii`,
  and SPGI-linked1's `Concept Id` / `PAMPA Classification` columns, and from those computes the four
  counts the steps assert: the key-matched count (Step 2), the key-matched-and-`v ii` count (Step 6),
  the count that adds `inconclusive` (Steps 8 and 10), and the no-selection-link count (Steps 9 and
  11). Each is then bracketed strictly inside its frame before use, so a fixture edit that made one
  of them 0 or the whole frame fails Setup instead of silently retuning what "the link works" means.
  In particular the Step 11 baseline is derived, not read back out of the same product state the
  assertion re-reads.
- **Every propagation step also reads the frames that must NOT have moved**, twice, a settle apart,
  so a bystander clause is a sustained zero rather than a snapshot taken before contamination could
  arrive. Without that half, a FILTER_TO_FILTER link running in both directions — exactly the
  cross-contamination link direction exists to prevent — leaves every step green.
- **Measurement limitation, declared:** link correctness cannot be asserted directly, so it is
  proved by the links' EFFECT at Steps 2, 6 and 8. A separate "link graph established" precondition
  step was removed as unverifiable; its reachable part (three named table views open with their
  expected row counts) is in Setup.
- **Selecting the top 5 rows of spgi-100 (Step 1)** is a canvas gesture: Shift+click row index 0,
  then Shift+click row index 4 — or, programmatically,
  `df1.selection.setAll(false); for (let i = 0; i < 5; i++) df1.selection.set(i, true)`.
- **Category selection (Steps 5 and 8) is driven programmatically** through
  `fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'link column 3', selected: ['v ii']})`
  and the same shape for Step 8. Because this is not a gesture on the card, each of these steps
  reads the criterion back off the surviving rows, so a criterion that silently failed to apply
  cannot pass as a count.
- **Both actuation paths were measured for Steps 5 and 8, and they agree.** A card click raises one
  notification the API path does not — `FILTER_CRITERIA_CHANGED` on the frame's own event bus — but
  a link never listens to it: it keys off `onFilterChanged` / `onSelectionChanged` and contributes
  through `onRowsFiltering` (`ddt/lib/src/data_frame/tools/link.dart:86-106,159-186`), which both
  paths raise. Measured on dev 2026-08-18 with the criteria set by clicking the category-name rows
  instead: SPGI-linked2 148 rows and SPGI-linked1 9 → 5 for Step 5, SPGI-linked1 2 rows for Step 8
  with SPGI-linked2 unmoved at 148 and spgi-100 unmoved at 100 filtered / 5 selected across a 4 s
  window. These are the counts the steps already assert, so the substitution measures the mechanism
  the user's gesture travels and not a weaker one.
- **The "rows remain visible" checks (Steps 2, 6, 8)** read
  `grok.shell.tv.dataFrame.filter.trueCount` on the **SPGI-linked1** view — the count of rows
  passing the composed filter, which is what that grid renders.
- **Readiness is never a fixed delay** (per GROK-19152): the Filter Panel barrier keys on the
  panel element appearing in the DOM, and a linked table's recount is polled for.
- **Changing a link's type (Step 11) has no JS API path** — the only path is the
  **Data > Link Tables…** dialog.
- **Teardown:** close the Link Tables dialog if still open, then close all probe views.
