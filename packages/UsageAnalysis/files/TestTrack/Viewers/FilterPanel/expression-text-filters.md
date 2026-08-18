---
feature: filters
realizes_atlas:
  - filters.cp.expression-and-text-ui
realizes:
  - viewers.filters
  - viewers.filters.free-text
priority: p2
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-20242
    status: fixed
scope_reductions:
  - id: SR-01
    check: E-TRACE-02
    rationale: |
      Scenario 2 expects the Filter Panel to already be presenting an Aroma
      text-filter card when it opens on beer.csv, and Step 15 asks for the
      two terms to be cleared in place. The spec instead reaches that state
      through two other panel affordances: the hamburger menu's Remove All,
      then the panel header's column selector on Aroma — which builds the
      SAME card, because a categorical column carrying the Text semantic
      type gets a text filter as its default (filter_control.dart:156).
      Nothing about the card is assumed: the spec waits until the Aroma card
      carries a `.d4-text-filter` body, its Fuse update indicator has been
      taken down (the liveness condition Step 12 asks for, in place of the
      fixed delay it forbids) and both its search input and Fuzzyness slider
      are in the DOM, and it then asserts the full beer row count before
      typing anything.
    verdict_status: SCOPE_REDUCTION
  - id: SR-02
    check: E-TRACE-02
    rationale: |
      Step 8 quotes the free-text expression as `HEIGHT < 160`. That string
      is NOT the card's expression language and is silently rejected: the
      card itself builds `${HEIGHT} < 160` (structure_filter.dart:165), and
      typing the bare column name leaves the row count exactly where the
      committed form rule already had it — which is why an equal-count check
      there passes without the input having done anything. The spec types
      the real syntax, asserts the input is not flagged invalid, and drives
      a DIFFERENT predicate first: `${HEIGHT} < 150`, ANDed with the
      committed HEIGHT < 160 rule (the card is still in the AND logic Step 6
      left it in), which must move the count off the single-rule value onto
      a second count derived off the HEIGHT column.
      The "return to the original" half is driven from a cleaned card rather
      than on top of the wider predicate, because Enter COMMITS a free-text
      query as a further term and a second typing can therefore never take
      the count back down on the same card. The spec clears the criteria
      with the panel header's reset icon, proves the card is carrying
      nothing (the full demog row count) and is still in free-text mode, and
      types the original predicate back in — so the returned count is
      measured from a proven, different pre-state.
    verdict_status: SCOPE_REDUCTION
expected_results:
  - anchor: Step 3
    expectation: >-
      df.filter.trueCount drops below the full demog row count after the first
      expression rule is committed, and equals the AGE > 50 row count computed
      independently off the AGE column — a drop alone does not say the rule was
      the criterion. fg.getFilterSummary() names one rule.
  - anchor: Step 5
    expectation: >-
      After the second rule is committed in AND mode, df.filter.trueCount drops
      further (or stays, depending on row overlap). fg.getFilterSummary() names
      two rules.
  - anchor: Step 6
    expectation: >-
      The two-rule card reads OR in its header before the toggle is clicked, and
      reads AND after — the toggle is confirmed through the header text, not
      assumed. df.filter.trueCount is greater than or equal to the AND value: the
      ordering trueCount_AND <= trueCount_OR holds. For this particular pair of
      overlapping-but-distinct rules the inequality is also STRICT
      (trueCount_AND < trueCount_OR), so a toggle that did nothing at all could
      not satisfy it.
  - anchor: Step 7
    expectation: >-
      After Remove Query on the first rule, df.filter.trueCount rises and
      fg.getFilterSummary() names one rule.
  - anchor: Step 8
    expectation: >-
      The card's own expressionMode reads `free-text` after the icon-italic
      toggle. After entering an equivalent expression, df.filter.trueCount
      matches the form-built single-rule trueCount recorded before the mode
      switch.
  - anchor: Step 9
    expectation: >-
      The card's own expressionMode reads `expression` again after the second
      icon-italic click, and df.filter.trueCount is unchanged.
  - anchor: Step 10
    expectation: >-
      GROK-20242 — with the rule built on the string column `USUBJID` (the regexp
      operation is offered on string columns only, so it cannot be built on the
      numeric AGE column), pasting a comma-separated list with the regexp
      operation selected rewrites the Value input to the pipe alternation
      `5|15|25`, and the committed rule text carries that alternation.
      df.filter.trueCount is a real subset (strictly between 0 and 5850).
      Appending a trailing comma leaves the rewritten value and trueCount
      identical.
  - anchor: Step 11
    expectation: >-
      Each of the three non-numeric operators builds the row set it names:
      `SEX equals F`, `RACE contains an` and `STARTED after 01/01/1991` each
      commit exactly one rule and leave df.filter.trueCount equal to the row
      count computed independently off that column for that predicate, with the
      filter summary naming the column. Each expected count is asserted to be
      strictly between 0 and 5850 BEFORE the comparison, so an operator that did
      nothing could not satisfy it.
  - anchor: Step 12
    expectation: >-
      The perturbation is proved before the restore: the card reads OR before the
      toggle and AND after, and in AND mode the two-rule count is strictly below
      the count the AGE rule gives on its own, so the suspension has something to
      undo. Unchecking a rule's own checkbox then suspends it without removing
      it: the card still lists two rules, the card's gridValues read [true, false], and
      df.filter.trueCount returns to exactly the count the remaining rule gives
      on its own. Re-checking restores the exact two-rule AND count. This is what
      separates a suspended rule from a removed one (Step 7).
  - anchor: Step 13
    expectation: >-
      With Fuzzyness pinned at 0, typing a term into the Aroma text-filter card
      search input and pressing Enter drops df.filter.trueCount below the full
      beer.csv row count and fg.getFilterSummary() reflects the active term. The
      surviving rows are also the RIGHT rows, counted in both directions: at
      least one passing row's Aroma contains the term, no passing row's Aroma
      lacks it, and no excluded row's Aroma contains it.
  - anchor: Step 14
    expectation: >-
      The two-term Aroma card reads OR in its header before the toggle and AND
      after it, and trueCount in AND mode (trueCount_and_beer) is strictly less
      than trueCount in OR mode (trueCount_or_beer) for the same two terms —
      trueCount_and_beer < trueCount_or_beer.
  - anchor: Step 15
    expectation: >-
      At fuzziness 0, the one-character near-miss term `maltx` matches no Aroma
      value and df.filter.trueCount == 0. After raising the Fuzzyness slider,
      df.filter.trueCount rises above 0 and grows further as the threshold
      increases. The term must be a near-miss of a real Aroma term: a far-off
      term such as `zzzzz` sits at 0 at every slider position, which would make
      the second half of this check unable to fail.
realized_as:
  - expression-text-filters-spec.ts
---

# Filters — Expression filter and Text filter driven through their own UI

## Setup

1. Open demog.csv and open a table view.
2. Open the Filter Panel via the toolbar.
3. Note the total row count shown in the Filter Panel header before any filter is applied (this is the starting row count for demog).

## Scenarios

### Scenario 1: Expression filter — form builder, AND/OR toggle, free-text mode, GROK-20242 paste

Steps:
1. Open the Filter Panel hamburger (context) menu and choose **Add filter > Expression filter** to add an Expression filter card.
2. In the Expression filter card, set Column to `AGE`, Operation to `>`, Value to `50`, and click the **(+)** Add-filter button to commit the rule.
3. Verify that the filtered row count dropped below the starting demog row count and that it equals
   the number of rows with `AGE > 50` worked out independently from the AGE column — a drop on its
   own does not say the rule was the criterion. Verify the filter summary names one rule. Note this
   filtered row count.
4. In the same card, set Column to `HEIGHT`, Operation to `<`, Value to `160`, and click **(+)** to commit the second rule in AND mode. Note the filtered row count.
5. Verify that the filter summary names two rules and the filtered row count dropped (or held) relative to the count after the first rule.
6. Verify the card header reads **OR** — a two-rule card opens in OR logic — and note the filtered
   row count in that mode. Click the header toggle and verify the header now reads **AND**; note the
   filtered row count again. Verify that the AND filtered row count is less than or equal to the OR
   filtered row count — the AND result is always an equal or stricter subset of the OR result for the
   same two rules — and that for this particular pair it is STRICTLY less, so a toggle that did
   nothing could not pass.
7. Right-click the first expression rule row in the card and select **Remove Query** from the context menu. Verify that the filtered row count rises back toward the single-rule result and the filter summary names one rule.
8. Click the italic/list icon on the Expression filter card header to switch to **free-text mode**. Verify the card reports itself as being in free-text mode. Type the expression `HEIGHT < 160` in the search input and press **Enter**. Verify that the filtered row count matches the count noted after the first rule — the mode switch is lossless for an equivalent expression.
9. Click the italic/list icon again to return to **form mode** (round-trip revert). Verify the card reports itself as being back in form mode and the filtered row count is unchanged.
10. In the form, set Column to `USUBJID` and Operation to **regexp** — the regexp operation is
    offered on string columns only, so it cannot be built on the numeric `AGE` column used in the
    steps above. Paste the string `5,15,25` into the Value input. Verify the Value input now reads
    the pipe alternation `5|15|25` and that the committed rule text carries that alternation.
    Verify that the filtered row count equals the number of rows whose `USUBJID` matches any of
    `5`, `15` or `25`, and that it is a real subset — neither zero nor the whole table. Then paste
    the same list with a trailing comma (`5,15,25,`) and verify the rewritten value and the
    filtered row count do not change — per GROK-20242, pasting a comma-separated list with regexp
    selected rewrites the list into a pipe-separated alternation and trims any trailing comma, so
    the trailing comma must be a no-op.

11. The operation dropdown carries string and date operators as well as the numeric comparisons used
    above, and those reach different evaluators. Build one rule for each on a fresh card and check
    it against a row count worked out independently from the column itself: `SEX equals F`,
    `RACE contains an`, and `STARTED after 01/01/1991`. For each, verify the card committed exactly
    one rule, the filtered row count equals the independently computed count, and the filter summary
    names that column. Before each check, confirm the independently computed count is neither zero
    nor the whole table — otherwise the comparison would be satisfied by an operator that does
    nothing at all.
12. A rule can also be switched off from its own checkbox, which is a different thing from removing
    it. On a fresh card commit `AGE > 50` (note the count it gives on its own), then `HEIGHT < 160`,
    and switch the card to AND so both rules bite — verify the header read OR before the click and
    AND after it, and that the two-rule count is strictly below the AGE-alone count, so the
    suspension below has something to undo; note that two-rule count. Uncheck the second
    rule's checkbox in the rule list. Verify: the card still lists two rules — the rule was
    suspended, not deleted — its own enabled flag reads off, and the filtered row count is exactly
    the count the AGE rule gave on its own. Re-check it and verify the two-rule count comes back
    exactly.

Expected:
- Each (+) commit drops the filtered row count and updates the filter summary rule count.
- Each of `equals`, `contains` and `after` selects exactly the row set computed independently for
  it, and that row set is neither empty nor the whole table.
- Unchecking a rule leaves the rule in the list, marks it as off, and returns the row count to what
  the remaining rule gives alone; re-checking restores the two-rule count exactly.
- The AND filtered row count is less than or equal to the OR filtered row count for the same pair of rules.
- After Remove Query, the summary names one rule and the filtered row count rises.
- Free-text mode and form mode yield equal filtered row counts for equivalent expressions; round-trip revert is lossless.
- Pasting `5,15,25,` and `5,15,25` produce identical filtered row counts.

### Scenario 2: Text filter — OR/AND modes and fuzzy search on beer.csv

Steps:
11. Open beer.csv and open a table view. Open the Filter Panel on it.
12. Wait for the Text filter card to become interactive before proceeding — the fuzzy-search engine initialises asynchronously; do not rely on a fixed delay.
13. Note the total row count shown in the Filter Panel header before any filter is applied (the starting row count for beer). Set **Fuzzyness** to 0 first, so the match is exact and the content check below is meaningful. In the **Aroma** column text-filter card search input, type `malt` and press **Enter**. Verify that the filtered row count drops below the starting beer row count and the filter summary reflects the active term. Then check WHICH rows survived, not only how many: every row still passing must have `malt` in its Aroma value, and no row containing `malt` may have been excluded. A row count on its own is satisfied by any filter that happens to keep the same number of rows.
14. Add a second term by typing `hop` and pressing **Enter**. Verify the card header reads **OR** — a two-term card opens in OR — and note the filtered row count in that mode. Switch the AND/OR control in the card to **AND** mode and verify the header now reads **AND**. Note the filtered row count in AND mode. Verify that the AND filtered row count is strictly less than the OR filtered row count — a single-mode reading proves nothing about the control being wired.
15. Clear both terms. Set **Fuzzyness** to 0 and type `maltx` — a one-character near-miss of the real Aroma term `malt` — and verify the filtered row count is 0. Raise the **Fuzzyness** slider to a non-zero position. Verify the filtered row count rises above 0 and increases further as the slider moves higher. The term has to be a near-miss rather than a far-off string such as `zzzzz`, because a term no slider position can recover would hold the count at 0 throughout and make the "raising fuzziness raises the count" half unable to fail — it would pass on a broken fuzzy engine.

Expected:
- Typing a term and pressing Enter drops the filtered row count; the filter summary reflects the term; at fuzziness 0 the surviving rows are exactly the rows whose Aroma contains the term, counted in both directions.
- The AND filtered row count is strictly less than the OR filtered row count for the same two terms.
- At fuzziness 0 a near-miss term yields a filtered row count of 0; raising fuzziness raises the filtered row count above 0 and grows it further.

## Teardown

Close every open table view (demog and beer) so the next scenario starts from an empty shell.

## Automation notes

- SIGNAL CAUTION: the Expression filter rule list and the Text filter term list are painted by an internal GridCore canvas — they are not DOM nodes. Assert them through `df.filter.trueCount` and `fg.getFilterSummary()`; a settle-gated canvas diff of the card region plus a no-error floor is an acceptable secondary signal. Never assert rule/term presence by DOM node count.
- Remove Query in Scenario 1 Step 7 is a CLASS-2 CANVAS GESTURE — how the context menu resolves the
  row it acts on, and the actuation that replaces it, are in
  `.claude/skills/grok-browser/references/viewers/filters.md`.
- Fuse.js applies to the TEXT HALF ONLY — the expression half neither loads Fuse nor has a Fuzzyness slider, so waiting for Fuse there would wait forever.
- Dataset split is intentional: the text half cannot run on demog (no `Text`-semType categorical column) and beer.csv's Aroma column qualifies. The expression half runs on demog.
- GROK-20242 (paste for regexp) is exercised in Scenario 1 Step 10, grounded on the string column `USUBJID` — the regexp operation is not offered on numeric columns, so it cannot be driven on AGE. The invariants asserted (alternation rewrite, real subset, trailing comma is a no-op) are the ones the ticket names.
- Step numbering: Scenario 2 continues the original single-run numbering and starts at 11, so its
  first two steps carry the same numbers as Scenario 1 Steps 11 and 12. The `expected_results`
  anchors `Step 11` and `Step 12` refer to the SCENARIO 1 steps (the operator sweep and the per-rule
  checkbox); anchors `Step 13`-`Step 15` refer to Scenario 2. Resolve an anchor by its content, not
  by position.
- Step 15 term choice: the term must be a one-character near-miss of a real Aroma value (`maltx`
  for `malt`), which yields 0 at fuzziness 0 and grows as the threshold rises. A far-off term such
  as `zzzzz` cannot be recovered at any slider position, so the count would stay at 0 throughout
  and the "raising fuzziness raises the count" half would be unfalsifiable.
