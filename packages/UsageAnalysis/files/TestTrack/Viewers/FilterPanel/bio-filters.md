---
feature: filters
realizes_atlas:
  - filters.cp.chem-and-bio-filters
realizes:
  - viewers.filters
priority: p2
target_layer: playwright
coverage_type: regression
related_bugs: []
realized_as:
  - bio-filters-spec.ts
scope_reductions:
  - id: SR-01
    check: E-TRACE-02
    rationale: |
      Step 3 asks for the Macromolecule column to be checked BY NAME in the
      "Select columns..." picker. That picker renders its column list as a
      canvas grid: per-row checkbox clicks do not register (measured across 4
      x-offsets x 28 y-positions — the "<n> checked" counter never moved), and
      the picker's own search box does not narrow what the "All" link checks.
      The only DOM-drivable actuation is All -> OK. The spec therefore drives
      All -> OK and proves the column binding on the far side instead: the
      counter is read immediately before and after the All click and must move
      from exactly 0 to at least 1, and the card the OK creates is then
      asserted on the CARD channel — exactly one card in the panel whose own
      caption is the Macromolecule column and which carries a Substructure
      input. The filter object's own column name is NOT the channel: a Bio
      substructure filter in the Filters group reports an empty column name and
      an empty type on this build (see the measurement limitation in the
      Automation notes), so filterColumnName cannot carry the claim. The group
      is read as a second channel by SIZE only — exactly one more filter after
      the menu leaf than before it — which is what a card drawn with no filter
      behind it fails. So "the Macromolecule column was offered and checked" is
      asserted on the created card's caption rather than on the picker row, and
      the binding is NOT asserted through the filter object.
    verdict_status: SCOPE_REDUCTION
expected_results:
  - anchor: "Scenario 1, Step 3"
    expectation: >-
      Any pre-existing substructure card is first removed through its own X
      icon and the panel is confirmed to hold zero Bio substructure cards bound
      to the Macromolecule column, so the menu drive that follows cannot be
      credited with a card it did not create. Add Filter > Bio Substructure
      Filter... is then actuated unconditionally through the panel hamburger;
      its "Select columns..." picker is read as exactly 0 checked before the
      All click and at least 1 after, and confirming it with OK takes the panel
      to exactly one Bio substructure card CAPTIONED with the Macromolecule
      column and carrying a Substructure input, while the Filters group grows by
      exactly one filter. The filter object's own column name is empty on this
      build and carries nothing (scope reduction SR-01), so the caption is the
      channel the binding is read on and the group is read for size only.
      Entering the subsequence T-T-Y-K-N-Y-V in that
      card's Substructure input and applying it drops df.filter.trueCount below
      the full peptides row count while keeping at least one peptide row — a
      subsequence that matched nothing would make the Step 4 round-trip
      unreadable.
  - anchor: "Scenario 1, Step 4"
    expectation: >-
      Immediately before the input is cleared the filter is confirmed to still be
      narrowing (trueCount strictly below the full peptides row count), so the
      restore below is falsifiable. After clearing the Bio substructure filter's
      input, df.filter.trueCount returns to the full peptides row count exactly,
      confirming the filter releases rows when cleared.
---

# Filters — Bio package filters

## Setup

The Bio half needs no Chem dataset and no Chem readiness decision. It opens its own table and
takes its own package-availability decision before any step runs.

That decision has **three distinguishable outcomes**, and they must not collapse into fewer:

- **Bio is not installed** — no `bioSubstructureFilter` function owned by the `Bio` package is
  registered within 20 s. That is a declared precondition failure: the scenario is skipped, because
  there is no subject to test.
- **Bio is installed but `seqHelper` has not come up yet** — the function IS registered and
  `apply()` throws `Package Bio .seqHelper is not initialized`. That is a startup race, not a
  broken build: the readiness probe keeps polling and only gives up when the window runs out. It is
  never a failure on first sighting.
- **Bio is installed but its substructure filter is broken** — the function IS registered and
  `apply()` throws anything OTHER than the not-initialized message, or keeps resolving to nothing
  until the readiness window expires. That is the regression this scenario exists to catch, so it is
  a plain assertion failure carrying the raised error, never a skip and never a retry.

The distinction between the second and third outcomes is **time**, not error class: the same
not-initialized message is retryable early and terminal at the deadline. The failure message
therefore reports how long the probe waited, how many not-initialized sightings it tolerated, and
whether it stopped on a fatal error or ran the window out, so a future reader can tell a genuinely
uninitialized build from a window that was simply too short.

**The readiness window is 60 s, measured on dev 2026-08-18.** Across 14 sequential runs the factory
was registered at 1-5 ms every time and became applicable at a median of 365 ms, worst case
3 856 ms (full set, ms: 3, 1 961, 3, 259, 471, 597, 3 856, 2 924, 1 015, 6, 58, 100, 5, 484). Every
one of those runs recorded zero retries of any kind, so the multi-second cases are a single slow
`apply()` call, not a package that was still coming up. The
`seqHelper` race did not reproduce in those 14 runs, so the window is not set from its distribution
— it is set from the two facts that are known: healthy readiness costs under 4 s, and the one
observed race was still uninitialized after the previous 20 s window had run out, so 20 s is
demonstrably too short. 60 s is 3x the window known to be insufficient and ~15x the worst healthy
case, and it bounds a genuinely broken build to a 60 s failure well inside the 600 s test timeout.
Registration is gated separately at 20 s, so an absent package still skips quickly instead of
burning the full readiness window.

## Scenarios

### Scenario 1: Bio substructure filter — subsequence filtering on peptides

Steps:

1. Open the file **System:DemoFiles/bio/peptides.csv**. Note the total row count shown
   in the table header.

2. Open the Filter Panel on this table by clicking the funnel icon in the ribbon toolbar.
   Wait for the filter cards to appear.

3. Ensure the Bio package is loaded. If the panel already carries a substructure filter card,
   remove it with the card's own **X** icon and confirm the panel now holds none — the menu gesture
   below must be seen to create the card, not to find one. Then open the hamburger menu →
   **Add Filter** → **Bio Substructure Filter...**, press **All** in the **Select columns...**
   dialog and press **OK**. Wait for the bio substructure filter card to appear.
   - In the **Substructure** input field of the bio filter card, enter `T-T-Y-K-N-Y-V`
     and confirm the input (press Enter or click **Apply** if a button is present).
   - Verify the **Select columns...** picker read exactly zero columns checked before **All** was
     pressed and at least one after. A picker confirmed with nothing checked creates no card at all,
     and the step would then fail on a missing input rather than on the filtering under test.
   - Verify the panel went from zero to exactly one bio substructure card, and that the card is
     bound to the Macromolecule column.
   - Wait for the Filter Panel header row count to stabilize.
   - Verify the Filter Panel header shows fewer rows than the total row count noted in
     Step 1 (the subsequence filter is active and reduces the visible peptide rows), and that at
     least one peptide row still passes — a subsequence matching nothing would make the Step 4
     round-trip unreadable.

4. Clear the bio substructure filter input (select all text and delete, or click a clear
   button if present). Wait for the Filter Panel header row count to stabilize.
   - Immediately BEFORE clearing, verify the filter is still narrowing (row count below the total).
     Without that, "all rows came back" is unfalsifiable — they may never have left.
   - Verify the Filter Panel header shows the same total row count as noted in Step 1
     (all peptide rows returned when the filter is cleared — the round-trip revert).

Expected:
- Step 3: the panel went from zero bio substructure cards to exactly one bound to the Macromolecule
  column, the picker moved from 0 checked to at least 1 across the **All** click, and the Filter
  Panel header shows fewer rows than the full peptide count — but more than zero — after the
  subsequence is entered.
- Step 4: the filter is confirmed to still be narrowing immediately before the input is cleared, and
  the Filter Panel header shows the full peptide row count after clearing the filter.

## Automation notes

The Bio card anatomy, the two readiness gates, the peptides row counts and the column-picker
actuation are in `.claude/skills/grok-browser/references/viewers/filters.md`.

- **This scenario and `chem-filters.md` are a deliberate keep-both pair.** Both realize
  `filters.cp.chem-and-bio-filters`, so a coverage pass will propose consolidating them. It must
  not: Gate B has no skip channel, and while both halves lived in one spec file a skipped Chem
  half beside a passing Bio half was recorded as a PASS — which is how the Chem scenario ran zero
  steps for weeks while its file was green. Separate files make either half's skip visible as a
  skip. The disposition for any redundancy finding on this pair is **keep both**.
- **Measurement limitation, declared:** a Bio substructure filter object in the Filters group
  reports an EMPTY column name and an EMPTY type on this build (measured on dev 2026-08-18 — the
  group held exactly one filter and it read `{"column":"","type":""}`), so the filter's column
  binding cannot be asserted through that channel. The binding is asserted on the card instead: the
  card's own caption must be the Macromolecule column and the card must carry the Substructure
  input. The group is still read as a second channel, as a count delta — it must hold exactly one
  more filter after the menu leaf than before it, so a card drawn without a filter behind it fails.
- **The table is opened through the shared open-table helper** (which closes previously open views,
  sets the animation-suppression flag and waits out the Bio/Chem semantic-type settle), but the
  Filter Panel is NOT opened through the shared open-panel helper: that helper's readiness barrier
  is the first filter card, and peptides.csv opens the panel with **zero** cards until the Bio
  package finishes loading — measured, the barrier burns its full 15 s and hard-fails before any
  readiness check can report the real condition. The barrier here is the panel element plus the
  Filters viewer's own props, followed by a card-count stability barrier.
- **Both readiness gates run BEFORE the panel is opened, and this ordering is load-bearing**
  (measured on dev 2026-08-18). Loading Bio is what makes the panel auto-populate — peptides.csv
  gets an ID card, an IC50 card and an **auto-added bio substructure card on the Macromolecule
  column**. Gate first and the panel is built once, deterministically, with all three cards present
  before the gesture; open the panel first and the auto-population lands at an arbitrary later
  moment. In the observed failure it landed between the two steps and inserted a SECOND, empty
  AlignedSequence card ahead of the one the scenario created, so the step that clears "the"
  Substructure input addressed the empty one and the round-trip could not complete. Three runs in
  ten failed that way before the reorder; four in four passed after it. This auto-added card is also
  exactly why the card-adding gesture must never sit behind an "is there already a substructure
  input?" guard.
- **Teardown:** close probe tables in `finally`, even on failure — the `finally` wraps everything
  from the table open onward, so the skip path and any hard throw both still close the table.

### Spec must keep

1. **The column-picker flow in the card-adding step, driven UNCONDITIONALLY** — remove any
   pre-existing substructure card, assert zero, then leaf → picker → All → OK → picker confirmed
   gone → exactly one card bound to the Macromolecule column. Not `fg.updateOrAdd` as a substitute
   for the menu gesture, and never behind an "is there already an input?" guard: a build that
   auto-adds the card would then never actuate the menu leaf at all, and a total regression of that
   leaf would be invisible. `input[placeholder="Substructure"]` is not Bio-specific — Chem's
   substructure filter uses the same placeholder — so card counts are scoped to the Macromolecule
   column, never to the placeholder alone.
2. **Dialogs addressed unambiguously** — the picker by `name=`. No bare `.d4-dialog` anywhere.
3. **The package-readiness guard stays OUTSIDE every `softStep`.** A `test.skip()` inside a
   step wrapper is swallowed and the file reports a non-failure with zero assertions executed. It
   also stays a guard on ONE condition: "no Bio-owned function registered" skips, while a registered
   function that never becomes applicable fails. A guard that swallows both retires the regression
   under test as a non-failure. The not-initialized retry must stay matched on that ONE message: a
   blanket `catch` that retried every error would make "broken" and "still starting" indistinguishable
   again, which is exactly the collapse this guard exists to prevent. Verified on dev 2026-08-18 by
   injection — a foreign error message fails at 4 ms with no retries, a persistent not-initialized
   message retries 118 times across 60 s and then fails, and an unregistered package skips at 20 s.
4. **The positive half of Step 4 is asserted immediately before the input is cleared.** Without it,
   "all rows came back" is satisfied by rows that never left.
