---
feature: chem
realizes_atlas:
  - chem.cp.substructure-search-top-menu
realizes:
  - chem.cp.substructure-search-top-menu
priority: p0
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: smoke
related_bugs: []
realized_as:
  - chem-substructure-search-top-menu-spec.ts
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-08-21-chem-new-02
    timestamp: "2026-08-21T00:00:00Z"
    failure_keys: []
    review_round: 1
  b:
    verdict: PASS
    cycle_id: direct-gate-b-2026-08-22-chem-substructure-search-top-menu-r2
    timestamp: "2026-08-22T14:26:00Z"
    spec_runs:
      - spec: chem-substructure-search-top-menu-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 22
        failure_keys: []
        run_mode: headless-cold
  e:
    verdict: PASS
    cycle_id: direct-gate-e-2026-08-22-chem-substructure-search-top-menu-r3
    timestamp: 2026-08-22T14:31:55Z
    failure_keys: []
expected_results:
  - anchor: "Scenario 1, Step 2"
    expectation: >-
      Driving Chem | Search | Substructure Search from the top menu does NOT
      search by itself: it adds exactly one substructure filter card for the
      molecule column to the Filter Panel, seeded EMPTY (its saved molBlock is
      DG.WHITE_MOLBLOCK and the filter reports it is not filtering), and opens
      that card's sketcher for the user to draw in. The card is confirmed
      present and its query empty before any structure is entered, and
      df.filter.trueCount still equals the Setup row count. On this FIRST
      invocation no column-picker dialog appears and exactly one sketcher opens,
      because the Molecule-typed column list at that moment is exactly
      [canonical_smiles] — pinned in Setup, not assumed — and the editor
      branches on that count. That holds for the first invocation only: the same
      list is re-measured in Scenario 2, Step 1 and is longer by then.
  - anchor: "Scenario 1, Step 3"
    expectation: >-
      After the query SMILES is entered into the filter card's sketcher and
      committed, the table filters to matching rows. The filtered row count
      (df.filter.trueCount) is strictly greater than zero and strictly less than
      the total row count — the filter is active and at least one row passes and
      at least one row is excluded — and no row is deleted. The surviving set is
      exactly the substructure predicate, checked row by row: every surviving
      row contains the query fragment and no excluded row contains it, with
      every molecule in the column parsing.
  - anchor: "Scenario 1, Step 4"
    expectation: >-
      After a non-matching query is applied in the same sketcher,
      df.filter.trueCount is exactly zero while the dataset row count is
      unchanged (no rows deleted). Zero is confirmed to be the correct answer
      rather than a silently failed search: checked row by row, no molecule in
      the column contains the query fragment.
  - anchor: "Scenario 2, Step 1"
    expectation: >-
      Re-invoking Chem | Search | Substructure Search after a search has already
      run opens the column-picker dialog titled "Substructure search" rather
      than the sketcher. The reason is measured, not assumed: the Molecule-typed
      column list is re-read immediately before the re-invocation and is longer
      than the single-entry list Setup pinned, because the searches already run
      cached a Molecule-typed helper column (see Automation notes). Confirming
      the picker leaves the filter on the user-visible canonical_smiles rather
      than on the cached helper column, and re-prepares it empty: the Filter
      Panel still holds exactly one substructure filter card for the column —
      the count is unchanged against a card count that is confirmed non-zero, so
      re-invoking replaces rather than stacks — and df.filter.trueCount returns
      to the full row count before the new query is entered.
  - anchor: "Scenario 2, Step 2"
    expectation: >-
      After the picker is confirmed and a new query entered that matches a
      different subset of rows, the SAME filter card is reused and updated
      rather than a second one added: df.filter.trueCount reflects the new match
      count, strictly greater than zero and strictly less than the total row
      count, and the new match count differs from the count produced by the
      first query in Scenario 1 Step 3 — demonstrating that re-invoking the
      command from the top menu replaces the previous filter rather than
      stacking it. Row by row, every surviving row contains the new fragment and
      no excluded row contains it, so the previous filter is gone rather than
      intersected with the new one.
---

# Chem — Substructure Search via top menu

## Setup

1. Open the file **System:DemoFiles/chem/smiles.csv**. Note the full row count displayed
   in the table header (this is the unfiltered baseline). Confirm the `canonical_smiles`
   column reports the semantic type **Molecule** — the substructure search command is only
   available on a Molecule column; if the semantic type is missing, stop and report which
   column was found and what type it carried.

## Scenarios

### Scenario 1: Query with matches, then query with no matches

Tests the fundamental top-menu substructure search path: the menu command prepares an empty
substructure filter and opens its sketcher; a matching query then narrows the table, and a
non-matching one empties it.

Steps (numbering continues from the Setup above, which is Step 1):

2. In the top menu, navigate to **Chem | Search | Substructure Search**. Verify that a
   substructure filter card for the molecule column appears in the Filter Panel with an
   **empty** query, that its sketcher opens for input, and that the row count is still the
   Setup baseline — the command prepares the search, it does not perform one.
3. Enter the probe SMILES `c1ccccc1` (benzene ring) into that filter card's sketcher and
   commit it. Verify the table filters to a non-empty, proper subset of the original rows:
   - The filtered row count is strictly greater than zero.
   - The filtered row count is strictly less than the total row count noted in Setup.
   - Every row that survives the filter contains the benzene ring, and no row that was
     filtered out contains it.
4. Without closing the sketcher, replace the query in it with one that matches no row in
   the dataset — use `[Au]` (a gold atom, absent from every molecule in smiles.csv).
   - Verify the table shows a row count of zero rather than reverting to the full
     unfiltered table.
   - Verify that no molecule in the column actually contains the query fragment, so the
     empty result is the right answer and not a search that failed silently.
   - Verify the total row count is unchanged from the Setup baseline — no rows were deleted,
     the table is merely filtered to an empty set.

Expected:
- Step 3: the table header reports a narrowed row count, strictly between zero and the full
  row count; every row that passes the filter genuinely contains the benzene ring substructure.
- Step 4: zero rows match, and no molecule in the dataset contains the query fragment; the
  total dataset row count is unchanged.

### Scenario 2: Re-invoke replaces, not stacks

Tests that calling Chem | Search | Substructure Search again replaces the current filter
rather than intersecting with it.

Steps:

1. Starting from the state at the end of Scenario 1 (table filtered to zero rows), in the
   top menu choose **Chem | Search | Substructure Search** once more. This time a
   **Substructure search** dialog appears first, asking which molecule column to search —
   the searches already run have cached a hidden helper column of canonical SMILES that is
   itself typed Molecule, so the command no longer sees a single candidate column. The
   picker offers only the user-visible `canonical_smiles`; confirm it with **OK**. Verify
   the filter is prepared empty again — the table returns to the full Setup row count and
   the sketcher opens.
2. Enter the probe SMILES `C(=O)O` (carboxylic acid) into the filter card's sketcher and
   commit it. Verify the same card was reused — the Filter Panel still holds exactly one substructure
   filter for the molecule column — and that the table now shows a non-empty set of rows:
   - The filtered row count is strictly greater than zero and strictly less than the total.
   - The count differs from the benzene-ring count produced in Scenario 1 Step 3, showing
     the new query was applied rather than the old one retained.
   - Verify the result is not equal to the full table row count — the filter is active.

Expected:
- Step 2: the row count is strictly between zero and the full row count; it differs from the
  benzene-ring Scenario 1 count; every surviving row contains the carboxylic acid fragment
  and no filtered-out row contains it.

### Throughout both scenarios

The run logs no console error of its own. The dev home page's PowerPack widget failures are
ambient and are excluded by name; any other console error is a failure of this path.

## Automation notes

- **Top-menu path:** the substructure search command is registered via `//top-menu: Chem |
  Search | Substructure Search` in `public/packages/Chem/src/package.ts` around L769. Drive
  it through the actual top-menu DOM hierarchy, not via `grok.functions.call` — the menu
  gesture is what this path tests.
- **Sketcher dialog:** the dialog opened by this command contains a SMILES input field. Enter
  the probe by typing into that field (`fill` on the SMILES input) and confirming with **OK**.
  Do not attempt to draw atoms by clicking on the canvas — there is no stable DOM target for
  individual canvas pixels inside the sketcher widget.
- **Readiness barrier:** after clicking OK, poll `df.filter.trueCount` (through
  `grok.shell.t.filter.trueCount`) until it settles rather than reading it in the same turn
  as the OK click. The filter computation is asynchronous and a snapshot taken immediately
  after OK may still reflect the previous state.
- **Empty result — there is no DOM channel for it.** The grid renders no empty-state element
  and no "0 rows" text node when a filter excludes everything; `.d4-grid-empty-message` does
  not exist in the product (live recon on dev, recorded in
  `references/decisions/chem.yaml:1136`). Do not write a locator for it — it would never match
  and never fail. What tells a correctly empty filter apart from one that failed silently is
  the data: `df.filter.trueCount === 0` with `df.rowCount` unchanged, plus a row-level
  re-derivation showing that no molecule contains the query fragment.
- **Row-level substructure check:** `grok.functions.call('Chem:getRdKitModule')` returns the
  RDKit module in-page; `get_qmol(query)` plus `get_mol(smiles).get_substruct_match(pattern)`
  per row re-derives the predicate independently of the filter. Over the 1000 rows of
  smiles.csv it takes ~300 ms, so it is affordable on every query rather than only as a
  sample. Delete each mol after use, and count unparseable molecules and assert that count is
  zero — otherwise the row-level result silently covers fewer rows than it claims.
- **Probe selection:** each probe must match some but not all molecules, or the
  strictly-less-than bound is satisfied by an inactive filter. On smiles.csv `c1ccccc1`
  matches 924 of 1000 and `C(=O)O` matches 314; `[Au]` matches none. If two probes that
  should differ produce the same count, treat it as a probe-selection failure, not a product
  pass.
- **Why the second invocation goes through a column picker:** running a search caches helper
  columns on the dataframe (`saveColumn`, `Chem/src/chem-searches.ts:181`), one of which,
  `~canonical_smiles.canonicalSmiles`, holds canonical SMILES and is itself detected as
  Molecule. `columns.bySemTypeAll(MOLECULE)` therefore returns two columns from then on, and
  `searchSubstructureEditor` (`Chem/src/package.ts:744`) switches from calling the function
  directly to showing `ui.dialog({title: 'Substructure search'})`, which the platform names
  `[name="dialog-Substructure-search"]`. An untitled dialog — the sketcher — is named
  `dialog-`, so the two are told apart by name. Verified on dev with a standalone probe:
  cancelling the first invocation without entering a query caches nothing, and the second and
  third invocations then open the sketcher directly.
- **Filter Panel positive control:** nothing here opens the Filter Panel explicitly — the
  command's own `getFiltersGroup()` call materialises it, and the card count observed right
  after the first invocation is 1. Any assertion comparing card counts must first assert the
  count is non-zero, or an unrendered panel satisfies it with 0 === 0.
- **Scenario 2 "differs from" assertion:** pin the Scenario 1 benzene count to a local
  variable before entering Scenario 2, then assert `trueCount !== pinnedCount`.
- **Console errors:** the dev home page logs PowerPack widget failures
  (`PowerPack:MostRecentEntities`, `PowerPack:RecentlySharedWithMe`, both from a broken
  `datagrok_reader` DB password) unrelated to this path, each followed by a separately logged
  translated stack trace keyed by the `ID = <id>` its message carries. Exclude those by name
  and by id, then require zero remaining — a filter narrow enough to match only the word
  "substructure" passes no matter what else the run logged.
- **Probe numbers are logged before their assertions** (`[probe] ...` lines): the Molecule-column
  list at Setup and again before the re-invocation, the sketcher/picker counts, the filter card
  state, and each query's kept/excluded/unparseable partition. On a green run Playwright prints
  no assertion messages, so without the logs the run records that the checks held but not the
  values they held at — and the Molecule-column count is the number that explains why the two
  invocations behave differently.
- Teardown: close the test table in a `finally` block, even on failure.
