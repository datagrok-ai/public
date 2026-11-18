@footguns @batch-2 @serial
Feature: Leaks and speed: what repeated use and pathological input must never do (Crux mode)
  Memory and speed under repeated use and pathological input (reactions read, drawn and freed in rounds, a grid
  scrolled, a 1,000-residue reaction, a 60-step route, a Transform of 50,000 molecules): these scenarios say what
  the user must be able to count on. Every threshold comes from engine-level numbers (measured in node on Chem's
  vendored chem-wasm) times a stated margin; where a stand's numbers differ, the threshold is set again from that
  measurement.
  - Engine memory: a scheme read, drawn and freed 300 times a round, 6 rounds, for a template, a route, a specific
    reaction and an unread text: the wasm memory grows 2.36 → 4.52 MB in round 1 and not at all after. Mapping 300 a
    round: flat at 16.5 MB. So any growth that keeps going after a warm-up is a leak. Threshold: after 5 warm-up rounds,
    at most 1 MiB more over the next 50 (16 wasm pages: allocator slack; today's engine-level growth is 0).
  - Live objects: back to their count before the rounds, exactly (no margin: a count is not noisy).
  - Cells: a scheme read 0.16 ms (p95 0.49), drawn 0.80 ms (p95 2.0, SVG; the canvas draw is not slower); about 15 new
    rows a scrolled page: 15 × 2.5 ms = 38 ms at p95, plus the grid's own work. Threshold: no frame over 100 ms (2.6×),
    and no worse than twice the same grid scrolled over a molecule column.
  - A 1,000-residue reaction (20,033 characters) draws in 7.0 s today, the tab frozen; 300 residues 0.84 s. Threshold:
    no main-thread task over 1 s on its account, and the cell shows the picture or words saying it is too large.
  - A 60-step route 157 ms, 20 steps 59 ms: threshold 500 ms (3×) for the 60-step one.
  - Transform: 0.093 ms a molecule on the engine (4.7 s for 50,000); it yields every 25 ms (reactions.ts:114). Threshold:
    no main-thread task over 200 ms while it runs, done within 60 s with its counts (10×: the canonicalisation it does
    first runs on Chem's workers and is not in the engine number).

  Background:
    # instrumentation, before the login's navigation: DevTools could not list the page's memories
    Given the page records the WebAssembly memories and Crux sketchers it makes
    And user is logged in
    And the package autostarts have completed
    And the Chemical Engine is "Crux"

  @catalogue-PERF-01
  Scenario: Opening and closing the molecule cell editor 50 times leaves Crux where it started
    Given user opens a table "molecules" with a SMILES column "molecule" holding "c1ccc2c(c1)cc(C(=O)NCC)n2"
    And user has opened and closed the cell editor of "cell 1 of molecule" 5 times
    And user notes the live Crux sketcher objects and the wasm memory of every Crux instance
    When user opens and closes the cell editor of "cell 1 of molecule" 50 times, drawing a bond each time
    Then the live Crux sketcher objects are as many as noted
    And the wasm memory of every Crux instance grew by at most 1 MiB
    And the page holds no more detached sketcher roots than before

  @catalogue-PERF-01
  Scenario: Switching the current reaction cell 50 times leaves the Reaction panel's objects where they started
    Given user opens the reaction corpus "files/tests/reactions-crux.csv" of Chem
    And user has made rows 1 to 10 current in turn
    And user notes the live reaction schemes, the Reaction panel's canvases and their size listeners
    When user makes rows 1 to 10 current in turn, 5 times
    # exact, so a leak of one scheme a switch fails (the rows' schemes were all made in the Given)
    Then the live reaction schemes are as many as noted
    And the Reaction panel's canvases and size listeners are as many as noted
    And the wasm memory of the main thread's Crux grew by at most 1 MiB

  @catalogue-PERF-02
  Scenario: A 2,000-row reaction column scrolls without a frame over 100 ms
    Given user opens a table of 2,000 reactions and routes made from the reaction corpus
    And the reaction column of the grid is 600 px wide and its rows 150 px tall
    When user scrolls the grid from the first row to the last, one page at a time
    Then no frame took more than 100 ms
    And the worst frame is at most twice the worst frame of the same grid scrolled over a molecule column
    And the live reaction schemes are at most the cell cache's bound

  @catalogue-PERF-04
  Scenario: A 60-step route draws within half a second
    Given user opens a table "route" with a reaction column "reaction" holding a 60-step chain "CO>>CCO>>CCCO>>…"
    When the grid draws "cell 1 of reaction"
    Then the cell's first draw took at most 500 ms

  # the mapping worker's wasm memory is out of the page's reach; it is checked below the page, at engine level, by
  # crux-core's memcheck, its reaction-map section (the same chem-wasm with RXNMapper's weights, mapping repeated,
  # memory flat after warm-up; run on every change to crux-core and every re-vendor) and by a PERF-05 sweep (16.5 MB
  # flat over 1,800 mappings). The two worker-memory steps are gone; the page checks the answers' time.
  @catalogue-PERF-05
  Scenario: Mapping 20 reactions in the Reaction panel answers each within 2 s
    Given user opens the reaction corpus "files/tests/reactions-crux.csv" of Chem
    And user has mapped 5 reactions in the Reaction panel
    When user maps 20 more reactions in the Reaction panel, one after another
    Then each mapping answered within 2 s

  @catalogue-PERF-06
  Scenario: Transform on 50,000 molecules keeps the tab responsive and reports its counts
    Given user opens a table of 50,000 molecules "smiles_200K" cut to 50,000 rows
    When user runs the transformation "Aromatic Bromination" on the molecule column
    # the products column's creation is the platform's (a 50,000-value string column's category
    # dictionary built and sorted, ~270 ms), judged on its own
    Then no main-thread task took more than 200 ms while it ran, the products column's creation aside
    And the platform created the products column in one task of at most 500 ms
    And it finished within 60 s with a message counting products, errors and rows with no reaction
