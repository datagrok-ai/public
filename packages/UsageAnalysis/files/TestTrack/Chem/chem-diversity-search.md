---
feature: chem
realizes_atlas:
  - chem.cp.diversity-search-viewer
realizes:
  - chem.cp.diversity-search-viewer
priority: p1
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: smoke
realized_as:
  - chem-diversity-search-spec.ts
related_bugs: []
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-08-20-chem-new-12
    timestamp: "2026-08-20T00:00:00Z"
    failure_keys: []
    review_round: 1
  e:
    verdict: PASS
    cycle_id: direct-gate-e-2026-08-22-chem-diversity-search
    timestamp: 2026-08-22T12:33:45Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: direct-gate-b-2026-08-22-chem-diversity-search
    timestamp: "2026-08-22T11:44:00Z"
    spec_runs:
      - spec: chem-diversity-search-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 21
        failure_keys: []
        run_mode: headless-cold
expected_results:
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      The Chem Diversity Search viewer is present in the table view and renders
      exactly one molecule card per selected structure: the card count equals
      the viewer's current limit and equals the length of its diverse id set.
      The viewer header names the metric and fingerprint then in force
      (Tanimoto, Morgan).
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      After the Distance Metric is changed to Cosine, the viewer re-runs and
      arrives at a different diverse id SET — compared order-independently, so a
      re-ordering of the same molecules does not satisfy it — while still
      filling the current limit. The header names the newly applied pairing
      (Cosine, Morgan), and the re-run logs no matching console error.
  - anchor: "Scenario 1 Step 6"
    expectation: >-
      After changing the Number of Structures (limit) to a value different from
      the initial value, the count of molecule cards shown in the viewer matches
      the new limit value.
  - anchor: "Scenario 2 Step 3"
    expectation: >-
      After the Fingerprint is changed from Morgan to MACCS, the viewer re-runs
      and arrives at a different diverse id SET — compared order-independently —
      while still filling the current limit. The header names the newly applied
      fingerprint (MACCS), and the re-run logs no matching console error.
  - anchor: "Scenario 2 Step 4"
    expectation: >-
      Setting Size to normal gives EVERY molecule card a 200x100 canvas, and
      setting it to large then gives every card a 300x150 canvas — measured
      across the whole card set, not a sample, so a viewer that resized one tile
      and left the rest fails. The card count survives the resize, and large is
      strictly wider and taller than normal. The small tier is not exercised.
  - anchor: "Scenario 2 Step 5"
    expectation: >-
      After Row Source is set to Filtered, the diverse id set holds EXACTLY
      min(filtered row count, limit) ids — not "at most", which a viewer showing
      nothing would also satisfy — every displayed id belongs to the filtered
      row set, and the viewer still renders cards.
---

# Chem — Diversity Search viewer

## Setup

1. Open the demo dataset at `System:DemoFiles/chem/smiles.csv` via the top menu
   File | Open or the Browse sidebar under Demo Files > chem.
2. Verify that the molecule column (canonical_smiles) is recognised as semType Molecule
   (cells render as 2D structures, not raw text).

## Scenarios

### Scenario 1: Add viewer and verify metric and limit controls change the displayed set

Steps:
1. With `smiles.csv` open, navigate to the top menu Chem | Search | Diversity Search.
2. Confirm that a Chem Diversity Search viewer appears in the table view layout.
3. Note the row ids of the molecule cards currently shown in the viewer (the initial
   diverse subset; record at least the first three ids visible in the viewer panel).
4. Verify that the viewer renders exactly one molecule card per selected structure — the
   card count equals the viewer limit — and that the header names the active metric and
   fingerprint (Tanimoto, Morgan).
5. Open the viewer's property panel (gear icon or right-click > Properties on the
   viewer header). Change the Distance Metric from its current value to a different
   metric: from Tanimoto to Cosine. Wait for the viewer to re-compute. Verify that the SET
   of molecule ids shown in the viewer differs from the set recorded in step 3 — a different
   ordering of the same molecules does not count.
6. In the property panel, change the Number of Structures (limit) from the default
   value (12) to 6. Wait for the viewer to re-render. Verify that the viewer now shows
   exactly 6 molecule cards.

Expected:
- The Chem Diversity Search viewer is added to the table view on step 2.
- The viewer renders exactly one card per selected structure on step 4; the console check
  belongs to step 5, where the re-run happens.
- Changing the Distance Metric causes the displayed molecule ids to change (step 5).
- Changing the limit to 6 results in exactly 6 molecule cards displayed (step 6).

### Scenario 2: Fingerprint, size, and row-source controls

Steps:
1. With the Chem Diversity Search viewer from Scenario 1 still open, open its property
   panel (gear icon or right-click > Properties on the viewer header).
2. Record the current set of molecule ids shown in the viewer.
3. Change the Fingerprint type from Morgan to MACCS. Wait for the viewer to finish
   computing. Verify that the SET of molecule ids changes relative to the set recorded in
   step 2 — a different ordering of the same molecules does not count.
4. Change the Size property from Normal to Large. Verify that EVERY molecule card grows in
   displayed dimensions, not merely the first one.
5. Apply a column filter to the table (for example, use the Filter Panel to narrow the
   row set to fewer than 50 rows). In the viewer property panel set Row Source to
   Filtered. Verify that the viewer re-renders showing only structures from the
   filtered row set, and that the count of displayed cards is exactly the filtered row
   count clamped to the current limit.

Expected:
- Changing the Fingerprint type produces a different diverse subset (step 3).
- Changing Size visually resizes the molecule cards (step 4).
- Setting Row Source to Filtered limits the displayed molecules to the filtered rows
  (step 5).

## Automation notes

- **Actuation channel.** All five controls (Distance Metric, Number of Structures, Fingerprint,
  Size, Row Source) are driven through the viewer's property panel, opened with the title-bar gear.
  The Misc category is collapsed when the panel opens, so it is expanded before any of its rows is
  touched; the enum rows are driven by clicking the value cell and selecting in the `select` the
  click creates, and Limit through `input.property-grid-slider-textbox` (its value cell is
  zero-height and not clickable). The table is narrowed through the Filter Panel: the card for a
  low-cardinality column is switched to categorical and its in-card search box is typed into, with
  the card's "filter by search results" box asserted on first.
- **Size map:** normal tiles are 200x100 and large tiles 300x150, read from `canvas.style.width` /
  `style.height`. The canvas `width`/`height` ATTRIBUTES are multiplied by `window.devicePixelRatio`
  (`render-molecule.ts:23-25`), so asserting them only holds at ratio 1. Measure across
  `.chem-diversity-search > div` children, not `children[0]` — the claim is about every card. The
  small tier is not exercised by this spec.
- **Set comparison, not sequence.** `renderMolIds` order is not the claim; the claim is which
  molecules were selected. Compare sorted copies, or a metric change that merely re-sorts the
  same subset — the exact regression worth catching — passes.
- **Header text** is the readable channel for the active pairing: `.chem-similarity-header`
  contains "Tanimoto, Morgan" initially, then "Cosine, Morgan", then "MACCS". Scope it to
  `[name="viewer-Chem-Diversity-Search"]` — the class is shared with Chem Similarity Search.
- **Probe numbers are logged before their assertions** (`[probe] ...` lines). On a green run
  Playwright prints no assertion messages, so without the logs the run records that the checks
  held but not the values they held at.
- **Error channel:** the viewer reports failures through `grok.shell.error(...)`
  (`chem-diversity-viewer.ts:53-55`) and `malformedDataWarning` through a balloon, so the "no error"
  checks read `.d4-balloon.error`, not the console. Balloons auto-hide after 5 s, so absence is read
  instantaneously right after the re-run settles rather than waited on.
