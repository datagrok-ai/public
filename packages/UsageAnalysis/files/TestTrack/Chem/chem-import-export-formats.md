---
feature: chem
realizes_atlas:
  - chem.cp.import-export-formats
realizes: []
priority: p0
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: smoke
related_bugs: []
scope_reductions:
  - id: SR-03
    check: E-TRACE-02
    rationale: |
      Fixture substitution, benign and recorded for traceability: the scenario
      names mol1K.sdf, which is absent from the server, so the SDF import steps
      run against System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf instead. The
      assertions (record count from "$$$$" terminators, Molecule semType,
      rendered cells) are fixture-independent, so the substitution does not
      weaken them.
    verdict_status: SCOPE_REDUCTION
  - id: SR-04
    check: E-TRACE-02
    rationale: |
      Fixture substitution for the V3000 check, recorded for traceability. Step 5
      names ApprovedDrugs2015.sdf as the V3000 sample, but SR-03 already spends that
      file on steps 1-4, and no other V3000 SD file ships with the platform. The step
      therefore builds one: it reads the V3000 molfile at
      System:DemoFiles/chem/mol/v3000.mol, appends the "$$$$" record terminator that
      makes it an SD file, opens the result, and deletes it. Opening the bare molfile
      instead — what the spec did before 2026-08-20 — measured the fixture rather than
      V3000 support: the SD parser ends a record on "$$$$" regardless of format
      version, so the V2000 aspirin.mol yields 0 records with its terminator stripped
      and v3000.mol yields 1 with one added (openchemlib 7.5.0, checked directly).
      This substitution strengthens the step rather than weakening it; what the
      scenario asserts about ApprovedDrugs2015.sdf's own record version is untested
      and not claimed here.
    verdict_status: SCOPE_REDUCTION
  - id: SR-05
    check: E-LAYER-ALIGN-01
    rationale: |
      Actuation reduced to the registered function, on evidence. The scenario
      originally named a UI route — right-click the InChI column header, then
      Convert > InChI to Molecule — and the spec called Chem:InchiToMol directly
      instead. Reconnaissance on dev (2026-08-20) established that the UI route
      does not exist: the column Actions panel on a string `inchi` column offers
      nine actions (string, default, Add filter, Remove, Rename..., Extract
      numbers, Extract, Change Type..., Hide) and no Convert entry, while a
      molecule column does show Convert Notation... in the same place; the
      header context menu was expanded in full, including its Add submenu, and
      the cell context menu was checked separately; running
      grok.data.detectSemanticTypes leaves the column at semType null and the
      menu unchanged. The function itself is registered (role: converter,
      inputRegexp (InChI\=.+)), so the scenario's premise about the converter
      holds — only the menu path is absent. The step text and the Scenario 3
      Step 2 anchor were lowered to the function route accordingly, so the
      scenario no longer claims a gesture the product does not offer.
      Uncovered as a consequence: the cell-level parse-failure marker, and only
      that — it could only be observed from a rendered column produced through
      the UI. The semType clause of the Scenario 3 Step 2 expectation is not
      affected, because it needs no UI: the step reads the semantic type off the
      output parameter of the converter's own call
      (call.outputParams['smiles'].property.semType), which the platform builds
      from the script header declaration #output: string smiles { semType:
      Molecule } (Chem/scripts/inchiToMol.py:7) rather than from anything the
      spec writes.
    verdict_status: SCOPE_REDUCTION
  - id: SR-06
    check: E-ASSERT-STRENGTH-01
    rationale: |
      Declared coverage gap on Scenario 1 Step 5: the step makes no
      structure-versus-raw-text claim about the drawn V3000 cell. The
      colour-spread discriminator the SDF and .mol steps use counts a pixel as
      structure ink when its channel spread reaches 30, which detects RDKit's
      heteroatom palette. The pinned fixture
      System:DemoFiles/chem/mol/v3000.mol is a pure hydrocarbon (M V30 COUNTS
      19 19, all nineteen atoms carbon) and Chem maps carbon and hydrogen to
      [0,0,0] (Chem/src/utils/chem-common-rdkit.ts:34-37), so a correctly drawn
      structure scores zero: the reading was a constant, and the Gate B FAIL of
      2026-08-22 (0/3, "Expected: > 0, Received: 0") was that constant, not a
      product signal. The assertion was removed rather than re-thresholded. The
      substitutes considered were rejected on evidence, not preference: swapping
      in a heteroatom-bearing molecule would contradict the fixture this
      scenario names in Step 5; a spatial rule (ink extent measured against this
      grid's own text control, or a two-state contrast against a forced text
      fallback) carries no magic threshold but rests on measurements that could
      not be taken, and the built SD file's likely control column carries the
      molfile's empty title line, so it may hold no ink to compare against. What
      remains uncovered in Step 5 is "the V3000 cell paints a structure rather
      than raw text".

      Corrected 2026-08-22, first revision of this record. It closed by saying
      V3000 parsing itself stayed covered by the conversion assertions, "which
      fail if the molblock is not parsed". That sentence is refuted by the
      product source and is withdrawn: grok.chem.convert does not throw on a
      parse failure. _convertMolNotation pre-seeds its result with the sentinel
      'MALFORMED_INPUT_VALUE' and returns it
      (Chem/src/utils/convert-notation-utils.ts:60-62), and js-api forwards
      whatever the function call produced (js-api/src/chem.ts:917-923). With
      V3000 parsing wholly broken the sentinel satisfied all four assertions the
      step made: no throw, a string, non-empty, no "V30" in it. So while SR-06
      stood, Step 5 had no working V3000 signal at all — the ink claim was gone
      and the conversion claim was a sentinel pass-through. This is repaired
      rather than declared: the step now also asserts that the conversion result
      parses as a structure, via grok.chem.checkSmiles, which routes through
      Chem:isSmiles (RDKit get_mol, Chem/src/package.ts:2289) — a different
      platform function from Chem:convertMolNotation, so the sentinel cannot
      satisfy both. V3000 parsing is covered again, by that guarded assertion
      and not by the four that preceded it.

      Third claim that lapsed silently with this record, now declared
      separately: Scenario 1 Step 3's per-cell rendering claim. The
      colour-spread reading is an existential over a bounded strip, not a
      per-cell check, so it does not support "every cell renders as a
      structure". That narrowing is recorded in SR-07 rather than here, because
      it concerns a different step and a different fixture; SR-06 names it only
      so this record's accounting is complete.
    verdict_status: SCOPE_REDUCTION
  - id: SR-07
    check: E-EXPECT-COVERAGE-01
    rationale: |
      Declared narrowing on Scenario 1 Step 3 (and the same reading in Step 6).
      The expectation states that every cell of the imported Molecule column
      renders as a structure with no raw-text fallback cells. What the check
      actually establishes is weaker, in two ways that are now written down.

      Covered: on the molecule column's own screen rectangle, below the column
      header and capped at 400 px of height, at least one pixel carries RDKit's
      heteroatom colour, while the monochrome control column of the same grid
      measured in the same pass scores exactly zero. That discriminates "the
      molecule renderer drew structure ink somewhere in the visible strip" from
      "the column is monochrome raw text throughout", with no pixel *count*
      threshold involved on either side — the molecule side asserts "> 0" and
      the control side "= 0". There is one threshold, but it is per pixel
      rather than per count: a pixel scores as coloured when its channel spread
      reaches 30 of 255, the margin the spec's own comment states sits far
      above antialiasing noise and far below real heteroatom ink.

      Not covered: the per-cell claim. ApprovedDrugs2015.sdf is a
      thousands-row table and the strip spans only the first screenful of rows,
      so cells below it are never read; and within the strip the assertion is
      existential, so the state where most cells fall back to raw text while a
      single heteroatom-bearing cell renders passes. The obvious strengthening —
      require coloured ink inside every visible cell's own row band — was
      rejected on the evidence in SR-06, not on preference: carbon and hydrogen
      map to [0,0,0] (chem-common-rdkit.ts:34-37), so any all-carbon structure
      in the fixture scores zero coloured pixels while rendering perfectly, and
      a per-cell rule would red on the product being correct. A per-cell rule
      that does not depend on the palette (per-cell ink extent against a
      same-row text control) rests on measurements this run could not take.
      Nothing here is claimed about cells outside the measured strip.
    verdict_status: SCOPE_REDUCTION
  - id: SR-08
    check: E-EXPECT-COVERAGE-01
    rationale: |
      Declared narrowing on Scenario 1 Step 7, the .smi step, and on that step
      only — it is deliberately not folded into SR-06 or SR-07, which cover
      Steps 5, 3 and 6 and never named Step 7. The step's text asks for three
      things: exactly three rows, the Molecule semantic type, and "every cell
      renders a structure". The first two are asserted directly. The third had
      no observation behind it at all until 2026-08-22: the step counted
      non-empty values in the molecule column, which are the three SMILES the
      step itself had just written to the fixture, so the count could not fail
      on any rendering regression.

      Covered now: the Molecule cell renderer is bound to the imported column
      (read from the current view's dataframe, and guarded to be the table this
      step opened); the grid of that table painted non-blank pixels, with the
      floor being the paint barrier's own exit condition rather than a tuned
      number; and no console.error call was made during this step, the buffer
      having been cleared at the step's start so the claim covers this step and
      not an earlier one.

      Not covered: the structure-versus-raw-text distinction, and the per-cell
      part of the claim. Neither is a preference. The colour-spread ink
      discriminator the SDF and .mol steps use cannot be pointed at this table,
      on two facts checked in the product source rather than assumed. First,
      _importSmi prepends a single "SMILES" header line and parses the result
      as CSV (Chem/src/file-importers/smi-importer.ts:4-12), so the imported
      table has exactly one column; the helper's control search skips columns
      whose semType is Molecule, finds nothing else, and falls through to the
      row-number-gutter control, which the recon note at the head of the spec
      records as never having been exercised on any measured table — so a
      reading taken there is neither known good nor known bad. Second, the
      fixture's middle line c1ccccc1 is benzene, every atom carbon, and Chem
      maps carbon and hydrogen to [0,0,0] (chem-common-rdkit.ts:34-37); a
      per-cell rule over this fixture would therefore score zero on a
      correctly drawn cell, the same constant-discriminator trap SR-06
      documents for the V3000 fixture. Changing the fixture to heteroatom-only
      SMILES would contradict the three SMILES the scenario names in Step 7.
      What remains uncovered in Step 7 is "each of the three cells paints a
      structure rather than raw text"; what is covered is that the structure
      renderer is bound, that the grid painted, and that the renderer raised
      no error while doing so.
    verdict_status: SCOPE_REDUCTION
  - id: SR-09
    check: E-EXPECT-COVERAGE-01
    rationale: |
      Declared narrowing on Scenario 2 Step 4, the round-trip. The expectation
      states that the re-imported Molecule column holds the same structures as
      the original rows. What is checked is three rows of a thousand — first,
      middle and last — so the reading is existential standing for a universal,
      the same shape SR-07 records for the per-cell rendering claim. The
      scenario body authorises exactly this ("the canonical SMILES of at least
      three spot-check rows (first, middle, last) must match"), so this record
      is a declaration rather than a repair; it exists because the narrowing was
      nowhere written down. Not covered: the 997 rows between the spot checks,
      where a structure could be silently altered or dropped by the export or
      the re-import without any assertion noticing. The row-count equality
      asserted alongside it does bound one failure mode — rows lost outright —
      but says nothing about identity per row. Strengthening to all thousand
      rows was not attempted here: it would mean a thousand round trips through
      Chem:convertMolNotation inside one step, inside a 300 s test budget that
      the preceding export step (Scenario 2 Step 3) can already spend up to
      120 s of on its download wait.
    verdict_status: SCOPE_REDUCTION
  - id: SR-11
    check: E-EXPECT-COVERAGE-01
    rationale: |
      Declared narrowing on the console clause, at all five read points that
      carry it — the SDF, V3000, .mol and .smi imports and the InChI
      conversion. It is recorded here rather than left to the exclusion
      subsection in the Automation notes because that subsection describes the
      filter, while this records what the assertion is thereby narrowed to.
      Nothing about the filter changes with this record.

      Covered: at each read point, no console.error call was made during that
      step other than entries the named-exclusion rule drops — an entry naming
      the known home-widget function PowerPack:RecentlySharedWithMe, an entry
      naming the datagrok_reader credential together with the authentication
      failure, or a logger stack-trace continuation belonging to an
      already-excluded event. The code computes counted = matched + other and
      captured = counted + excluded, and every read point asserts counted = 0
      (chem-import-export-formats-spec.ts:353-357), so what the assertion
      consumes is the captured total MINUS the named exclusions, not the raw
      captured total. The excluded entries are not invisible: each one is
      counted and printed in full in the run log beside the pre-exclusion
      total, so growth in the exclusion is readable from the run.

      Not covered: a chem-side console.error whose own text happens to name
      PowerPack:RecentlySharedWithMe, or to name datagrok_reader together with
      an authentication failure, is dropped with the ambient noise and does not
      fail the assertion. A chem-side error that merely follows an excluded one
      is unaffected — its head line is not a trace continuation, so it is kept
      and counted. Also not covered, and by design: the stand's datagrok_reader
      authentication failure itself, which this spec is silent about
      deliberately and which is flagged to the operator separately.

      This record exists because the narrowing was stated in the exclusion
      subsection and in the spec, but three outward carriers still stated the
      pre-exclusion form: the Scenario 1 Step 3 anchor, the Step 3 body prose,
      and the SR-10 retirement paragraph. All three now name the exclusion
      limit. The exclusion rule itself is correct as written and is unchanged.
    verdict_status: SCOPE_REDUCTION
gate_verdicts:
  b:
    verdict: PASS
    cycle_id: direct-gate-b-2026-08-22-chem-import-export-formats-r9
    timestamp: 2026-08-22T15:05:00Z
    spec_runs:
      - spec: chem-import-export-formats-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 55
        run_mode: headless-cold
        failure_keys: []
  a:
    verdict: PASS
    cycle_id: 2026-08-20-chem-new-09
    timestamp: "2026-08-20T12:00:00Z"
    failure_keys: []
    review_round: 2
  e:
    verdict: PASS
    cycle_id: direct-gate-e-2026-08-22-chem-import-export-formats-r7
    timestamp: 2026-08-22T15:20:36Z
    failure_keys: []
expected_results:
  - anchor: "Scenario 1 Step 3"
    expectation: >-
      The opened table's Molecule column carries semType=Molecule and every cell
      renders a structure thumbnail — no raw-text fallback cells, and no console
      error is raised while rendering. The console clause is checked over the
      console.error channel only: the capture hook patches console.error, so
      window.onerror, unhandled promise rejections and console.warn are outside
      what "no console error" is asserted to mean here. Narrowed, see SR-11:
      what is asserted is the captured total minus the named ambient exclusions
      — the stand's datagrok_reader authentication failure and the known
      home-widget function — not the raw captured total. Narrowed, see SR-07:
      what is checked is that structure ink appears somewhere in the molecule
      column's first screenful of rows while the monochrome control column of
      the same grid scores zero. Not covered: the per-cell part of the claim —
      cells below the measured strip are never read, and within it a single
      rendered cell satisfies the check even if the rest fall back to raw text.
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      Row count shown in the table header matches the file's record count for
      each format: mol → 1 record, mol2 → 1 record, smi → row count equals the
      number of SMILES lines in the file, sdf → row count equals the number of
      records delimited by "$$$$" in the file.
  - anchor: "Scenario 2 Step 3"
    expectation: >-
      The export writes an .sdf artefact without raising an error or a warning
      notification, and the record count of that artefact (its number of "$$$$"
      terminators) equals the number of rows exported.
  - anchor: "Scenario 2 Step 4"
    expectation: >-
      Re-importing the artefact the export wrote produces a table whose row
      count equals the source row count and whose Molecule column holds the same
      structures as the original rows — the round-trip is lossless at the
      molecule-identity level (canonical SMILES match before and after).
      Narrowed, see SR-09: the identity match is checked on three spot-check
      rows (first, middle, last) of a thousand, so rows between them are not
      covered.
  - anchor: "Scenario 3 Step 2"
    expectation: >-
      The converter's own output parameter is typed semType=Molecule — read from
      the call's output metadata, which the platform builds from the script
      header, not from any column the check writes (SR-05) — and its output for
      a well-formed InChI string is a parseable structure differing from the raw
      InChI input, and its failure on a malformed InChI surfaces as the
      conversion failing rather than as a silently valid structure — the reverse
      direction of To InChI is exercised. Not covered: the cell-level
      parse-failure marker (a crossed-out structure or a red cell indicator).
      The step observes the converter's own outcome, not the rendered cell. The
      marker stays uncovered: the reconnaissance recorded in SR-05 found no UI
      route to this converter, so there is no rendered-cell step to observe it
      from.
realized_as:
  - chem-import-export-formats-spec.ts
---

# Chem — Import/Export Formats

## Setup

1. Open the Datagrok platform and navigate to Browse > Files > System:DemoFiles > chem.
2. Confirm the following demo files are accessible: `smiles.csv`,
   `sdf/ApprovedDrugs2015.sdf`, `mol/aspirin.mol`, `mol/v3000.mol`, and
   `System:DemoFiles/bio/ngl-formats/adrenalin.mol2`. (These are the sources for
   the format-specific sub-steps below.) `mol1K.sdf`, which earlier revisions of
   this scenario named, is **not** on the server — see SR-03; the SDF steps run
   against `ApprovedDrugs2015.sdf` instead.

Actuation channel: files are opened by handing the path to the platform's own
file-handler dispatch (the `OpenFile` function), not by navigating Browse and
double-clicking. Only the Browse navigation gesture is skipped — the importers
that are the subject of this scenario (`importSdf`, `importMol`, `importMol2`,
`importSmi`) are the ones the dispatch selects and runs, exactly as they would
be for a user's double-click.

The export side is different, and the gap there is real: the SDF export is
invoked by calling the `Chem:saveAsSdf` function directly, so the route a user
takes to it — the table hamburger menu / File > Export > "As SDF…" — is not
exercised. Everything from the export dialog onwards (the dialog itself, its OK,
and the file it writes) is covered; only that menu gesture is not.

## Scenarios

### Scenario 1: Import — four file formats recognised as Molecule semType

Steps:
1. Confirm that `mol1K.sdf`, the SDF this scenario used to name, is genuinely
   missing from the server. This is a note about the fixture, not a claim about
   the product: it supports SR-03 and nothing else, and it keeps the
   substitution visible in the run rather than only in the frontmatter.
2. Open `System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf` (SDF format) — the
   substitute SR-03 records.
3. Verify that the opened table's molecule column carries semType `Molecule`,
   that the `Molecule` cell renderer is bound to it, that the grid canvas is
   painted rather than blank, and that no `console.error` call is made while
   rendering — the captured total minus the named ambient exclusions is
   asserted (SR-11 records what that exclusion does and does not leave
   covered), and the rdkit renderer's
   own error patterns only classify what survived exclusion. Only
   `console.error` is hooked, so `window.onerror`, unhandled promise rejections
   and `console.warn` fall outside the claim; the capture hook is proved live
   before its silence counts (see Automation notes). (The check reads the column model; the
   column-header badge itself is not inspected.) The structure-versus-raw-text part is checked over the
   molecule column's first screenful of rows and is satisfied by ink anywhere in
   that strip, not cell by cell — SR-07 records what that does and does not
   cover.
4. Verify that the row count of the opened table equals the number of records in
   the file — for an SDF, the number of `$$$$` terminators the file contains.
   The count is read from the file rather than assumed, so no fixture-specific
   number is written down here.
5. Close this table. Check V3000 support. `ApprovedDrugs2015.sdf` is already
   spent on steps 1–4 and no other V3000 SD file ships with the platform, so
   build one: take the V3000 molfile `System:DemoFiles/chem/mol/v3000.mol`,
   append the `$$$$` record terminator that makes it an SD file, open the
   result, then delete it (SR-04 records the substitution).

   Verify that it imports as one row, that its column is typed `Molecule`, that
   the cell still holds the V3000 molblock, and — the one assertion here that
   depends on V3000 actually being parsed — that converting that cell yields a
   SMILES which itself parses as a structure (`grok.chem.checkSmiles`). That
   parseability check is the load-bearing one and must not be dropped as
   redundant: `grok.chem.convert` does not throw when the molblock fails to
   parse, it returns the sentinel `'MALFORMED_INPUT_VALUE'`
   (`Chem/src/utils/convert-notation-utils.ts:60-62`, forwarded unchanged by
   `js-api/src/chem.ts:917-923`), so the four checks around it — no throw, a
   string, non-empty, no `V30` in it — are satisfied in full with the molblock
   never parsed. The row count and the
   semType are kept only as fixture-health checks: the record ends where this
   step wrote `$$$$`, and `_importSdfString` stamps the semType without ever
   parsing the molblock, so neither of them can fail on a V3000 regression.
   The grid must paint and no `console.error` call must be made. Whether the
   V3000 cell draws a structure rather than falling back to raw text is **not**
   checked here — see SR-06.
6. Close this table. Open the single-structure `.mol` file
   `System:DemoFiles/chem/mol/aspirin.mol` (V2000, and it carries its `$$$$`
   terminator). Verify: one-row table, molecule column typed `Molecule`, cells
   painted, no `console.error` calls.
7. Close this table. No `.smi` file ships with the platform, so create one for
   this check: write three SMILES, one per line (`CCO`, `c1ccccc1`, `CC(=O)O`),
   to a temporary `.smi` file in a writable share, then open it the same way a
   user would. Verify: the table has exactly three rows, the molecule column
   carries the Molecule semantic type, and every cell renders a structure.
   Of the rendering part, what is checked is that the Molecule cell renderer is
   bound to the column, that the grid painted rather than staying blank, and
   that no `console.error` call is made — not that each of the three
   cells draws a structure rather than raw text; SR-08 records what that does
   and does not cover, and why the colour-spread reading used in Steps 3 and 6
   cannot be pointed at this one-column, all-carbon fixture.
   Delete the temporary file afterwards — this check must leave nothing behind.

   The file must go through the `.smi` handler, not the CSV parser: `.smi` has
   its own importer (`Chem/src/package.ts:2055`, `ext: 'smi'`), so asserting on
   `smiles.csv` instead would exercise a different code path and miss any
   regression in `importSmi`.
8. Close this table. Open the single-molecule `.mol2` file
   `System:DemoFiles/bio/ngl-formats/adrenalin.mol2`, routed through its own
   handler (`importMol2`, `Chem/src/package.ts:2069`, `ext: 'mol2'`). Verify:
   the table has one row per `@<TRIPOS>MOLECULE` block in the file — the
   expected count is read from the file, not assumed — and the molecule column
   is typed `Molecule`.

### Scenario 2: Export as SDF and round-trip re-import

Steps:
1. Open `smiles.csv` from System:DemoFiles/chem. The table loads with a
   Molecule column (semType=Molecule).
2. Invoke the SDF exporter (a user reaches it from the table hamburger menu or
   File > Export > "As SDF…"; that menu gesture is the one gap declared in
   Setup). The SDF export dialog appears with a Molecules column selector.
3. Accept the defaults and click OK. Verify that the dialog closes, that the
   export writes out an `.sdf` file, and that no error or warning notification
   appears. Verify that the written file's record count — its number of `$$$$`
   terminators — equals the number of rows exported.

   The warning channel is proved live before it is used as evidence of absence:
   raise a known warning, watch its balloon appear, wait for it to clear, and
   only then assert that the export raises none.
4. Re-import the file the export wrote. Verify that the re-imported table has
   the same row count as the source table and that its molecule column holds the
   same structures — the canonical SMILES of at least three spot-check rows
   (first, middle, last) must match the corresponding original rows.

### Scenario 3: InChI-to-molecule converter (reverse direction of To InChI)

Steps:
1. Prepare a small CSV file with a column of InChI strings (for example,
   copy the InChI values produced by Chem > Calculate > To InChI on
   `smiles.csv`, then open that result as a new table, or use the
   `System:DemoFiles/chem/smiles.csv` table and run To InChI first to
   produce the InChI column). The table must have a column whose cell
   values begin with `InChI=`.
2. Run the registered `Chem:InchiToMol` converter over the InChI column and add
   its output to the table as a new column, typed with the semantic type the
   converter itself reports for that output. The converter is registered
   with `meta.role=converter` and `inputRegexp` `(InChI\=.+)`, but on dev
   (2026-08-20) it is not reachable from the column-header menu, the cell menu,
   or the column Actions panel — see SR-05. Invoking the registered function is
   therefore the only available route, and it is what this step covers.
3. Verify: the converter's output is typed `Molecule` — read from the semantic
   type the platform attaches to the output parameter of the converter's own
   call, which comes from the script's `#output: string smiles { semType:
   Molecule }` declaration and not from anything this check writes. Verify also
   that the output for each well-formed InChI row is a parseable structure that
   differs from the raw InChI string it came from, and that a malformed InChI
   does not silently produce a valid structure — the conversion fails instead.

   Not covered: the cell-level parse-failure marker (a crossed-out structure or
   a red cell indicator). What is checked is the converter's own outcome on
   malformed input, which is a different thing from what the cell then draws.
   The marker stays uncovered: the reconnaissance recorded in SR-05 found no
   UI route to this converter, so there is no rendered-cell step to observe
   it from.

## Automation notes

Out of scope by operator decision (2026-08-20): **opening a `.mol` file that carries
no `$$$$` record terminator.** Two handlers claim the extension — `importSdf`
(`ext: 'sdf,mol'`, Chem/src/package.ts:2041) and `importMol` (`ext: 'mol'`, :2082) —
and the SDF one wins, so such a file opens as an empty table with no warning. That is
not under test here and no ticket is being raised for it. Do not add a check for it,
and do not treat the empty result as a regression if it surfaces.

This is a limit on the routing case only. V3000 support itself IS under test: Scenario 1
Step 5 builds a genuine SD file from the V3000 molfile and asserts that the molblock
converts to a SMILES that parses as a structure, which is the signal that depends on
V3000 being parsed. Non-emptiness and the absence of `V30` are not that signal and
must not be taken as sufficient on their own: `grok.chem.convert` returns the sentinel
`'MALFORMED_INPUT_VALUE'` rather than throwing when parsing fails
(`Chem/src/utils/convert-notation-utils.ts:60-62`), and that sentinel is itself a
non-empty, `V30`-free string.
The `.mol` format is likewise still covered, by `aspirin.mol` in Step 6, which carries
its terminator.

### Console-error assertion: SR-10 retired on measurement (2026-08-22)

The console clause used to assert a *filtered* count: captured `console.error` calls were
reduced to five rdkit-renderer patterns (`rdkit[-_]cell[-_]renderer`, "method not found",
`'gS'`, `cellRenderer\.render`, `NullError`) and only those were required to be zero, so
anything outside the patterns was invisible. That narrowing was recorded as SR-10 and was
explicitly conditional: the reason for not asserting the whole buffer was that the stand's
ambient `console.error` baseline had never been measured, and an exclusion list written from
expectation rather than observation would be the same defect mirrored. Every read point was
instrumented to log the captured total, the matched count and the first unmatched messages so
one clean run would produce that baseline.

The baseline was measured. Gate B (cycle `direct-gate-b-2026-08-22-chem-import-export-formats-r6`)
ran three attempts; all five read points logged `captured=0 rdkit-matched=0 other=0 []` on
every attempt. The capture hook was armed six times per attempt and its planted probe was
verified caught each time, so the empty buffer is an observation, not an unarmed silence.
An exclusion list derived from that run would have had zero entries, so there is nothing to
exclude and no reason to keep the filter on the assertion. SR-10 is therefore retired rather
than dropped: the assertion stopped being filtered by the rdkit pattern set, which survives only
as diagnostic classification in the log line. As written that day the assertion consumed the
whole captured count; it no longer does. The named-exclusion subsection below, added later the
same day on an observation, subtracts the named ambient entries before the assertion, so what is
consumed is **captured minus named-excluded** — see SR-11 for what the clause covers now. The
Scenario 1 Step 3 anchor is restored to the general claim it originally made, less that
exclusion.

What the baseline does and does not license, carried forward so the widened claim is not
overstated:

1. **Channel.** Only `console.error` is hooked. `window.onerror`, unhandled promise
   rejections and `console.warn` never enter the buffer, so every zero-assertion here means
   "no `console.error` call", not "no error". This limit is written into the Step 3 anchor and
   into the spec's step prose.
2. **Windows.** The claim covers only the six windows `beginErrorCapture` arms. Login and
   platform boot happen before the first arming; Scenario 1 Step 8 (`.mol2`) arms nothing and
   claims nothing; Scenario 2 Step 3 (export) uses the balloon channel, not this one. Those
   stretches are unmeasured.
3. **Sample.** One stand, three consecutive runs, one browser profile. A zero baseline on
   another stand, another profile, or a noisier day is not established by this measurement; a
   future non-zero reading is a reason to name the noise, not to re-narrow the assertion.

#### First named exclusion: the stand's `datagrok_reader` DB-auth noise (2026-08-22)

The non-zero reading point 3 anticipated arrived immediately, and was handled the way point 3
prescribes. What is excluded is one stand-side fault — the server-side database connection made
with the `datagrok_reader` credential failing authentication — surfacing through whichever
home widget happens to be retrying while a capture window is open. It is **not** a rule about
one named widget, and two different widgets have now raised it.

Gate B cycle `direct-gate-b-2026-08-22-chem-import-export-formats-r7` ran three attempts. On
**attempt 1 only**, and at the **InChI read point only**, the buffer held
`captured=3 rdkit-matched=0 other=3`. The three entries are one platform event, not three:
the `PowerPack:RecentlySharedWithMe` home widget failed its server-side database connection
("FATAL: password authentication failed for user `datagrok_reader`", with the
`org.postgresql...ConnectionFactoryImpl.doAuthentication` stack), the platform surfaced it
through its logger into `console.error`, and the logger's own stack-trace continuation
("Translating stack trace... ID = pGCYd", then "Stack trace pGCYd") followed.

The next cycle, `direct-gate-b-2026-08-22-chem-import-export-formats-r8`, saw the same fault
again — again on **attempt 1** and again at the **InChI read point** — this time under **two**
widget names: `PowerPack:RecentlySharedWithMe` and `PowerPack:MostRecentEntities`, together
reading `captured=4 counted=0 named-excluded=4`. The second widget is nowhere in the rule's
text; it was dropped by the credential clause ("password authentication failed for user" together
with `datagrok_reader`), which is the rule behaving exactly as written and is why the rule is
described here by the credential rather than by the widget. In both cycles nothing in the chem
import/export path is involved, the event landed inside that step's capture window only because
the window happened to be open when a widget retried, it did not recur on attempts 2 or 3, and
the other four read points measured zero on every attempt. The grid readings were byte-identical
across attempts, so the chem render path itself raised nothing on any run.

The exclusion is narrow. An entry is excluded when it names a known widget function
(`PowerPack:RecentlySharedWithMe`) or names `datagrok_reader` together with the authentication
failure — the conjunction is required, so a chem-side database or renderer error is not
swallowed with it. The logger's stack-trace continuation carries no identity of its own and its
id is minted per occurrence, so it cannot be matched by a rule naming an id: a continuation is
excluded only when it belongs to an already-excluded event, either by directly following one in
the capture buffer or by carrying an id that an already-excluded continuation announced. A stack
trace that follows a **kept** error is kept and counted with it.

Excluded messages are counted and printed, not made to vanish. The log line at every read point
prints the full captured total before exclusion, the number excluded by name, and **every**
excluded message — uncapped, since the growth case is precisely the one a capped list would hide
(in r8 a three-entry slice would have concealed the fourth entry, the one that revealed the
second widget). The assertion consumes `captured − named-excluded`. So an exclusion that grows,
or a new message the list does not name, is visible in the run log rather than silently dropped.

Still open, and deliberately not settled here: whether a `datagrok_reader` password-authentication
failure on dev is something this spec should be silent about at all. It surfaced in two
consecutive cycles, and in the second under two widget names in a single run, so it is recurring
rather than a one-off. **This spec is silent about it by design, not because it is resolved.**
The exclusion makes the chem assertion stop reporting a stand-credential problem it has no
standing to adjudicate; it does not make the problem acceptable, and it is flagged to the
operator separately.
