---
slug: spectronaut-malformed-mislabel
phase: 12-spectronaut-input-coverage
status: resolved
goal: find_root_cause_only
created: 2026-05-15T18:30:09Z
source: 12-UAT.md gap "Spectronaut precursor import reports only genuinely malformed lines as 'malformed'" (test 1, also observed test 2)
---

# Debug: Spectronaut streaming import mislabels correctly-filtered rows as "malformed"

## Symptoms (pre-filled from UAT)

- **Truth (expected):** The Spectronaut precursor import surfaces only genuinely
  malformed/unparseable lines as "malformed"; rows dropped by the by-design
  filters (decoy/contaminant, q-value > threshold) are normal filtering, not errors.
- **Actual:** Importing the committed fixture (`files/demo/spectronaut-hye-precursor.tsv`)
  showed `grok.shell.info: "Spectronaut import: skipped 30 malformed line(s)"`,
  and importing `files/demo/spectronaut-hye-mix.tsv` showed `"skipped 159 malformed
  line(s)"` — yet every data line in both files is structurally well-formed.
- **Errors:** None (info-level message; import otherwise succeeds with correct output).
- **Reproduction:** Test 1 and Test 2 in 12-UAT.md.
- **Timeline:** Discovered during Phase 12 UAT (2026-05-15).

## Root Cause (established, with empirical corroboration)

`parseSpectronautStream` in `src/parsers/spectronaut-parser.ts` collapses **four
semantically distinct** line outcomes into a single `skipped` counter, then reports
that counter to the user as **"malformed line(s)"**.

`handleFields` (`src/parsers/spectronaut-parser.ts:307-344`) returns `false` for:

1. `f.length < expectedFields` (line 308) — **genuinely malformed** (truncated row).
2. `!protein` (line 311) — empty protein group (degenerate; treat as malformed/unusable).
3. `protein.startsWith('CON__') || protein.startsWith('REV__')` (line 312) —
   **by-design decoy/contaminant filter** (duckdb-parity, intentional).
4. `q !== null && q > qValueThreshold` (line 316) — **by-design q-value filter**
   (duckdb-parity, intentional).
5. `qty === null && q === null` (line 320) — **genuinely unusable** row (both casts null).

The caller (lines 392-393, 417-418) does `if (handleFields(...)) rowCount++; else
skipped++;` — every `false`, regardless of reason, increments the **same** `skipped`.
Line 426 then reports: `grok.shell.info(\`Spectronaut import: skipped ${skipped}
malformed line(s)\`)`. Categories 3 and 4 are correct, intentional filtering — not
"malformed".

### Empirical corroboration (two independent files, exact filter-count matches)

- **Test 1 — `files/demo/spectronaut-hye-precursor.tsv`:** `wc -l` = 493 (492 data
  lines); field-count distribution = `492 lines x 11 fields` (header is 11 cols) →
  **zero structurally malformed lines**. Reported "30 malformed" = the fixture's
  decoy + q>0.01 + degenerate filter-branch rows (the fixture is purpose-built to
  exercise every R2 branch).
- **Test 2 — `files/demo/spectronaut-hye-mix.tsv`:** rows with numeric
  `EG.Qvalue > 0.01` (col 16) = **exactly 159**; CON__/REV__ rows = 0. Reported
  "159 malformed" == **exactly** the 159 q-filtered rows. Decisive: the "malformed"
  count is the by-design q-filter drop count, verbatim.

### Impact

Data integrity is **not** affected — aggregation output is provably correct
(matches the committed duckdb golden via Plan 03's `streaming output equals duckdb
golden` test; Test 2's result is numerically equivalent to the pre-Phase-12 text
path because PG.IBAQ is group-constant across all 741 groups). The defect is the
**user-facing message**: on the 2.6 GB headline file the streaming path drops the
vast majority of input rows by design (millions of precursor rows collapse to
protein groups; q>0.01 and decoys filtered), so the user would see
*"skipped ~millions of malformed line(s)"* on a perfectly correct import — a false
corruption signal on the exact use case Phase 12 exists to deliver.

## Files Involved

- `src/parsers/spectronaut-parser.ts`
  - `handleFields` (307-344): returns bare `boolean`, conflating malformed vs
    by-design-filtered.
  - caller loops (392-393, 417-418): single `skipped++` for all `false`.
  - progress message (402-405): `"${skipped} malformed skipped"`.
  - completion message (425-426): `grok.shell.info("... skipped N malformed line(s)")`.
- `src/parsers/spectronaut-parser.ts` `pivotSpectronaut` (text path, 47-80) applies
  the **same** q-filter + decoy drop **silently** — the correct UX precedent to match
  for the by-design category.

## Suggested Fix Direction (for gsd-planner --gaps)

Make `handleFields` return a discriminated outcome instead of a bare boolean
(e.g. `'kept' | 'malformed' | 'filtered'`):
- `'malformed'` ← `f.length < expectedFields`; `!protein`; `qty === null && q === null`.
- `'filtered'` ← `CON__/REV__` decoy/contaminant; `q > qValueThreshold`.
- `'kept'` ← aggregated.

Track two counters (`malformed`, `filtered`). User-facing surface:
- Report only `malformed` as "malformed/unparseable line(s)" (keep current wording
  for that genuine category).
- For `filtered`: either stay **silent** (matches `pivotSpectronaut`/text-path UX,
  least surprising) or, if a count is desired, phrase it neutrally as
  "filtered N low-confidence / decoy rows" — never "malformed".
- Update the in-progress message (402-405) symmetrically so it does not say
  "malformed" for filtered rows.

No change to filtering logic, aggregation, or output shape — message/categorization
only. Existing Spectronaut tests must stay green; add a test asserting that a
fixture containing only decoy + q>threshold rows yields `malformed == 0`.

## Resolution

Root cause found and fully evidenced. Handed off to gsd-planner (gap_closure) for
a targeted fix plan. No code changed in this debug session (goal:
find_root_cause_only).
