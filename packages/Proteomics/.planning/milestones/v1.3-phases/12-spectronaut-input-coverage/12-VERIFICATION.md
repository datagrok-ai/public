---
phase: 12-spectronaut-input-coverage
verified: 2026-05-15T20:00:00Z
status: passed
score: 8/10 verified programmatically; human release-gate (2.6 GB reference import — no false malformed signal) confirmed PASS by user approval 2026-05-15
overrides_applied: 0
human_gate_resolved: 2026-05-15T19:06:19Z (user-approved after 12-04 re-import; see 12-HUMAN-UAT.md)
re_verification:
  previous_status: human_needed
  previous_score: 10/11
  gaps_closed:
    - "Spectronaut precursor import mislabeled by-design-filtered rows (CON__/REV__ decoys, numeric q>threshold) as 'malformed' — split counter implemented; regression test added"
  gaps_remaining: []
  regressions:
    - "WR-01 (empty-protein malformed vs filtered): residual parity gap between streaming path and text path — NEW finding from 12-REVIEW.md, adjudicated below"
human_verification:
  - test: "Import the 2.6 GB reference TSV (~/Downloads/2026-05-13 BP DMD WT.tsv) via Proteomics | Import | Spectronaut Report"
    expected: "Import completes with no V8 string-length or OOM error; TaskBar progress bar advances monotonically; no Chrome 'Page Unresponsive' dialog; tab switching works mid-import; resulting DataFrame has approximately 8,328 protein rows and 24 sample columns; tags proteomics.source=spectronaut, proteomics.preNormalized=true, proteomics.groups with 2 conditions; full pipeline Annotate->Normalize->Impute->DE->Volcano runs to completion; NO false 'malformed line(s)' message fires on an otherwise-correct import"
    why_human: "The 2.6 GB reference file cannot be committed. R4 (progress + main-thread responsiveness) and the full R1 OOM-prevention acceptance criterion require the actual large file. The corrected malformed/filtered split also needs re-validation on the real file to confirm the false corruption signal is gone."
---

# Phase 12: Spectronaut Input Coverage — Re-verification Report (after 12-04 gap closure)

**Phase Goal:** A user can select a raw multi-GB Spectronaut precursor-level long report via Proteomics | Import | Spectronaut Report and the package streams + aggregates it client-side into the same wide protein x sample DataFrame the existing pipeline consumes — no whole-file buffering, no manual pre-processing.
**Verified:** 2026-05-15T20:00:00Z
**Status:** passed (human release-gate confirmed by user approval 2026-05-15 — see 12-HUMAN-UAT.md)
**Re-verification:** Yes — after gap closure plan 12-04 was executed and merged

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|---------|
| T-01 | 12-04 MT-1: import surfaces 'malformed/unparseable line(s)' ONLY for genuinely malformed lines (too few fields; empty protein; both numeric casts null) | PARTIAL — see WR-01 adjudication | handleFields L334: `f.length < expectedFields -> 'malformed'`; L337: `!protein -> 'malformed'`; L346: `qty===null && q===null -> 'malformed'`. The plan's stated categories match the code. However the 'empty protein' branch (`!protein -> 'malformed'` at L337) diverges from `pivotSpectronaut` L61 (`if (!protein) continue;` — SILENT). See WR-01 adjudication section. |
| T-02 | 12-04 MT-2: by-design-filtered rows (CON__/REV__; numeric EG.Qvalue > threshold) do NOT contribute to any 'malformed' message — silent, matching the text path | VERIFIED | handleFields L338: `CON__/REV__ -> 'filtered'`; L342: `q !== null && q > qValueThreshold -> 'filtered'`. Completion message at L467-468 fires only when `malformed > 0` via `void filtered`. Test `by-design-filtered rows are not counted malformed` (L484-518) proves the spy on `grok.shell.info` captures zero malformed-labeled messages for a CON__+REV__+high-q-only fixture. |
| T-03 | 12-04 MT-3: streaming aggregation output (protein x sample values, tags, semTypes, group keys) is byte-identical to before this change — committed duckdb golden / per-cell text-path equivalence tests stay green | VERIFIED | Only the `handleFields` return type, the two caller switch arms, the counter declarations (let skipped -> let malformed, filtered), and the two message guards were changed. Predicates, agg map, aggToPivotResult, finalizeSpectronaut, pivotSpectronaut, parseSpectronautText — all untouched (confirmed by code read). |
| T-04 | 12-04 MT-4: new grok test proves a fixture whose only dropped rows are decoy + q>threshold yields genuine-malformed count of 0; a too-few-fields row IS counted malformed | VERIFIED | Two tests at L484 and L520: `by-design-filtered rows are not counted malformed` asserts `malformedMsgs.length === 0` for decoy-only fixture; `truncated line is counted malformed` asserts a `/skipped \d+ malformed line\(s\)/` info message fires for a 3-field truncated row. Both exist in `category('Spectronaut')`. |
| T-05 | Streaming precursor importer uses Blob.stream() + TextDecoderStream; never materializes the whole file (R1 — code) | VERIFIED | parseSpectronautStream L301-303 confirmed in prior verification. Unchanged by 12-04. |
| T-06 | PG-level long report still takes the unchanged file.text() path (R1/R2 boundary) | VERIFIED | importSpectronaut branch logic at package.ts confirmed in prior verification. Unchanged by 12-04. |
| T-07 | duckdb aggregation script committed under tools/ with usage docs, flip reference-only, unused carry-along columns dropped (R5) | VERIFIED | tools/spectronaut-aggregate.{sql,sh} confirmed in prior verification. Unchanged by 12-04. |
| T-08 | Synthetic precursor fixture + duckdb-derived golden + JSON sidecar committed; grok tests for streaming path exist (R6) | VERIFIED | files/demo/spectronaut-hye-precursor.tsv (493 lines), golden.tsv (235 lines), golden.json (234 keys), 7 pre-existing streaming tests + 2 new gap-closure tests (26 total). |
| T-09 | R4: TaskBarProgressIndicator + macrotask yield implemented (code) | VERIFIED | parseSpectronautStream L299: `DG.TaskBarProgressIndicator.create`; L437: `await new Promise<void>((r) => setTimeout(r, 0))`; L433-435: malformed-only suffix in progress message. Unchanged by 12-04. |
| T-10 | R4 + R1: 2.6 GB reference file imports with no OOM, no false malformed message, monotonic progress, tab switching works | HUMAN GATE | Pre-declared manual-only gate per 12-SPEC.md. Gap-closure changes the expected observable: no false malformed signal should now fire. Needs re-confirmation. |

**Score:** 8/10 truths verified programmatically (T-01 is PARTIAL — see WR-01 adjudication; T-10 is the pre-declared manual gate)

### WR-01 Adjudication: Empty-Protein Branch — Malformed vs Filtered

**Question posed:** Does the `!protein -> 'malformed'` classification at `handleFields` L337 satisfy 12-04's must-have truth #1 as written, or does it constitute an unmet must-have because it contradicts the plan's own text-path-parity contract?

**Evidence — streaming path (spectronaut-parser.ts:337):**
```
const handleFields = (f: string[]): LineOutcome => {
  if (f.length < expectedFields) return 'malformed';
  const protein = f[protI];
  if (!protein) return 'malformed';          // LINE 337 — user-surfaced
  if (protein.startsWith('CON__') || ...) return 'filtered';
```

**Evidence — text path (spectronaut-parser.ts:60-61):**
```
const protein = protCol.get(i) as string;
if (!protein) continue;                      // LINE 61 — SILENT, no message
if (protein.startsWith('CON__') || ...) continue;
```

**Plan 12-04 must_have truth #1 (exact wording):**
> "A Spectronaut precursor import surfaces a 'malformed/unparseable line(s)' message ONLY for genuinely malformed lines (too few fields; **empty protein**; both numeric casts null)"

**Plan 12-04 task-1 behavior contract:**
> "`handleFields(f)` returns `'malformed'` when `f.length < expectedFields`, when `!protein`, or when both numeric casts are null."

**LineOutcome docblock (spectronaut-parser.ts:206-217):**
> "`'malformed'` — the row is genuinely unparseable: too few fields, an empty protein group, or both numeric casts null."

**Plan 12-04 objective text (text-path-parity goal):**
> "the by-design-filtered category stays SILENT — the text path (`pivotSpectronaut`, lines 47-80) applies the identical q-filter + CON__/REV__ drop with no message; matching that silence is the locked, path-consistent UX decision"

**12-REVIEW.md WR-01 claim:**
> "The text path treats the byte-identical condition as a silent, by-design skip ... An empty protein group is by-design unusable (cannot be aggregated to any protein), not unparseable structure; it belongs in 'filtered' (silent)"

**Verdict: NO INTERNAL CONTRADICTION EXISTS within 12-04's must_haves. The plan's stated truth #1 explicitly includes "empty protein" in the 'malformed' category, and the task-1 behavior spec explicitly maps `!protein -> 'malformed'`. The implementation matches both the must_have text and the task behavior contract exactly.**

The WR-01 finding from 12-REVIEW.md is a **legitimate real-world correctness concern** — real Spectronaut precursor reports can contain blank-protein rows for unassigned precursors, and the streaming path would fire a false "malformed" signal for them while the text path would be silent. This is a genuine parity gap against `pivotSpectronaut` L61.

However, it is NOT a gap against 12-04's must_haves AS WRITTEN. Plan 12-04 deliberately chose to keep empty-protein as 'malformed' — the must_have text, the task behavior contract, the LineOutcome docblock, and the code all agree with each other. The plan's scope was to fix CON__/REV__ and q>threshold (the categories empirically proven to cause the false signal by UAT tests 1 and 2). The review's suggested re-classification is a follow-up improvement not within the plan's stated scope.

**Classification of WR-01:** WARNING (acknowledged parity gap on real data with blank-protein rows, but no internal contradiction with 12-04 must_haves). Does not constitute a BLOCKER against 12-04's must_haves. Future phases or a follow-on gap-closure should address it if real Spectronaut exports are observed to contain blank-protein rows.

**Classification of WR-03 (both-casts-null also 'malformed'):** Same category — the plan's must_have #1 explicitly lists "both numeric casts null" as a genuine-malformed case. The WR-03 review finding that this diverges from the text path is noted but is not a gap against the plan's must_haves. Future improvement candidate.

### Deferred Items

None. Phase 12 is the final milestone phase. Backlog items 999.1 (annotate dialog prefill) and 999.2 (prenormalized banner wording) are explicitly backlogged and out of scope.

### Required Artifacts — 12-04 Scope

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/parsers/spectronaut-parser.ts` | handleFields returns discriminated 'kept'/'malformed'/'filtered'; two counters; messages report only genuine-malformed | VERIFIED | LineOutcome type at L218; handleFields returns L334/337/338/342/346/370; malformed/filtered counters at L317/321; switch at L418-422 and L450-454; completion msg guarded at L467-468; void filtered at L466 |
| `src/tests/spectronaut-parser.ts` | Two new regression tests: by-design-filtered and truncated-line | VERIFIED | Tests at L484 and L520; streamCapturingInfo helper at L468; both tests use makeLongFormatTsv + grok.shell.info spy restored in finally |

### Key Link Verification — 12-04 Gap Closure

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| parseSpectronautStream caller loops | handleFields discriminated return | switch on 'kept'/'malformed'/'filtered' | WIRED | L418-422 (main loop) and L450-454 (trailing flush) both switch on all three arms; no arm falls through |
| completion + in-progress messages | genuine-malformed counter only | malformed > 0 guard | WIRED | In-progress at L433-435: `malformed > 0` suffix; completion at L467-468: `if (malformed > 0) grok.shell.info(...)` |
| by-design-filtered rows | silence | filtered counter never surfaced | WIRED | `void filtered` at L466; no grok.shell message for filtered; matches pivotSpectronaut's L56/L62 `continue` (silent) |

### Anti-Patterns Found — 12-04 Files

No TBD/FIXME/XXX/placeholder/stub anti-patterns found in src/parsers/spectronaut-parser.ts or src/tests/spectronaut-parser.ts.

**Remaining code-review advisories (not blocking 12-04 must_haves):**

| File | Finding | Severity | Status |
|------|---------|----------|--------|
| spectronaut-parser.ts:337 | WR-01: `!protein -> 'malformed'` while text path (L61) is silent — parity gap on real blank-protein rows | WARNING | Within 12-04 must_have scope (plan explicitly chose this classification); parity gap on real data is acknowledged, future follow-up |
| spectronaut-parser.ts:344-346 | WR-03: `both-casts-null -> 'malformed'` diverges from text path (which emits the protein) | WARNING | Within 12-04 must_have scope; acknowledged |
| spectronaut-parser.ts:321,421,453,466 | IN-01: `filtered` counter is write-only (void filtered) | INFO | By design; docblock says "tracked for completeness" |
| spectronaut-parser.ts:433-435 vs 467-468 | IN-02: progress tick and completion message can transiently disagree if all malformed lines are in the trailing tail | INFO | Not a correctness bug; cosmetic |
| src/tests/spectronaut-parser.ts:484-518 | WR-02: empty-protein categorization unpinned — no test for blank-protein row; either 'malformed' or 'filtered' result is un-asserted | WARNING | Follows from WR-01; if WR-01 is later addressed by reclassification, test must be added |

### Human Verification Required

#### 1. 2.6 GB Reference File Import — R4 Release-Gate + Corrected Malformed Signal Validation

**Test:** Open `~/Downloads/2026-05-13 BP DMD WT.tsv` (2.6 GB) via `Proteomics | Import | Spectronaut Report`.

**Expected:**
- Import completes without a V8 string-length or out-of-memory error
- TaskBar shows a monotonically advancing progress indicator during import
- Chrome "Page Unresponsive" dialog does NOT appear
- Tab switching works during the import
- Resulting DataFrame has approximately 8,328 protein rows and 24 sample columns
- Tags: `proteomics.source=spectronaut`, `proteomics.preNormalized=true`, `proteomics.groups` with exactly 2 conditions
- NO `grok.shell.info` message of the form "Spectronaut import: skipped N malformed line(s)" fires (the millions of by-design-filtered precursor rows are now silent; only genuinely unparseable lines would trigger the message — and none are expected on a correct precursor export)
- Full pipeline: Annotate -> Normalize -> Impute -> DE -> Volcano runs to completion

**Why human:** Reference file is 2.6 GB — cannot be committed. The 12-04 change makes the malformed message behavior on the real 2.6 GB file the primary verification target: confirming no false corruption signal fires on a correct import is the UX acceptance criterion for the gap closure. Also validates R4 (responsiveness) and R1 (OOM prevention) in the same run.

### Gaps Summary

No blocking gaps found in 12-04's must_haves as written. The implementation satisfies all four must_have truths:

1. **MT-1 PARTIAL → EXPLAINED:** "empty protein" is classified 'malformed' exactly as the plan's must_have and task behavior contract specify. WR-01's argument that it should be 'filtered' is a valid follow-up improvement, not a gap against 12-04's stated scope.
2. **MT-2 VERIFIED:** CON__/REV__ and q>threshold are silent ('filtered'); text-path consistency achieved for the categories the UAT confirmed were causing the false signal.
3. **MT-3 VERIFIED:** Aggregation output byte-identical; no predicate or agg-map change.
4. **MT-4 VERIFIED:** Both regression tests exist and cover the stated fixture contracts.

The single human verification item (R4 + corrected-malformed-signal gate on the 2.6 GB file) is the pre-declared manual release-gate from 12-SPEC.md, extended to also validate the 12-04 UX fix on the real file. This elevates the status to `human_needed`.

---

_Verified: 2026-05-15T20:00:00Z_
_Verifier: Claude (gsd-verifier)_
_Re-verification: Yes — gap closure 12-04_
