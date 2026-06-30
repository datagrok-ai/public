---
created: 2026-05-15T19:54:42.755Z
title: Tighten streaming↔text parity for empty-protein Spectronaut rows
area: import
files:
  - src/parsers/spectronaut-parser.ts:337 (handleFields: !protein → 'malformed')
  - src/parsers/spectronaut-parser.ts:60-61 (pivotSpectronaut: !protein → silent continue)
  - src/tests/spectronaut-parser.ts (regression net — no empty-protein categorization test)
---

## Problem

Code review of Phase 12 gap-closure 12-04 raised WR-01 (and WR-03 as a narrower
instance). After 12-04, `parseSpectronautStream`'s `handleFields` returns a
discriminated `'kept' | 'malformed' | 'filtered'` outcome so by-design-filtered
rows (CON__/REV__ decoys, numeric q > threshold) are silently dropped, matching
the text path `pivotSpectronaut`.

However the **empty-`PG.ProteinGroups`** branch was deliberately left in the
user-surfaced `'malformed'` category (`spectronaut-parser.ts:337`,
`if (!protein) return 'malformed'`), whereas the text path drops the
byte-identical row **silently** (`pivotSpectronaut:61`, `if (!protein) continue;`).
A correct Spectronaut precursor report containing unassigned-precursor rows
(blank protein cell) would therefore fire a "skipped N malformed line(s)"
message on the streaming path but not on the text path.

This is NOT a defect against Phase 12 / plan 12-04: 12-04's `must_haves.truths`
line 1 explicitly lists "empty protein" as a genuine-malformed case, and the
`LineOutcome` docblock + task behavior spec match the code. The Phase 12
verifier adjudicated it as a non-blocking follow-up. It is a *parity-tightening*
improvement — full streaming↔text silence parity — not a correctness fix.
WR-03 (`both numeric casts null → 'malformed'`) is the same class of
text-path-parity question and also currently untested.

## Solution

TBD — decide the intended contract first (this is a deliberate behavior choice,
not a mechanical bug):

- Option A: reclassify empty-protein (and re-examine both-casts-null) from
  `'malformed'` to `'filtered'` in `handleFields` so the streaming path is
  byte-silent exactly where `pivotSpectronaut` is, then update the
  `LineOutcome` docblock and add a regression test pinning empty-protein →
  silent (the one by-design drop the 12-04 net misses).
- Option B: keep empty-protein as `'malformed'` (treating it as a genuinely
  unusable row — no protein to attribute quant to) but make the text path
  consistent, OR explicitly document the intentional asymmetry.

Either way, add the missing regression test so the empty-protein and
both-casts-null categorization is pinned in whichever direction is chosen.
Whichever path: keep aggregation output byte-identical (the committed duckdb
golden / per-cell text-path equivalence tests must stay green — predicates and
the `agg` map must not change, only the outcome categorization/messages).
Reference: Phase 12 `12-REVIEW.md` WR-01/WR-03, `12-VERIFICATION.md` WR-01
adjudication section.
