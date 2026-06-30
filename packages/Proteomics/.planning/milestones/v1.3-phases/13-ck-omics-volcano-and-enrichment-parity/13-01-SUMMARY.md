---
phase: 13-ck-omics-volcano-and-enrichment-parity
plan: 01
subsystem: testing
tags: [semantic-type, detectors, spectronaut, userDataStorage, scaffolding]

requires:
  - phase: 12-spectronaut-input-coverage
    provides: spectronaut-candidates-parser, spectronaut-parser autoPopulateGroups
provides:
  - SEMTYPE.SUBCELLULAR_LOCATION + mirrored detectors.js detector (R1 prerequisite)
  - 13-WAVE0-FINDINGS.md resolving A2 (Candidates sign), A3 (D-09 site), A4 (cache)
  - registered SubcellularLocation and Volcano test categories
affects: [13-03, 13-04, 13-05, 13-06]

tech-stack:
  added: []
  patterns:
    - "Duplicated-literal contract: SEMTYPE const value mirrored verbatim in detectors.js"
    - "Wave-0 de-risking findings doc consumed by later-wave plans"

key-files:
  created:
    - .planning/phases/13-ck-omics-volcano-and-enrichment-parity/13-WAVE0-FINDINGS.md
    - src/tests/subcellular-location.ts
    - src/tests/volcano.ts
  modified:
    - src/utils/proteomics-types.ts
    - detectors.js
    - src/package-test.ts

key-decisions:
  - "Candidates AVG Log2 Ratio = log2(Numerator/Denominator); no hidden inversion — R3 flip is conditional/per-row only"
  - "D-09 fix is direction-only at differential-expression.ts:289-292 (DE dialog default), not the parser"
  - "Subcellular cache = per-key userDataStorage map (not one blob) with __schema_v invalidation"

patterns-established:
  - "SEMTYPE literal duplicated by contract into detectors.js, grep-gated"
  - "Wave-0 FINDINGS.md de-risks downstream waves before they plan their edits"

requirements-completed: [R1, R3]

duration: 67min
completed: 2026-05-17
---

# Phase 13 Plan 01: Wave-0 De-risking + Scaffolding Summary

**Added SEMTYPE.SUBCELLULAR_LOCATION + mirrored detector, resolved Candidates sign / DE-direction / cache assumptions in 13-WAVE0-FINDINGS.md, and registered the SubcellularLocation + Volcano test categories.**

## Performance

- **Duration:** ~67 min (includes orchestration recovery: committing untracked plans, retiring blocked worktree subagents, switching to inline sequential)
- **Started:** 2026-05-17T12:53Z (phase init)
- **Completed:** 2026-05-17T13:59Z
- **Tasks:** 3
- **Files modified:** 6 (3 created, 3 modified)

## Accomplishments
- `SEMTYPE.SUBCELLULAR_LOCATION = 'Proteomics-SubcellularLocation'` appended (existing 5 keys unchanged/ordered) + `detectSubcellularLocation` in detectors.js mirroring `detectGeneSymbol` (name-hinted, verbatim literal).
- 13-WAVE0-FINDINGS.md: A2 proves the Candidates sign is internally consistent (DMD row −3.98 = log2(126.58/1991.24), biologically correct) → locks the conditional per-row R3 flip rule (ported from CKomics_tool2.py:1625-1662); A3 pins the D-09 fix to `differential-expression.ts:289-292` as direction-only; A4 selects a per-key userDataStorage map with `__schema_v` invalidation (cites js-api/src/dapi.ts:716-757).
- SubcellularLocation + Volcano test categories registered and wired into package-test.ts for 13-04/13-06 to extend.

## Task Commits

1. **Task 1: Add SEMTYPE + detector** - `a20b9a072e` (feat)
2. **Task 2: Resolve A2/A3/A4 → 13-WAVE0-FINDINGS.md** - `2243d25644` (docs)
3. **Task 3: Register test scaffolds** - `c8c2f256f3` (test)

**Plan metadata:** this commit (docs: complete plan)

## Files Created/Modified
- `src/utils/proteomics-types.ts` - +SUBCELLULAR_LOCATION literal
- `detectors.js` - +detectSubcellularLocation (name-hinted)
- `13-WAVE0-FINDINGS.md` - A2/A3/A4 resolutions with file:line evidence
- `src/tests/subcellular-location.ts`, `src/tests/volcano.ts` - category scaffolds
- `src/package-test.ts` - +2 test imports

## Decisions Made
See key-decisions frontmatter — all three are documented with evidence in 13-WAVE0-FINDINGS.md.

## Deviations from Plan

None - plan executed exactly as written. (The orchestration-level recovery — committing the
previously-untracked phase 13 planning artifacts and switching from blocked parallel-worktree
subagents to inline sequential execution — happened outside this plan's tasks; the three tasks
themselves ran exactly as specified.)

## Issues Encountered
- Pre-task orchestration: phase 13 plans were untracked, so worktree-isolated executors could
  not see them; spawned subagents were additionally denied Bash. Resolved by committing the
  planning artifacts (`8cbe5aab4a`) and executing inline (user-approved). Captured in memory
  `feedback_proteomics_execute_phase_inline`.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Wave 1 sibling 13-02 (enrichment) is independent and can proceed.
- 13-03/13-04/13-05 (Wave 2) are unblocked: R3 rule, cache mechanism, and D-09 site are now locked in 13-WAVE0-FINDINGS.md.
- Runtime `grok test --category SubcellularLocation|Volcano` discovery is deferred to the Wave-1 build/test gate (localhost:8080 server confirmed up, HTTP 200).

---
*Phase: 13-ck-omics-volcano-and-enrichment-parity*
*Completed: 2026-05-17*
