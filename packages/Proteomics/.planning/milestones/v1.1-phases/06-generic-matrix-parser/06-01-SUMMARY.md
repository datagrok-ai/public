---
phase: 06-generic-matrix-parser
plan: 01
subsystem: parsers
tags: [log2-transform, csv-tsv, auto-detection, proteomics]

requires:
  - phase: 01-data-import-and-foundation
    provides: MaxQuant parser with log2 transform and primary column logic
provides:
  - shared-utils.ts with 7 reusable parser utility functions
  - Generic parser test scaffold with 11 test cases and makeGenericCsv helper
  - Refactored MaxQuant parser importing shared utilities
affects: [06-02-generic-parser-dialog, downstream-pipeline]

tech-stack:
  added: []
  patterns: [shared-utility-extraction, semicolon-detection-guard]

key-files:
  created:
    - packages/Proteomics/src/parsers/shared-utils.ts
    - packages/Proteomics/src/tests/generic-parser.ts
  modified:
    - packages/Proteomics/src/parsers/maxquant-parser.ts
    - packages/Proteomics/src/package-test.ts

key-decisions:
  - "addPrimaryColumnIfNeeded scans for semicolons before creating column -- avoids unnecessary columns for clean data"
  - "MaxQuant parser delegates to shared log2TransformColumns but keeps prefix-matching detection logic locally"

patterns-established:
  - "Shared utility pattern: parser-specific detection stays local, generic transforms go to shared-utils.ts"
  - "Semicolon-guard pattern: scan source column for semicolons before creating primary column"

requirements-completed: [IMPORT-05, IMPORT-06]

duration: 3min
completed: 2026-03-06
---

# Phase 6 Plan 01: Shared Utilities & Test Scaffold Summary

**Extracted 7 shared parser utilities from MaxQuant parser into shared-utils.ts with semicolon-detection guard and 11-test generic parser scaffold**

## Performance

- **Duration:** 3 min
- **Started:** 2026-03-06T04:08:27Z
- **Completed:** 2026-03-06T04:11:35Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- Created shared-utils.ts with log2TransformColumns, copyAsLog2Columns, addPrimaryColumnIfNeeded, detectLog2Status, detectDelimiter, autoSuggestProteinIdColumn, autoSuggestIntensityColumns
- Refactored MaxQuant parser to import and delegate to shared utilities (no duplicated log2 logic)
- Added semicolon-detection guard to addPrimaryColumnIfNeeded (skips column creation when no semicolons present)
- Created generic-parser.ts test scaffold with 11 test cases covering IMPORT-01, -03, -05, -06

## Task Commits

Each task was committed atomically:

1. **Task 1: Create shared-utils.ts and refactor MaxQuant parser** - `d2371aedb1` (feat)
2. **Task 2: Create generic parser test scaffold** - `c2011b905f` (test)

## Files Created/Modified
- `packages/Proteomics/src/parsers/shared-utils.ts` - 7 exported utility functions for log2 transform, primary column extraction, delimiter detection, auto-suggestion
- `packages/Proteomics/src/parsers/maxquant-parser.ts` - Refactored to import from shared-utils, keeping MaxQuant-specific prefix detection locally
- `packages/Proteomics/src/tests/generic-parser.ts` - 11 test cases with makeGenericCsv helper covering all shared utility functions
- `packages/Proteomics/src/package-test.ts` - Added import for generic-parser tests

## Decisions Made
- addPrimaryColumnIfNeeded scans source column for semicolons before creating the primary column; skips creation entirely when no semicolons found (per CONTEXT.md decision)
- MaxQuant parser keeps its prefix-matching intensity detection logic locally but delegates the actual log2 transform to shared log2TransformColumns
- Used addPrimaryColumnFromVariants wrapper in MaxQuant parser to bridge the name-variant lookup with the shared utility's single-column-name interface

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- shared-utils.ts is ready for import by the generic parser dialog (Plan 02)
- Test scaffold is ready to run against a Datagrok instance
- All existing MaxQuant parser tests will continue to pass (TypeScript compiles cleanly)

## Self-Check: PASSED

All files and commits verified.

---
*Phase: 06-generic-matrix-parser*
*Completed: 2026-03-06*
