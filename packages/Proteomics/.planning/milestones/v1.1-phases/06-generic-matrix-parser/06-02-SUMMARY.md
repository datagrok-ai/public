---
phase: 06-generic-matrix-parser
plan: 02
subsystem: parsers
tags: [column-mapping, dialog, auto-suggestion, log2-detection, csv-tsv, generic-import]

requires:
  - phase: 06-generic-matrix-parser
    provides: shared-utils.ts with log2 transform, auto-suggestion, delimiter detection
provides:
  - Generic matrix import dialog with column mapping, auto-suggestion, log2 detection, live preview
  - Menu entry 'Proteomics | Import | Generic Matrix...'
  - MaxQuant importer backported with source tags and file-based naming
affects: [downstream-pipeline, 07-qc-dashboard]

tech-stack:
  added: []
  patterns: [file-picker-dialog, reactive-column-mapping, bigint-column-support]

key-files:
  created:
    - packages/Proteomics/src/parsers/generic-parser.ts
  modified:
    - packages/Proteomics/src/package.ts
    - packages/Proteomics/src/parsers/shared-utils.ts

key-decisions:
  - "Used DG.Utils.openFile() for file selection followed by separate mapping dialog"
  - "Log2 toggle value = !isLog2 (ON means needs transform, OFF means already log2)"
  - "All numeric columns shown as available intensity columns, not just keyword matches"
  - "Added autoSuggestGeneNameColumn() for gene name auto-suggestion"

patterns-established:
  - "File import pattern: openFile -> parse -> show mapping dialog -> assign semTypes -> addTableView"
  - "BigInt column support: Number() conversion required for BIG_INT column values in shared utilities"

requirements-completed: [IMPORT-01, IMPORT-02, IMPORT-03, IMPORT-04, IMPORT-05, IMPORT-06]

duration: 16h
completed: 2026-03-06
---

# Phase 6 Plan 02: Generic Import Dialog Summary

**Column mapping dialog with auto-suggestion, log2 detection, live preview, and BigInt support for importing any CSV/TSV proteomics matrix into the full analysis pipeline**

## Performance

- **Duration:** ~16h (across checkpoint pause for human verification)
- **Started:** 2026-03-06T04:14:00Z
- **Completed:** 2026-03-06T20:08:00Z
- **Tasks:** 3
- **Files modified:** 3

## Accomplishments
- Created generic-parser.ts with showGenericImportDialog() featuring file picker, column mapping dialog, auto-suggestion, log2 detection, live preview grid, and full import handler
- Added 'Proteomics | Import | Generic Matrix...' menu entry in package.ts
- Backported MaxQuant importer with df.name from filename and proteomics.source = 'maxquant' tag
- Fixed BigInt column support in shared-utils.ts (Number() conversion for BIG_INT values)
- Added autoSuggestGeneNameColumn() for gene name auto-suggestion
- All numeric columns (including BIG_INT) shown as available intensity columns

## Task Commits

Each task was committed atomically:

1. **Task 1: Build generic-parser.ts with column mapping dialog** - `eee98e17aa` (feat)
2. **Task 2: Add menu entry and backport MaxQuant importer** - `110d9fd4a0` (feat)
3. **Task 3: Verification fixes (BigInt support, auto-suggest, preview)** - `7bb0e9ab26` (fix)

## Files Created/Modified
- `packages/Proteomics/src/parsers/generic-parser.ts` - Generic matrix import dialog with file picker, column mapping, auto-suggestion, log2 detection, live preview, and import handler
- `packages/Proteomics/src/package.ts` - Added Generic Matrix menu entry, backported MaxQuant source tags and file naming
- `packages/Proteomics/src/parsers/shared-utils.ts` - BigInt column support (Number() conversion), added autoSuggestGeneNameColumn()

## Decisions Made
- Used DG.Utils.openFile() for file selection followed by a separate column mapping dialog (not an integrated file-and-map dialog)
- Log2 toggle semantics: ON = needs transform (raw intensities), OFF = already log2-transformed
- All numeric columns shown as available for intensity selection, with keyword-matched columns pre-selected as suggestions
- Added autoSuggestGeneNameColumn() to auto-suggest gene name column (keywords: gene, symbol)
- Number() conversion added for BigInt column values throughout shared-utils to prevent type errors

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] BigInt column type support in shared-utils and dialog**
- **Found during:** Task 3 (verification)
- **Issue:** BIG_INT columns caused type errors -- log2TransformColumns, copyAsLog2Columns, and detectLog2Status did not handle BigInt values
- **Fix:** Added DG.COLUMN_TYPE.BIG_INT checks and Number() conversion for BigInt values in all shared utility functions
- **Files modified:** packages/Proteomics/src/parsers/shared-utils.ts
- **Committed in:** 7bb0e9ab26

**2. [Rule 1 - Bug] All numeric columns available as intensity candidates**
- **Found during:** Task 3 (verification)
- **Issue:** Only keyword-matched columns appeared as available intensity columns, but users need to select any numeric column
- **Fix:** Changed intensity columns input to show all numeric columns as available, with keyword matches as the pre-selected default value
- **Files modified:** packages/Proteomics/src/parsers/generic-parser.ts
- **Committed in:** 7bb0e9ab26

**3. [Rule 2 - Missing Critical] Gene name auto-suggestion**
- **Found during:** Task 3 (verification)
- **Issue:** Gene name column had no auto-suggestion, making it harder for users to identify the correct column
- **Fix:** Added autoSuggestGeneNameColumn() to shared-utils and wired it into the dialog
- **Files modified:** packages/Proteomics/src/parsers/shared-utils.ts, generic-parser.ts
- **Committed in:** 7bb0e9ab26

**4. [Rule 1 - Bug] Hint message visibility and preview sizing**
- **Found during:** Task 3 (verification)
- **Issue:** Hint message not visible, preview grid too constrained
- **Fix:** Fixed hint div styling and removed column limit for preview
- **Files modified:** packages/Proteomics/src/parsers/generic-parser.ts
- **Committed in:** 7bb0e9ab26

---

**Total deviations:** 4 auto-fixed (3 bugs, 1 missing critical)
**Impact on plan:** All fixes necessary for correct operation with real-world data. No scope creep.

## Issues Encountered
- BigInt column values required explicit Number() conversion -- Datagrok's BIG_INT type returns BigInt values that cannot be passed directly to Math.log2()
- The `value: null` parameter for ui.input.column required `undefined` instead (TypeScript type mismatch)

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 6 is now complete -- generic matrix import dialog is fully functional
- All IMPORT requirements (IMPORT-01 through IMPORT-06) are satisfied
- Ready for Phase 7 (QC Dashboard) which depends on the import infrastructure
- Any CSV/TSV proteomics matrix can now enter the full downstream pipeline (annotation, normalization, imputation, DE)

## Self-Check: PASSED

All files and commits verified.

---
*Phase: 06-generic-matrix-parser*
*Completed: 2026-03-06*
