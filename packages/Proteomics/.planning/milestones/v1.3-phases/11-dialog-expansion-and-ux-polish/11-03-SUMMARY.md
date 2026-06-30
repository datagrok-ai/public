---
phase: 11-dialog-expansion-and-ux-polish
plan: 03
subsystem: ui
tags: [viewers, volcano, pca, heatmap, dataframe-naming, ux]

requires:
  - phase: 10-spectronaut-parser-and-core-algorithms
    provides: viewer creation functions and parser infrastructure
provides:
  - Descriptive titles on volcano, PCA, and heatmap viewers
  - Filename-based DataFrame naming across all import paths
  - PCA auxiliary DataFrame named with source table context
affects: []

tech-stack:
  added: []
  patterns: [viewer-title-passthrough, parser-no-name-policy]

key-files:
  created: []
  modified:
    - packages/Proteomics/src/viewers/volcano.ts
    - packages/Proteomics/src/viewers/pca-plot.ts
    - packages/Proteomics/src/viewers/heatmap.ts
    - packages/Proteomics/src/package.ts
    - packages/Proteomics/src/parsers/maxquant-parser.ts
    - packages/Proteomics/src/parsers/spectronaut-parser.ts

key-decisions:
  - "Parsers should not set df.name; naming is the import handler's responsibility"
  - "Volcano title uses group2 vs group1 order matching existing DE comparison direction"
  - "Demo function sets df.name = 'proteinGroups' as fallback since no filename exists"

patterns-established:
  - "viewer-title-passthrough: viewer creation functions accept optional title, callers pass analysis context"
  - "parser-no-name-policy: parsers return unnamed DataFrames; import handlers set df.name from filename"

requirements-completed: [UX-01, UX-02]

duration: 3min
completed: 2026-03-07
---

# Phase 11 Plan 03: Viewer Titles and DataFrame Naming Summary

**Descriptive titles on volcano/PCA/heatmap viewers with filename-based DataFrame naming across all import paths**

## Performance

- **Duration:** 3 min
- **Started:** 2026-03-07T21:38:02Z
- **Completed:** 2026-03-07T21:40:39Z
- **Tasks:** 2
- **Files modified:** 6

## Accomplishments
- All three viewer types (volcano, PCA, heatmap) now display descriptive titles with analysis context
- DataFrame naming moved from parsers to import handlers, ensuring filenames are used as table names
- PCA auxiliary DataFrames include source table name (e.g., "PCA: proteinGroups")
- Demo function provides sensible default name when no filename is available

## Task Commits

Each task was committed atomically:

1. **Task 1: Add title parameter to viewer creation functions** - `819c693f9b` (feat)
2. **Task 2: Wire viewer titles in package.ts and fix DataFrame naming** - `0c1e5e1eaa` (feat)

## Files Created/Modified
- `packages/Proteomics/src/viewers/volcano.ts` - Added optional title to options, set via setOptions
- `packages/Proteomics/src/viewers/pca-plot.ts` - Added optional title parameter, set via setOptions
- `packages/Proteomics/src/viewers/heatmap.ts` - Added optional title to options, set via setOptions
- `packages/Proteomics/src/package.ts` - Pass titles to all viewer calls, set pcaDf.name, set demo df.name
- `packages/Proteomics/src/parsers/maxquant-parser.ts` - Removed df.name = 'proteinGroups'
- `packages/Proteomics/src/parsers/spectronaut-parser.ts` - Removed df.name = 'spectronaut'

## Decisions Made
- Parsers should not set df.name; naming is the import handler's responsibility (prevents overwriting filename-based names)
- Volcano title uses group2 vs group1 order matching the existing DE comparison direction
- Demo function sets df.name = 'proteinGroups' as a sensible fallback since no filename exists in the demo path

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 11 Plan 03 completes the viewer titles and DataFrame naming requirements (UX-01, UX-02)
- All viewer creation functions are backward-compatible (title parameters are optional)

---
*Phase: 11-dialog-expansion-and-ux-polish*
*Completed: 2026-03-07*
