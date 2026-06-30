---
phase: 13-ck-omics-volcano-and-enrichment-parity
plan: 02
subsystem: api
tags: [enrichment, g-profiler, wikipathways, viewers, dockmanager, tdd]

requires:
  - phase: 09-enrichment-visualization-volcano-integration
    provides: buildEnrichmentDf shape, wireEnrichmentToVolcano cross-link, dock idiom
provides:
  - splitGenesByDirection (up/down by fc sign, shared all-detected background)
  - buildEnrichmentDf Direction-label support; merged directional enrichment DataFrame
  - WikiPathways (WP) default-on g:Profiler source
  - side-by-side Up/Down dot+bar enrichment layout
affects: [13-06, 14]

tech-stack:
  added: []
  patterns:
    - "Pure splitGenesByDirection extraction → unit-testable directional split"
    - "DataFrame.append to merge per-direction enrichment frames; re-apply tag/name/color after append"

key-files:
  created: []
  modified:
    - src/analysis/enrichment.ts
    - src/viewers/enrichment-viewers.ts
    - src/tests/enrichment.ts

key-decisions:
  - "Strict fc>0 / fc<0 split (CK-omics run_gprofiler_analysis); background = all detected, shared by both calls"
  - "Skip a direction's g:GOSt call when it has zero genes (no empty-query request)"
  - "append() returns a new frame → re-apply proteomics.enrichment tag, name, FDR color-coding"
  - "Single-direction or no Direction column falls back to original 2-viewer layout"

patterns-established:
  - "Directional enrichment: one g:GOSt per non-empty direction, identical custom background"
  - "Phase-9 cross-link preserved verbatim; Direction is an additive categorical only"

requirements-completed: [R2, R5]

duration: 140min
completed: 2026-05-17
---

# Phase 13 Plan 02: Up/Down Enrichment Split + WikiPathways Summary

**runEnrichmentPipeline now runs one g:Profiler query per direction (up=fc>0 / down=fc<0) over a shared all-detected background, merges into one Direction-tagged DataFrame, renders Up/Down dot+bar side-by-side, and adds WikiPathways as a default-on source — Phase-9 cross-link untouched.**

## Performance

- **Duration:** ~140 min (inline sequential; includes full build + Enrichment test gate)
- **Started:** 2026-05-17T14:00Z
- **Completed:** 2026-05-17T16:19Z
- **Tasks:** 2 (Task 1 TDD: RED → GREEN)
- **Files modified:** 3

## Accomplishments
- `splitGenesByDirection()` — pure, unit-tested: strict `fc>0`/`fc<0`, shared all-detected background, null fc/p → background-only.
- `buildEnrichmentDf(..., direction?)` adds a categorical `Direction` column (bulk init); pipeline issues ≤2 g:GOSt calls (skips an empty direction), merges via `append()`, re-applies enrichment tag/name/FDR color.
- `WikiPathways` (`'WP'` literal) bool input, default true, pushed into `selectedSources`; gGOSt request body otherwise byte-unchanged.
- `openEnrichmentVisualization` docks Up (dot/bar) and Down (dot/bar) side-by-side when both directions present; single-direction → original layout. `wireEnrichmentToVolcano` signature/body unchanged.
- 5 new Enrichment tests; full suite green (12 Enrichment + 8 Enrichment Visualization, exit 0) on localhost after `npm run build`.

## Task Commits

1. **Task 1 RED: failing directional-split tests** - `ddb0c1e5d2` (test)
2. **Task 1 GREEN: split + Direction + WP source** - `e0bd8140d4` (feat)
3. **Task 2: side-by-side Up/Down charts + cross-link test** - `8c35d55bce` (feat)

**Plan metadata:** this commit (docs: complete plan)

## Files Created/Modified
- `src/analysis/enrichment.ts` - splitGenesByDirection, Direction param, two-call pipeline, WP input
- `src/viewers/enrichment-viewers.ts` - filterByDirection helper + split dock layout
- `src/tests/enrichment.ts` - 5 tests (split ×2, Direction ×2, cross-link regression)

## Decisions Made
See key-decisions frontmatter. All match CK-omics `run_gprofiler_analysis` (4450-4684) and `get_selected_databases` (478-486).

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- `timeout`/`gtimeout` not on PATH (macOS) — ran `grok test` under the tool's own timeout instead. No functional impact.
- Verification batched: the per-task `<verify>` is `npm run build && grok test --category Enrichment`; ran once after both tasks (one authoritative GREEN gate) rather than 4× inline, consistent with the execute-phase post-wave build/test gate. Result: build compiled successfully, all targeted tests pass.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Wave 1 complete (13-01 + 13-02). Wave 2 (13-03/04/05) unblocked.
- Visual 4-viewer side-by-side layout asserted logically + unit-tested for the cross-link; the actual docked visual arrangement is a HUMAN-UAT item (requires a real enrichment run against g:Profiler).

---
*Phase: 13-ck-omics-volcano-and-enrichment-parity*
*Completed: 2026-05-17*
