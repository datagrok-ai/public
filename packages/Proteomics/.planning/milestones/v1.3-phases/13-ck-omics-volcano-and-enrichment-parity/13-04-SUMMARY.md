---
phase: 13-ck-omics-volcano-and-enrichment-parity
plan: 04
subsystem: api
tags: [uniprot, subcellular-location, classifier, fetchProxy, userDataStorage, cache, tdd]

requires:
  - phase: 13-ck-omics-volcano-and-enrichment-parity
    provides: SEMTYPE.SUBCELLULAR_LOCATION (13-01), 13-WAVE0-FINDINGS.md A4 (cache mechanism)
provides:
  - parseSubcellularLocation + LOCATION_KEYWORDS + LOCATION_COLORS (verbatim CK-omics contract)
  - mergeStreamTsv (pure positional reviewed-priority parse)
  - getSubcellularLocations (cached, chunked fetchProxy, D-03 reviewed-by-gene fallback)
  - fetchUniProtEntry (cached single-accession JSON; panel delegates here)
affects: [13-06]

tech-stack:
  added: []
  patterns:
    - "Verbatim algorithm port with insertion-order tie-break preserved"
    - "Write-through userDataStorage map cache keyed by __schema_v"
    - "Shared fetch module; panel delegates (single UniProt fetch site)"

key-files:
  created:
    - src/analysis/subcellular-location.ts
  modified:
    - src/panels/uniprot-panel.ts
    - src/tests/subcellular-location.ts

key-decisions:
  - "Classifier/palette/keyword map ported byte-faithfully from CKomics_tool2.py:1384-1436 + locked README — no re-derivation"
  - "Cache = userDataStorage whole-map put/get with __schema_v='13-04-1' (13-WAVE0-FINDINGS A4); stale schema → discard"
  - "Taxonomy filter optional (A1) — omitted when organism unknown (accessions globally unique)"
  - "Unresolved accessions cached as Unknown so they are not re-queried every session"
  - "Panel JSON fetch cached in-memory in the shared module (closes folded cache-uniprot todo)"

patterns-established:
  - "All UniProt I/O via grok.dapi.fetchProxy — zero raw fetch()"
  - "Per-batch try/catch; never throw out of the fetch loop (V5 defensive)"

requirements-completed: [R1]

duration: 234min
completed: 2026-05-17
---

# Phase 13 Plan 04: Subcellular-Location Module + Panel Refactor Summary

**New shared src/analysis/subcellular-location.ts: verbatim CK-omics 11-category classifier + locked hex palette, pure positional reviewed-priority TSV merge, cached chunked fetchProxy stream fetch with the D-03 reviewed-by-gene fallback; uniprot-panel.ts delegates its fetch to the module's cache.**

## Performance

- **Duration:** ~234 min (inline sequential; verbatim-port care + build + test)
- **Started:** 2026-05-17T18:05Z
- **Completed:** 2026-05-17T21:59Z
- **Tasks:** 2 (Task 1 TDD RED → GREEN)
- **Files modified:** 3 (1 created, 2 modified)

## Accomplishments
- `parseSubcellularLocation` / `LOCATION_KEYWORDS` (11 categories, insertion order) / `LOCATION_COLORS` (12 ARGB) ported verbatim from `CKomics_tool2.py:1384-1436` + the locked README. Earliest-position-across-all-categories, subcell-before-GO, strict-`<` insertion-order tie-break, no pre-strip of the raw string.
- `mergeStreamTsv` — pure positional parse (header line skipped, `parts[0..4]`), within-text priority merge (reviewed non-Unknown overwrites; existing Unknown overwritten by non-Unknown).
- `getSubcellularLocations` — cache-first (userDataStorage map, `__schema_v` discard-on-mismatch), ~100/chunk `accession:` OR-query via `grok.dapi.fetchProxy` (optional taxonomy filter), D-03 second pass batching ~20 `gene_exact` genes with `reviewed:true` and back-assigning; per-batch try/catch, never throws out of the loop; write-through cache.
- `fetchUniProtEntry` — cached single-accession JSON; `uniprot-panel.ts` delegates to it (parseAccession still exported there, render path unchanged). Folded cache-uniprot todo closed.
- 10 SubcellularLocation tests green on localhost after `npm run build`; zero raw `fetch(` (fetchProxy only).

## Task Commits

1. **Task 1 RED: classifier/palette/merge tests** - `f720c3d70d` (test)
2. **Task 1 GREEN: verbatim module** - `844c7fdbcf` (feat)
3. **Task 2: panel delegation refactor** - `403cf27b90` (feat)

**Plan metadata:** this commit (docs: complete plan)

## Files Created/Modified
- `src/analysis/subcellular-location.ts` - new shared module (classifier, palette, fetch, cache, D-03)
- `src/panels/uniprot-panel.ts` - fetchUniProtData delegates to fetchUniProtEntry; grok import dropped (unused after delegation)
- `src/tests/subcellular-location.ts` - 9 classifier/palette/merge tests (+ the 13-01 SEMTYPE scaffold test retained)

## Decisions Made
See key-decisions frontmatter. Cache uses the whole-map `userDataStorage.put/get` form (the stored entity IS a per-accession key→value map — consistent with 13-WAVE0-FINDINGS A4 and RESEARCH §Cache Mechanism).

## Deviations from Plan

None - plan executed exactly as written. (`fetchUniProtEntry` typed as `unknown` rather than re-exporting the panel's `UniProtEntry` interface to avoid a cross-module type move; the panel casts at the boundary — render path and types unchanged.)

## Issues Encountered
- Dropped the now-unused `import * as grok` from uniprot-panel.ts after delegation (verified no residual `grok.` usage) — required for a clean lint/build.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- 13-05 (Wave 2 sibling) independent and ready; 13-06 (Wave 3) will consume `getSubcellularLocations` + `LOCATION_COLORS` for volcano coloring.
- Live UniProt stream behavior, real cache round-trip, and the D-03 fallback against the BP rat/mouse engagement data are HUMAN-UAT items (unit tests cover the pure classifier/palette/merge; network paths are not unit-exercised).

---
*Phase: 13-ck-omics-volcano-and-enrichment-parity*
*Completed: 2026-05-17*
