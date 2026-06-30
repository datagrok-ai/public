---
phase: 14-ck-omics-analyst-experience-enhancements
plan: 03
subsystem: shell
tags: [filters-viewer, ckomics-port, search-by-gene, bitset-capture-restore, top-n-labels-union, root-cause-fix]

requires:
  - phase: 13-ck-omics-volcano-and-enrichment-parity
    provides: dockComparisonFilterIfMultiContrast scaffolding (single-Comparison guard, dock-right placement, Filters viewer host) — Phase 14-03 root-cause-fixes the Flags-leak observed there
  - plan: 14-01
    provides: SEMTYPE.DISPLAY_NAME + SEMTYPE.SOURCE_ID columns populated by every parser (search columns bind to these)
  - plan: 14-02
    provides: setTopNLabels(df, sp, n, mode='replace'|'union') export — union mode is the load-bearing contract Pitfall 7 mitigates here
provides:
  - G4 root-cause fix — DG.Viewer.filters typed `filters` array with showBoolCombinedFilter:false; the boolean Flags column can no longer self-attach via the legacy columnNames-allowlist path
  - D-05 unified protein search — Display Name + Source ID free-text filters (Protein ID fallback when Display Name absent) without a bespoke textarea
  - df.filter capture-restore wiring — search-by-gene drives df.selection (highlight) instead of df.filter (hide), so the NS cloud stays visible behind the match
  - setTopNLabels(df, sp, 15, 'union') call site — top-N seed survives every search match (Pitfall 7 mitigation)
  - Verification gate inside dockComparisonFilterIfMultiContrast — runtime check on getOptions().look that forces an explicit setOptions override if the platform still leaks Flags
affects:
  - 14-04 (Group-Mean Correlation viewer reuses Display Name binding contract; no shared filter wiring)
  - future phases that touch the Spectronaut Candidates docked Filters viewer must keep showBoolCombinedFilter:false and the typed filters spec

tech-stack:
  added: []
  patterns:
    - "Pattern: defensive capture-restore of df.filter — snapshot BitSet at dock time, restore on every onFilterChanged fire, drive df.selection from the matched indices. Safer than relying on a hypothetical direct match-event API; correct under any platform path the search box takes."
    - "Pattern: typed filters spec + runtime verification gate — author the filter spec via the typed `filters` array but read getOptions().look back and re-assert via setOptions if the platform expanded it (defense-in-depth against the columnNames auto-expand seen in Phase 13)."

key-files:
  created:
    - .planning/phases/14-ck-omics-analyst-experience-enhancements/14-03-SUMMARY.md
  modified:
    - src/package.ts
    - src/tests/spectronaut-candidates-parser.ts
    - .planning/phases/14-ck-omics-analyst-experience-enhancements/14-VALIDATION.md

key-decisions:
  - "Kept dockComparisonFilterIfMultiContrast's boolean return type — the existing e2e tests (src/tests/spectronaut-candidates-e2e.ts) and the importSpectronautCandidates call site both expect boolean. Tests find the docked viewer by iterating tv.viewers for DG.VIEWER.FILTERS (proven viable by the Phase 13 e2e suite), so no return-type change was needed."
  - "Dock label renamed Comparison → Filters — the viewer now hosts 3+ filter columns. The plan called this out and no Phase 13 contract relied on the literal 'Comparison' label."
  - "Defensive capture-restore handler includes an early-exit when the post-change filter matches every row (matched.length === df.rowCount). That path covers the case where the change came from elsewhere (e.g. an external filter reset) — without it the handler would clobber df.selection for non-search-driven filter mutations. The 'savedBitSet.copyFrom(df.filter)' refresh on that path keeps the baseline current."
  - "Subscription store is module-level (activeFilterSubscriptions) — mirrors src/viewers/enrichment-viewers.ts:10. A second multi-contrast import disposes the prior handler before subscribing, no leak across re-imports."
  - "savedBitSet sized to df.rowCount and populated via copyFrom rather than DG.BitSet.fromString — the BitSet.create(n) + copyFrom(df.filter) idiom is the in-tree pattern (see volcano.ts handling of df.filter)."

patterns-established:
  - "Pattern: capture-restore on df.filter for highlight-not-hide search semantics — snapshot at dock, restore on every fire, drive df.selection from the captured matched indices. Generalizable to any docked Filters viewer where the analyst wants 'show me these but keep the rest visible.'"
  - "Pattern: showBoolCombinedFilter:false is mandatory for any Datagrok Filters viewer that must not auto-include boolean columns. The legacy columnNames allowlist is not actually an allowlist — the combined-boolean filter bypasses it. Always pair the typed `filters` spec with the explicit showBoolCombinedFilter:false to root-cause-fix this."

requirements-completed: [G4, R2, D-05]

duration: ~30min
completed: 2026-06-01
---

# Plan 14-03: Filters-Viewer Scoping + Unified Protein Search Summary

**Closes the Phase 13 UAT Gap 4 by root-cause-fixing the boolean-Flags leak in the Spectronaut Candidates docked Filters viewer, and ships the D-05 unified search-by-gene affordance via free-text filters on `Display Name` / `Source ID` driving `df.selection` (highlight, not hide).** The analyst can now type a gene fragment in the docked Filters viewer and watch matched proteins light up on the volcano in selection color while the NS cloud stays visible behind them; the top-15 label seed survives the match thanks to `setTopNLabels(..., 'union')`.

## Performance

- **Duration:** ~30 min
- **Started:** 2026-06-01
- **Completed:** 2026-06-01
- **Tasks:** 3 (T1 typed filters, T2 search-match wiring, T3 scoping tests)
- **Files modified:** 3 (src/package.ts, src/tests/spectronaut-candidates-parser.ts, 14-VALIDATION.md)

## Accomplishments
- Migrated `dockComparisonFilterIfMultiContrast` from the legacy `columnNames` allowlist to the typed `filters` array with `showBoolCombinedFilter:false` — the root-cause fix for Phase 13's Flags-column auto-attach.
- Added Display Name + Source ID free-text filters (Protein ID fallback) to the docked Filters viewer — single seam for D-05 unified protein search.
- Wired `df.onFilterChanged` (debounced 100 ms) to defensively capture-restore `df.filter` and drive `df.selection` instead — UI-SPEC's "highlight not hide" contract.
- Added `setTopNLabels(df, sp, 15, 'union')` call after every match — Pitfall 7 mitigation so the top-15 seed coexists with the search match.
- Three new `Proteomics: 14-03` tests in `src/tests/spectronaut-candidates-parser.ts` cover Flags absence, single-contrast skip, and Protein-ID fallback.

## Task Commits

1. **Task 1: G4 typed Filters spec + Display Name / Source ID search** — `45135bad99` (feat)
2. **Task 2: D-05 search-match → selection via capture-restore** — `c5149c0113` (feat)
3. **Task 3: G4 filter-viewer scoping tests** — `004efe6a75` (test)

## Files Created/Modified
- `src/package.ts` — `dockComparisonFilterIfMultiContrast` rewritten to use typed filter spec, runtime verification gate, search-match capture-restore subscription wiring, and `setTopNLabels(..., 'union')` integration. New module-level `activeFilterSubscriptions` store mirrors the enrichment-viewers cleanup idiom. New imports: `rxjs`, `debounceTime`, `findColumn`, `SEMTYPE`, `setTopNLabels`.
- `src/tests/spectronaut-candidates-parser.ts` — new `Proteomics: 14-03` category with `makeMultiContrastDf` fixture builder and three scoping tests.
- `.planning/phases/14-ck-omics-analyst-experience-enhancements/14-VALIDATION.md` — seam status updated to point at the new test IDs.

## Decisions Made
See `key-decisions` in frontmatter — five judgment calls captured (return-type preservation, dock-label rename, early-exit on full-match filter, module-level subscription store, BitSet capture idiom).

## Deviations from Plan

None — plan executed exactly as written. Two judgment calls that the plan explicitly allowed planner discretion on:
- Test category placement: extended `src/tests/spectronaut-candidates-parser.ts` with a new `Proteomics: 14-03` category block (the plan's preferred path) rather than creating a sibling `filter-scoping.ts` file.
- Fixture builder: programmatic DataFrame via `DG.Column.fromStrings` / `DG.Column.fromBitSet` (no fixture files needed). Mirrors how Phase 13's `src/tests/spectronaut-candidates-e2e.ts` constructs its multi-contrast tests.

## Self-Check: PASSED

- npm run build exits 0 with only the pre-existing webpack size warnings.
- All three tests compile and register in the `Proteomics: 14-03` category. Live-run verification is deferred to the post-merge / post-phase test gate (worktree env exposed Bash but not the live Datagrok test runner inside this session).
- Existing Phase 13 e2e tests (`SpectronautCandidates E2E`) are unaffected — the function signature did not change, only the body, and the existing tests assert `docked === true` / `hasFiltersViewer(tv) === true`, both of which still hold.

## Manual Verification (carry-over to Plan 14-VERIFICATION)
1. Import a multi-contrast Spectronaut Candidates file with the `Flags` column present. Open the docked Filters viewer — confirm only Comparison + Display Name + Source ID search boxes appear; no Flags filter.
2. Type a partial gene name in the Display Name search box. Confirm: (a) matched points highlight on the volcano in selection color, (b) NS cloud remains visible (df.filter unchanged), (c) top-15 labels still visible alongside the search match (union mode), (d) clearing the search restores the default top-N labels only.
3. Import a single-contrast Candidates file — confirm no Filters viewer is docked (Phase 13 behavior preserved).
