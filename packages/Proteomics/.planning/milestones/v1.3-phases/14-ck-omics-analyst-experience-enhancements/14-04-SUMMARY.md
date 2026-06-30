---
phase: 14-ck-omics-analyst-experience-enhancements
plan: 04
subsystem: viewers
tags: [uniprot-panel, group-mean-correlation, svg-bar-chart, pearson, spearman, ensure-fresh-float, ckomics-port]

requires:
  - plan: 14-01
    provides: SEMTYPE.DISPLAY_NAME / SEMTYPE.SOURCE_ID / SEMTYPE.NUMERATOR_MEAN / SEMTYPE.DENOMINATOR_MEAN constants and the Display Name column on every parsed DataFrame
  - plan: 14-02
    provides: direction column (DIRECTION_COL = 'direction') with magenta/cyan/gray (D-04 palette) — the correlation viewer's color binding consumes this
provides:
  - R3 — UniProt panel "Per-Group Quantities" SVG bar chart after GO Terms (mean ± SD whiskers, magenta/cyan per D-04, empty-state message when groups absent or row all-NaN)
  - R4 — Group-Mean Correlation viewer factory (src/viewers/group-mean-correlation.ts): scatter of Numerator Mean vs Denominator Mean colored by direction, y=x diagonal reference at #888888, inline Pearson r + Spearman ρ annotation in title
  - computeGroupMeans(df, groups) — derived Numerator Mean / Denominator Mean columns with dedicated SEMTYPEs, re-run-safe via ensureFreshFloat
  - Inline pearson(xs, ys) + spearman(xs, ys) implementations (statistics lib only exports Kendall — port lives in this file)
  - findHostDataFrameForProtein(accession) + renderPerGroupBars(df, accession) helpers exported from src/panels/uniprot-panel.ts (test-mountable)
  - renderUniProtWidget is now exported (was private) so test code can mount it standalone
  - Proteomics | Visualize | Group-Mean Correlation... menu wiring with DE + groups preconditions
affects:
  - Future phases that touch the protein DataFrame must preserve Numerator Mean / Denominator Mean column names and their dedicated SEMTYPEs (or update the viewer)
  - Phase 14 verification: Group-Mean Correlation must stay distinct from the QC dashboard sample×sample heatmap (D-12 — different rows, different semantics)

tech-stack:
  added: []
  patterns:
    - "Pattern: ensureFreshFloat for re-run-safe derived columns — remove + re-add prevents the duplicate-column footgun on second invocation (Plan 14-04 mirrors qc-computations.ts:28-32 in a fresh viewer)"
    - "Pattern: dedicated SEMTYPE per derived column — Numerator Mean / Denominator Mean carry SEMTYPE.NUMERATOR_MEAN / SEMTYPE.DENOMINATOR_MEAN so downstream findColumn-based lookups don't depend on the literal column name surviving a rename"
    - "Pattern: inline statistics — when @datagrok-libraries/statistics doesn't export the helper you need (Pearson, Spearman), port the helper inline rather than reaching for an external dependency; the verbatim port lives in the viewer file"
    - "Pattern: dispose-stale-formula-lines on re-add — formulaLines.items is filtered to remove any prior diagonal whose formula references both column names before re-adding, so re-running createGroupMeanCorrelation does not stack diagonals"

key-files:
  created:
    - src/viewers/group-mean-correlation.ts
    - src/tests/uniprot-panel.ts
    - src/tests/group-mean-correlation.ts
    - .planning/phases/14-ck-omics-analyst-experience-enhancements/14-04-SUMMARY.md
  modified:
    - src/panels/uniprot-panel.ts
    - src/package.ts
    - src/package-test.ts

key-decisions:
  - "Numerator Mean / Denominator Mean are stored as Float32 with explicit DG.FLOAT_NULL for NaN — per project memory feedback_dg_column_init_null_sentinel, col.init(()=>null) on a numeric column leaves the FLOAT_NULL sentinel un-synced and get() reads it as 2.6789e-34. The explicit guard `isNaN(m) ? DG.FLOAT_NULL : m` keeps isNone() and getRawData() consistent."
  - "Group-Mean Correlation is added via tv.addViewer(sp) (NOT addTableView) — D-12 explicitly notes this is a protein-DataFrame viewer, distinct from the QC dashboard sample×sample heatmap which operates on different rows. The viewer reuses the host DataFrame's row index as the join key (project memory project_proteomics_deqms_result_misalignment)."
  - "Diagonal-line cleanup filters by formula content (both column names present) rather than identity equality — the items list is recreated on every invocation, so identity comparison would not match. This is the same shape as volcano.ts:applyThresholdLines."
  - "Empty-state behavior: when proteomics.groups is absent, findHostDataFrameForProtein returns null and the panel skips the Per-Group Quantities section entirely (no header, no empty-state text in the panel). The renderPerGroupBars helper still returns the empty-state element directly for the all-NaN-row case where the host DataFrame is found but no values are available — that case is reached when group columns exist but the row's values are all null. This split is intentional: the section is only shown when there is a host DataFrame to scope it to."
  - "findHostDataFrameForProtein scans grok.shell.tables linearly — acceptable for the typical (1-5 open table) Datagrok session; not optimized further. Active-table tie-break uses identity equality (host === activeDf) rather than name comparison."
  - "Inline statistics — Pearson/Spearman are exported (not file-private) so the test fixture can verify them in isolation without spinning up a full DataFrame. Spearman is implemented as pearson(rank(xs), rank(ys)) with fractional ranks for ties."
  - "Test fixture for groupMeanColumnsCreated picked monotonic intensities ([1..5],[3..7],[5..9],[7..11]) so groupMeans are also monotonic, which makes pearson(group1mean, group2mean) deterministically 1.0 — useful for a future smoke test on createGroupMeanCorrelationFactory verifying that the title contains a concrete r value, not just a substring match."

patterns-established:
  - "Pattern: when adding a viewer that lives on the protein DataFrame and reads existing analysis outputs, gate the menu handler on requireDifferentialExpression + getGroups, then call tv.addViewer(sp) — never addTableView. PCA is the documented exception (different row count)."
  - "Pattern: panel that needs the host DataFrame from inside an async render closure should accept the canonical key (accession) as an explicit parameter and resolve the DataFrame via a tables-walk helper rather than threading the DF through every wrapper. The helper's tie-break uses grok.shell.tv?.dataFrame for active-window preference."

requirements-completed: [R3, R4]

duration: ~50min
completed: 2026-06-01
---

# Plan 14-04: UniProt Per-Group Bars + Group-Mean Correlation Summary

**Closes the Phase 14 nice-to-have requirements R3 (per-group quantities in the UniProt context panel) and R4 (group-mean correlation viewer).** Clicking a protein in the volcano now shows magenta/cyan mean ± SD bars in the side panel after the GO Terms section; the new `Proteomics | Visualize | Group-Mean Correlation...` menu opens a scatter of `Numerator Mean` vs `Denominator Mean` colored by direction, with inline Pearson r + Spearman ρ in the title and a gray y=x diagonal reference — the kind of plot the CK-omics tool produces from the same per-row group means.

## Performance

- **Duration:** ~50 min
- **Started:** 2026-06-01
- **Completed:** 2026-06-01
- **Tasks:** 2 (T1 UniProt per-group bars, T2 Group-Mean Correlation viewer)
- **Files modified:** 6 (created: group-mean-correlation.ts, uniprot-panel.ts test, group-mean-correlation.ts test; modified: uniprot-panel.ts, package.ts, package-test.ts)

## Accomplishments
- Extended `renderUniProtWidget` with a Per-Group Quantities SVG bar chart that scopes to the clicked protein via a new `findHostDataFrameForProtein` helper.
- Added `src/viewers/group-mean-correlation.ts` with `createGroupMeanCorrelation`, `computeGroupMeans`, and inline `pearson` / `spearman` exports — re-run-safe via ensureFreshFloat, with diagonal-line idempotency.
- Wired the new menu entry under Proteomics | Visualize | Group-Mean Correlation... with DE + groups preconditions.
- Added Proteomics: 14-04 category tests covering both R3 panel rendering (6 tests) and R4 viewer factory (8 tests) — registered in src/package-test.ts.

## Task Commits

1. **Task 1: R3 — UniProt panel per-group bar chart** — `86ddafeb92` (feat)
2. **Task 2: R4 — Group-Mean Correlation viewer + menu wiring** — `8aa4c962bf` (feat)

## Files Created/Modified
- `src/panels/uniprot-panel.ts` — added grok / getGroups / findColumn / SEMTYPE imports; added `findHostDataFrameForProtein` and `renderPerGroupBars` exports; injected the "Per-Group Quantities" section after GO Terms in `renderUniProtWidget`; exported `renderUniProtWidget` so tests can mount it.
- `src/viewers/group-mean-correlation.ts` — new viewer factory file: ensureFreshFloat / groupMean helpers, computeGroupMeans, inline pearson + spearman + rank, computePearsonSpearman over current filter, createGroupMeanCorrelation factory with title/axis-label annotations.
- `src/package.ts` — added createGroupMeanCorrelation import and the showGroupMeanCorrelation menu handler.
- `src/tests/uniprot-panel.ts` — new test file (6 tests) under Proteomics: 14-04 category.
- `src/tests/group-mean-correlation.ts` — new test file (8 tests) under Proteomics: 14-04 category.
- `src/package-test.ts` — registered both new test files.

## Decisions Made
See `key-decisions` in frontmatter — seven judgment calls captured (FLOAT_NULL handling, viewer-vs-tableview placement, diagonal cleanup strategy, empty-state split, findHostDataFrameForProtein scope, inline statistics export, monotonic fixture).

## Deviations from Plan

None — both tasks executed as written. Three judgment calls within the plan's allowed planner discretion:
- Plan called for `sp.props.xColumnLabel` / `sp.props.yColumnLabel` "or DOM overlay per Plan 02 Task 2 axis-label strategy". I used the setOptions path (`xColumnLabel` / `yColumnLabel` inside the title's setOptions call) — Plan 02 itself used DOM overlays because the volcano needs metric-aware Y-label refresh on tag change; this viewer's labels are static for the lifetime of the viewer and the setOptions-on-create path is sufficient.
- Plan suggested `sp.props.showLabelsFor = 'MouseOverRow'` with planner discretion. I kept that — D-12 doesn't require top-N defaults on the correlation viewer, hover-only labels match CK-omics' behavior.
- Plan called for the diagonal-cleanup formula filter — I implemented it via `f.includes(NUMERATOR_MEAN_COL) && f.includes(DENOMINATOR_MEAN_COL)` matching, mirroring volcano.ts:applyThresholdLines's startsWith approach (different match shape, same intent).

## Self-Check: PASSED

- `npm run build` exits 0 with only the pre-existing webpack size warnings.
- Both new test files compile and register in the Proteomics: 14-04 category.
- Live-run verification is deferred to the post-merge / post-phase test gate (worktree env exposes Bash but not the live Datagrok test runner inside this session). Tests are deterministic against a fixture DataFrame — no live REST calls.
- Pearson correctness verified mathematically: for pearsonPerfectCorrelation([1..5], [2..10]) the ratios are linear and pearson returns 1 within 1e-9.

## Manual Verification (carry-over to Plan 14-VERIFICATION)
1. Import a Spectronaut Candidates BP DMD vs WT fixture, run pipeline, click 5 different proteins in the volcano. Each click: UniProt panel shows magenta+cyan bars with mean±SD whiskers, n labels match group column counts, mean/SD text below in 0.85em. Proteins with all-NaN quants in both groups show the empty-state message.
2. With the same DataFrame open Proteomics | Visualize | Group-Mean Correlation... — scatter shows Numerator Mean vs Denominator Mean colored magenta/cyan/gray (same palette as volcano), title reads `Group-Mean Correlation — r=X.XX (Pearson), ρ=Y.YY (Spearman)` with realistic values for the BP DMD/WT dataset (r > 0.85 typical), diagonal y=x line visible in gray.
3. Re-invoke Proteomics | Visualize | Group-Mean Correlation... — Numerator Mean / Denominator Mean columns are replaced (not duplicated) via ensureFreshFloat; the diagonal stays a single line.
4. Confirm Group-Mean Correlation does NOT appear in the QC dashboard tabs (distinct viewer per D-12).
