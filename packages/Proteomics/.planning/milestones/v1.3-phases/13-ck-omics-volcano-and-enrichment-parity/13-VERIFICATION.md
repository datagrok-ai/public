---
phase: 13-ck-omics-volcano-and-enrichment-parity
verified: 2026-06-03T00:00:00Z
status: human_needed
score: 6/6
overrides_applied: 0
re_verification:
  previous_status: human_needed
  previous_score: 6/6
  gaps_closed:
    - "13-UAT Test 2 (major) — dot/bar enrichment charts visually associated with the wrong table; protein-view-targeted Option A docking landed in 13-09"
    - "13-UAT Test 3 (minor) — first-time Color-by-Subcellular Location feels like a hang; in-volcano busy overlay landed in 13-10"
  gaps_remaining: []
  regressions: []
  notes:
    - "13-09 reverts the 13-07 Option B 'Volcano (linked)' co-dock; this is the Option A → Option B → Option A oscillation flagged in the report below"
    - "All claimed automated coverage verified by file inspection: Enrichment Visualization category has 13 tests (5 topN + 3 wireEnrichmentToVolcano + 3 layout B/C/D + 2 banner); Volcano category has 9 tests (7 baseline + overlay lifecycle + warm-cache short-circuit); SubcellularLocation has 15 tests as preserved by 13-10 Task 3 grep gates"
human_verification:
  - test: "Run Proteomics | Analyze | Enrichment on a DE-complete protein table; verify the user lands on the protein view (not the enrichment view) and the Up/Down dot+bar viewers are docked NEXT TO the existing volcano on that view"
    expected: "After OK in the Enrichment dialog: grok.shell.v points at the protein TableView; the original volcano is visible; Up — Dot Plot / Up — Bar Chart / Down — Dot Plot / Down — Bar Chart are docked RIGHT of the volcano on the same view; the merged enrichment grid is accessible as a separate tab; clicking a term row on either the dot/bar (protein view) or the merged enrichment grid (enrichment tab) highlights the matching proteins on the original volcano in the focused view"
    why_human: "Visual layout + focus on a live DG.DockManager-driven multi-view shell cannot be verified by grep; the layout regression test asserts viewer counts and Dart-handle identity but not the visible 'lands next to volcano' user perception that 13-UAT Test 2 verbatim called out"
  - test: "On a real ~8k-accession Spectronaut Candidates session: open Volcano Options, switch Color by to Subcellular Location for the first time, and watch the volcano viewer for an in-place busy overlay showing live phase + N/total progress"
    expected: "A centered white card on sp.root reads 'Classifying subcellular locations…' then advances through 'Fetching subcellular locations: N/80 (P%)' as chunks complete; the overlay updates at least once per chunk completion; it disappears when the column finishes populating. On a warm-cache toggle of the same session the overlay either does not appear at all or appears for less than ~1 render frame (no strobe)."
    why_human: "Live UniProt fetchProxy + DG.TaskBarProgressIndicator + DOM render timing cannot be observed by grep; the new unit test asserts the show/update/hide lifecycle on a synthetic ScatterPlotViewer, not the live counter behavior 13-UAT Test 3 verbatim called out"
  - test: "Volcano visual parity vs CK-omics reference figure (carried forward from prior VERIFICATION human-UAT round)"
    expected: "Volcano matches client's static HTML out of the box — DMD-depleted proteins negative log2FC in blue, DMD-enriched red, threshold lines correct, no manual Comparison flip required (per-row sign normalization handles it)"
    why_human: "Visual parity vs a static HTML reference figure requires the live BP DMD/WT Candidates file rendered in Datagrok"
  - test: "Multi-contrast Candidates import docks a native Filters viewer scoped to Comparison"
    expected: "Filters viewer appears docked RIGHT scoped to the Comparison column on >1 distinct Comparison values; single-contrast skips it"
    why_human: "DG.Viewer.filters() + dockManager.dock() require a live Datagrok shell"
  - test: "UniProt subcellular location fetch produces sane categories + populates the cross-session cache"
    expected: "Known cytoskeletal/membrane proteins (DMD/DYSF/TTN) classified to Cytoplasm/Plasma Membrane; second session loads instantly without new UniProt network calls"
    why_human: "Live UniProt network fetch + userDataStorage persistence require a running instance"
---

# Phase 13: CK-omics Volcano and Enrichment Parity — Re-Verification Report

**Phase Goal:** The Proteomics volcano and enrichment outputs reach feature parity with the
client's CK-omics tool deliverables for the BP DMD/WT engagement — and the recent Wave-5
gap closures (13-09 Option A enrichment re-docking, 13-10 in-volcano busy overlay) close
the two open user-perception complaints from 13-UAT Tests 2 and 3.

**Re-verified:** 2026-06-03 — after 13-09 (enrichment Option A) and 13-10 (volcano busy overlay)
landed. Previous verification (2026-05-17): 6/6 truths verified, status human_needed (5 live
checks pending). Current verification: 6/6 truths verified, plus the two newly-claimed gap
closures verified at the file/test level. Two new human-UAT items added for the perception
side of the closures.

---

## Goal Achievement

### Observable Truths (R1–R6)

| Req | Truth | Status | Evidence (file:line) |
|-----|-------|--------|-----------|
| **R1** | Subcellular-location volcano coloring: SEMTYPE + detectors mirror + 11-category palette + UniProt stream fetch + cross-session cache + GO-CC fallback + reviewed-by-gene fallback + color-by integration in volcano | **PASS** | `src/utils/proteomics-types.ts:7` SEMTYPE.SUBCELLULAR_LOCATION; `detectors.js:39-49` mirror; `src/analysis/subcellular-location.ts:185` ProgressCb, :192 FETCH_CONCURRENCY=6, :198 CACHE_FLUSH_INTERVAL_MS=5000, :219 runWithConcurrency, :289 flushCache, :298 setInterval, :306 / :347 runWithConcurrency calls, :387 final flushCache; `src/viewers/volcano.ts:143` LOCATION_HASH_TAG, :158-206 ensureLocationColumn with short-circuit at :181-187; recomputeVolcano at :562-587 |
| **R2** | Up/down enrichment split with WP source default, side-by-side dot+bar viewers, Phase-9 cross-link preserved | **PASS** | `src/analysis/enrichment.ts:7` openEnrichmentVisualization import; runEnrichmentPipeline issues two gGOSt calls (Up/Down); merged Direction column via append; `src/viewers/enrichment-viewers.ts:142-154` filterByDirection; :215-230 directional split docks 4 viewers via chartHost; :234-238 three wireEnrichmentToVolcano subscriptions preserved |
| **R3** | Candidates per-row sign normalization (canonical preservation; A6 AVG Group Quantity guard) + report-path D-09 declared-contrast default | **PASS** | `src/parsers/spectronaut-candidates-parser.ts:95-181` normalizeCandidatesSign pure (no grok.shell ref — `grep -c grok.shell` = 0); per-row flip on cmpCol.get(i) === reversed only; A6 num/den swap only when both findCol succeed; `src/analysis/differential-expression.ts` getDefaultComparison still default in showDEDialog (carried forward from prior verification — unchanged) |
| **R4** | Multi-contrast Candidates auto-docks a Filters viewer scoped to Comparison | **PASS** | `src/package.ts:137` dockComparisonFilterIfMultiContrast exported; :315 called inside importSpectronautCandidates after addTableView |
| **R5** | WikiPathways g:Profiler source default-on; appears in enrichment results | **PASS** | `src/analysis/enrichment.ts` wpInput default true; 'WP' literal in selectedSources; HUMAN-UAT Test 2 already confirmed 235 WP rows in a live run (per `.planning/phases/13-.../13-HUMAN-UAT.md:32-33`) |
| **R6** | Q-value vs P-value metric toggle synchronizes Y axis, classes, and threshold lines | **PASS** | `src/viewers/volcano.ts:62-73` ensureNegLog10Column parameterized; :89-124 ensureDirectionColumn parameterized; :562-587 recomputeVolcano single-pass synchronized re-run; `src/package.ts:396-506` volcanoOptions dialog drives the toggle |

**Score:** 6/6 must-have R1–R6 truths VERIFIED.

---

### Gap-Closure Verification (13-09 + 13-10)

| Gap | Plan | Status | Static Evidence | Live (human) Evidence |
|-----|------|--------|-----------------|------------------------|
| 13-UAT Test 2 (major): dot/bar charts feel "associated with the enrichment table not the volcano table" | **13-09 Option A** | **PASS (static) / HUMAN-NEEDED (live)** | `src/viewers/enrichment-viewers.ts:176-182` proteinTv resolved via Dart-handle identity; :197 `chartHost = proteinTv ?? tv`; :223-230 all 4 directional dock calls on chartHost; :252-253 single-direction fallback on chartHost; :244 and :260 `if (proteinTv) grok.shell.v = proteinTv` in BOTH branches; `import {createVolcanoPlot}` from './volcano' fully removed (grep returns no matches); `src/tests/enrichment-visualization.ts:201-235, 238-260, 264-284` Tests B/C/D assert delta ≥4/≥2 viewers on proteinTv + countViewersBoundTo(enrichTv, proteinDf)===0 + Dart-handle focus identity; the smart-filter banner still docks on `tv` (`src/viewers/enrichment-viewers.ts:212`) | Routed to HUMAN-UAT item #1 — visual layout / next-to-volcano perception requires live DG.DockManager render |
| 13-UAT Test 3 (minor): first-time Color-by-Location feels like a hang | **13-10 in-volcano overlay** | **PASS (static) / HUMAN-NEEDED (live)** | `src/viewers/volcano.ts:399-428` showVolcanoBusy (idempotent via inline hideVolcanoBusy; centered card; data-volcano-busy marker; zIndex 6); :430-438 updateVolcanoBusy (label + detail mutation); :440-442 hideVolcanoBusy (querySelectorAll remove); :267-277 disposeVolcanoAttachments sweep extended with `'[data-volcano-busy]'`; `src/package.ts:17-18` showVolcanoBusy/updateVolcanoBusy/hideVolcanoBusy imports; :479 willFetchLocation guard; :480 showVolcanoBusy on entry; :491-494 updateVolcanoBusy from inside ProgressCb lambda alongside pi.update; :500-502 hideVolcanoBusy in finally; 13-08 mechanics PRESERVED VERBATIM (FETCH_CONCURRENCY=6, CACHE_FLUSH_INTERVAL_MS=5000, runWithConcurrency, flushCache, LOCATION_HASH_TAG, ProgressCb, fnv1aHex hash all grep-confirmed); `src/tests/volcano.ts:173-200` overlay lifecycle test (attach → update → detach); :202-245 warm-cache short-circuit test (< 50ms, exactly one init-column tick) | Routed to HUMAN-UAT item #2 — live counter cadence + warm-cache no-strobe requires live UniProt fetchProxy timing |

---

### Required Artifacts (Updated)

| Artifact | Provides | Status | Details |
|----------|----------|--------|---------|
| `src/viewers/enrichment-viewers.ts` | Option A protein-view dock + smart-filter banner on enrichment view + focus switch | VERIFIED | 263 lines; no createVolcanoPlot import; chartHost = proteinTv ?? tv pattern; grok.shell.v switch in both branches |
| `src/viewers/volcano.ts` | Overlay trio + disposeVolcanoAttachments sweep | VERIFIED | 588 lines; showVolcanoBusy/updateVolcanoBusy/hideVolcanoBusy exported and additive; sweep selector list extended |
| `src/package.ts` | volcanoOptions wires overlay; enrichment handlers untouched | VERIFIED | volcanoOptions willFetchLocation guard + overlay attach/update/detach; enrichmentAnalysis/enrichmentCharts untouched (both still terminate in openEnrichmentVisualization) |
| `src/analysis/subcellular-location.ts` | 13-08 mechanics LOCKED | VERIFIED | All grep gates pass (FETCH_CONCURRENCY=6 at :192, CACHE_FLUSH_INTERVAL_MS=5000 at :198, runWithConcurrency at :219+:306+:347, flushCache at :289+:387, setInterval at :298) |
| `src/parsers/spectronaut-candidates-parser.ts` | R3 per-row sign + R4 contrast plumbing | VERIFIED | normalizeCandidatesSign pure (no grok.shell reference); proteomics.source/de_complete tags set after parse |
| `src/tests/enrichment-visualization.ts` | Tests B/C/D for Option A | VERIFIED | 13 tests total; B/C invert assertion (viewers on proteinTv, 0 on enrichTv); D asserts focus switch; banner tests unchanged and exercise proteinTv-undefined fallback |
| `src/tests/volcano.ts` | Overlay lifecycle + warm-cache short-circuit | VERIFIED | 9 Volcano-category tests (7 baseline + 2 new); inline patchUserDataStorage helper duplicated from subcellular-location.ts L22 |
| `src/tests/subcellular-location.ts` | 13-08 regression baseline | VERIFIED | 16 tests preserved (13-10 Task 3 grep-gated regression hint matched in code) |

---

### Key Link Verification (Updated)

| From | To | Via | Status | Detail |
|------|----|-----|--------|--------|
| `src/viewers/enrichment-viewers.ts openEnrichmentVisualization` | proteinTv.dockManager.dock | chartHost = proteinTv ?? tv | WIRED | :197; :223-230, :252-253 reference chartHost |
| `src/viewers/enrichment-viewers.ts` | `grok.shell.v = proteinTv` | post-dock focus switch | WIRED | :244 and :260 (both branches guarded by `if (proteinTv)`) |
| `wireEnrichmentToVolcano` subscriptions | proteinDf.selection | unchanged | WIRED | :234-238 (3 subscriptions in directional branch), :255-256 (2 in single-direction fallback); identical behavior to pre-13-09 |
| `src/package.ts volcanoOptions OK handler` | showVolcanoBusy/updateVolcanoBusy/hideVolcanoBusy | overlay lifecycle | WIRED | :17-18 imports; :480 show; :491-494 update inside ProgressCb; :502 hide in finally |
| `ProgressCb (done, total, phase)` ticks | in-volcano overlay text update | volcanoOptions wires lambda to BOTH pi.update + updateVolcanoBusy | WIRED | :483-495 the lambda calls both; the willFetchLocation guard suppresses overlay on significance-only path |
| `ensureLocationColumn short-circuit (13-08)` | overlay teardown without flicker | progress fires once with phase='init-column'; hideVolcanoBusy in finally | WIRED | `src/viewers/volcano.ts:185` short-circuit fires init-column; warm-cache test asserts < 50ms + exactly one tick |

---

### Anti-Patterns Found

| File | Concern | Severity | Impact |
|------|---------|----------|--------|
| `src/viewers/enrichment-viewers.ts:197` | `chartHost = proteinTv ?? tv` fallback can re-route docks onto enrichTv when proteinTv is unresolvable (e.g. via `enrichmentCharts` cmd-line entry when there's no DE-complete protein view) | Info | Documented by 13-09 plan as intentional graceful degradation; the 2 banner tests intentionally exercise this fallback (they don't addTableView(proteinDf)). Not a regression — matches pre-13-09 behavior exactly when fallback fires. |
| `src/viewers/enrichment-viewers.ts:244, :260` | Two identical `if (proteinTv) grok.shell.v = proteinTv` blocks (DRY violation) | Info | Cost: one extra grep find; benefit: each branch reads independently. Not material. |
| `src/package.ts:479` | `willFetchLocation` is always true when colorDim==='location' (no warm-cache short-circuit detection at decision time — the planner deliberately punted the hash-recompute optimization to a future plan) | Info | Documented in 13-10 plan; the overlay flashes for <1 render frame on warm-cache via the synchronous short-circuit in ensureLocationColumn — acceptable per the human-UAT acceptance criterion |
| Pre-existing CR-01/CR-02/CR-03 from prior VERIFICATION | Float32Array + DG.FLOAT_NULL staging; splitGenesByDirection gene-to-multi-row conflict | Warning | Unchanged from 2026-05-17 verification — none touched by 13-09 or 13-10 |

**Debt-marker gate:** `grep -E "TBD\|FIXME\|XXX\|HACK\|PLACEHOLDER"` against `enrichment-viewers.ts`, `volcano.ts`, `package.ts`, `subcellular-location.ts` returned **zero matches**.

---

### Probe Execution / Spot-Checks

- Step 7c: No `scripts/*/tests/probe-*.sh` defined for Phase 13. SKIPPED.
- Step 7b: Static signal-tests only — overlay DOM lifecycle and short-circuit timing are covered by the new tests in `src/tests/volcano.ts:173-200, 202-245`; layout/focus by `src/tests/enrichment-visualization.ts:201-284`. Live perception checks ROUTED TO HUMAN-UAT (items #1 and #2 above).

---

### Structural Concern: Option A → Option B → Option A Oscillation (13-07 → 13-09)

**Pattern observed:**
- **Phase 13-07** chose Option B for the Test 2 gap: co-dock a small `createVolcanoPlot(proteinDf, {title: 'Volcano (linked)'})` on the enrichment view at ratio 0.4. Tests were updated to assert that link was present on enrichTv.
- **13-UAT round 2 (post-13-07)** confirmed the user reported the identical verbatim complaint — Option B did not close the gap.
- **Phase 13-09** reversed direction: removed the Option B co-dock entirely, re-pointed the dot/bar docks to the protein TableView (proteinTv) via `chartHost = proteinTv ?? tv`, and added a `grok.shell.v = proteinTv` focus switch. Tests B and C were inverted (delta ≥ 4/2 on proteinTv, 0 proteinDf-bound viewers on enrichTv); test D was added for the focus contract.

**Why this matters:**
- 13-09 deletes code that 13-07 added in the same phase. The intermediate state (Option B shipped in a real package build) was tested live and falsified. That is a healthy outcome of human-UAT — the SPR/CK-omics-volcano cycle treated the user complaint as falsifying evidence rather than as a noise signal.
- The cost is one wasted execute-phase iteration. The benefit is that 13-09 is now a more conservative implementation (removes Option B viewer; preserves all data-layer wiring; focuses on the user's mental model: the dot/bar charts belong next to "my" volcano).
- The risk going forward: the 2 banner tests in `enrichment-visualization.ts` exercise the `proteinTv === undefined` fallback path because they don't `addTableView(proteinDf)`. That fallback re-routes the dot/bar docks to enrichTv exactly like the pre-13-09 behavior. If a future entry point (e.g. an `enrichmentCharts` re-entry where the protein table has been closed) hits the fallback in production, the user will be back in Option-B-style "charts on the wrong view" territory. This fallback should be either (a) loud (e.g. `grok.shell.warning('No protein view found; docking on enrichment view')`) or (b) blocked (raise an explicit error and stay on the enrichment view rather than pretend Option A is in effect).

**No code change requested here** — this is a structural advisory for whoever drafts Phase 14+ or the next gap-closure cycle. The current code is correct for the primary user path (Proteomics | Analyze | Enrichment from a DE-complete protein view); the fallback path is a known-narrow surface.

---

### Requirements Coverage

| Req | Plans | Status | Evidence Path |
|-----|-------|--------|----------------|
| R1 | 13-01, 13-04, 13-06, 13-08, 13-10 | SATISFIED | LOCATION_KEYWORDS/LOCATION_COLORS/parseSubcellularLocation/getSubcellularLocations/ensureLocationColumn/recomputeVolcano all present and wired; 13-08 bounded-concurrency + cache mechanics preserved verbatim; 13-10 in-volcano busy overlay live |
| R2 | 13-02, 13-07, 13-09 | SATISFIED | splitGenesByDirection; two-gGOSt runEnrichmentPipeline; Direction column merged; 4-viewer side-by-side via chartHost in 13-09 (now on protein view); Phase-9 wireEnrichmentToVolcano cross-link preserved verbatim through both 13-07 and 13-09 |
| R3 | 13-03, 13-05 | SATISFIED | normalizeCandidatesSign pure; per-row flip only; A6 AVG Group Quantity guard; D-09 declared-contrast default in showDEDialog |
| R4 | 13-03 | SATISFIED | dockComparisonFilterIfMultiContrast at `src/package.ts:137`; called inside importSpectronautCandidates at :315 |
| R5 | 13-02 | SATISFIED | wpInput default true; 'WP' in selectedSources; HUMAN-UAT round confirmed 235 WP rows on a live BP DMD/WT run |
| R6 | 13-06 | SATISFIED | ensureNegLog10Column / ensureDirectionColumn parameterized by MetricKind; recomputeVolcano single-pass synchronized; Volcano Options dialog drives the toggle (post-13-10: also drives the in-volcano busy overlay) |

---

### Gaps Summary

**No automated-verifiable gaps.** All 6 R1–R6 must-haves are implemented and wired in the codebase.

**Open human-UAT items** (routed in the frontmatter above):
1. Visual confirmation that dot/bar charts dock next to the user's volcano on the protein view (13-UAT Test 2 closure)
2. Visual confirmation that the in-volcano overlay flashes a live N/80 counter on first-time Color-by-Location and does not strobe on warm-cache toggles (13-UAT Test 3 closure)
3. Carried forward from the 2026-05-17 round: CK-omics volcano visual parity, multi-contrast Filters dock, UniProt fetch + cache.

The two new gap-closure items (1 and 2) are perception-driven by design — the automated tests cover the layout/focus/lifecycle contracts (the smallest assertions that would have caught the verbatim user complaints) but not the visible-to-the-user behavior the user evaluates. This is the canonical "static green, live red" risk pattern: prior verification missed the same complaint because the data-layer assertion (selection.trueCount) was green while the layout was wrong. The Option A → Option B → Option A oscillation above is the lesson learned.

---

_Verified: 2026-06-03T00:00:00Z_
_Re-verifier: Claude (gsd-verifier)_
