---
phase: 14
slug: ck-omics-analyst-experience-enhancements
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-05-28
---

# Phase 14 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution. Datagrok TS package — `grok test` is the host-server suite; tests live in `src/tests/` and re-export through `src/package-test.ts`.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | `@datagrok-libraries/utils/src/test` (Datagrok in-host `grok test`) — runs the registered `tests` array against a live server (localhost or dev). No vitest/jest. |
| **Config file** | `package.json` (test scripts), `webpack.config.js` (test bundle entry `src/package-test.ts`) |
| **Quick run command** | `grok test --skip-build --catch-unhandled=false --category "Proteomics: Phase 14"` (re-build after any change in `src/tests/` — `--skip-build` reuses stale bundle and silently drops new tests; see project memory `feedback_grok_test_skipbuild_stale`) |
| **Full suite command** | `grok test --catch-unhandled=false` (rebuilds, runs every registered category) |
| **Estimated runtime** | ~60s quick / ~300s full (Datagrok host roundtrip dominates) |

The package targets a running Datagrok server. The user's working server is the localhost in `~/.grok/config.yaml`. Tests that call live endpoints (Ensembl REST, UniProt) are gated behind `category: 'Proteomics: External'` and skipped in CI; deterministic unit-style tests (parsers, helpers, smart-filter port) run unconditionally.

---

## Sampling Rate

- **After every task commit:** Run the quick command above scoped to the touched plan's category (`--category "Proteomics: 14-01"`, etc.). Resampling cadence per `references/tdd.md` heuristics — never 3 consecutive code tasks without an automated check.
- **After every plan wave:** Run the quick command across all Phase-14 categories (`--category "Proteomics: Phase 14"`).
- **Before `/gsd:verify-work`:** Full suite must be green; volcano + UniProt panel + enrichment + parser categories all pass.
- **Max feedback latency:** ~60 seconds for the quick command.

---

## Per-Task Verification Map

The planner will fill exact `{N}-{plan}-{task}` IDs once PLAN.md files exist. The skeleton below names the verifiable seams the research identified (`14-RESEARCH.md` validation architecture section) so each plan can attach `<automated>` blocks to them.

| Seam | Plan (target) | Wave (target) | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|------|---------------|---------------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| `geneLabelResolver.detectSpecies(id)` | 14-01 (R1) | 1 | R1 | — | unknown prefix → returns `null`, downstream keeps raw ID + `†` | unit (no host) | `grok test --category "Proteomics: 14-01"` | ❌ W0 (`src/tests/gene-label-resolver.ts` new) | ⬜ pending |
| `applySmartPathwayFilter(enrichDf, opts)` (port of CK-omics `apply_smart_pathway_filtering`) | 14-05 (R5) | 2 | R5 | — | drops generic parents when a child is in the result set; caps non-GO:BP sources at combined top-N; default-on with dialog opt-out | unit (fixture DataFrame) | `grok test --category "Proteomics: 14-05"` | ❌ W0 (`src/tests/smart-pathway-filter.ts` new) | ⬜ pending |
| Direction-column strings + color hex | 14-02 (D-04/G1) | 2 | R2/G1 | — | DataFrame contains `"Enriched in <g1>"` / `"Enriched in <g2>"` / `"Not significant"`; volcano color column carries `0xFFFF00FF` / `0xFF00FFFF` / `0xFFAAAAAA` | unit (DataFrame assertion) | `grok test --category "Proteomics: 14-02"` | ❌ W0 (extend `src/tests/volcano.ts`) | ⬜ pending |
| Live-counter recompute on `df.filter` / `df.selection` / viewer-property | 14-02 (R2) | 2 | R2 | — | After `df.filter.copyFrom(...)`, counter reads new totals within next animation frame; subscriptions disposed on viewer detach | unit + integration | `grok test --category "Proteomics: 14-02"` | ❌ W0 (extend `src/tests/volcano.ts`) | ⬜ pending |
| Filters viewer scope excludes `Flags` (G4) and respects `IFiltersSettings.filters` | 14-03 (D-05/G4) | 3 | G4 | — | docked viewer's `getOptions().look.columnNames` reflects exactly `[Comparison, Display Name, Protein ID]`; `Flags` absent | integration (in-host) | `grok test --category "Proteomics: 14-03"` | 14-03 T3 (`filtersScopingNoFlags`, `filtersScopingSingleContrast`, `filtersScopingFallbackToProteinIdWhenNoDisplayName` in `src/tests/spectronaut-candidates-parser.ts`) | ⬜ pending live run |
| UniProt panel per-group bar chart (D-11/R3) | 14-04 (R3) | 3 | R3 | — | panel root contains 2 (or N) bars with mean+SD text; updates on cell-row change | integration (panel mount) | `grok test --category "Proteomics: 14-04"` | ❌ W0 (extend `src/tests/uniprot-panel.ts` or create) | ⬜ pending |
| Group-mean correlation viewer derived columns + corr title (D-12/R4) | 14-04 (R4) | 3 | R4 | — | `Numerator Mean` + `Denominator Mean` columns present and unique; title annotation contains `r=` (Pearson) and `ρ=` (Spearman); no duplicate on re-run | unit (viewer factory) | `grok test --category "Proteomics: 14-04"` | ❌ W0 (`src/tests/group-mean-correlation.ts` new) | ⬜ pending |
| Ensembl resolver cache mirrors Phase-13 `__schema_v` pattern (D-09) | 14-01 (R1) | 1 | R1 | — | After resolve, `grok.userSettings` contains `proteomics.ensembl.cache.v{N}` keys; bumping `__schema_v` invalidates | unit + manual | `grok test --category "Proteomics: 14-01-cache"` | ❌ W0 (extend gene-label-resolver tests) | ⬜ pending |
| Volcano Options dialog preload reads `viewer.getOptions().look.*` + `df.tag('proteomics.volcano_metric')` (G2) | 14-02 (G2) | 2 | G2 | — | Reopening dialog after changing color column shows the now-applied column; OK→reopen idempotent | integration (dialog reopen) | `grok test --category "Proteomics: 14-02-dialog"` | ❌ W0 (extend `src/tests/volcano.ts`) | ⬜ pending |
| Color → Location switch progress + cache hit (G3) | 14-02 (G3) | 2 | G3 | — | First switch shows `DG.TaskBarProgressIndicator`; second switch hits cache (assert via spy that no network fetch fires) | integration | `grok test --category "Proteomics: 14-02-color"` | ❌ W0 (extend `src/tests/volcano.ts`) | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky · ❌ W0 = test file does not yet exist; Wave 0 of the relevant plan creates it.*

The planner will refine the `Plan (target)` and `Wave (target)` columns when assigning tasks. The intent is one row per verifiable seam from `14-RESEARCH.md`.

---

## Wave 0 Requirements

- [ ] `src/tests/gene-label-resolver.ts` — fixture IDs for ENSRNOG / ENSMUSG / ENSG / ENSDARG / MGP_ / LOC / RGD / AABR; covers R1 detection regex, three-level fallback (external_name → display_name → description), `*` and `†` marker rules, duplicate-description disambiguation
- [ ] `src/tests/smart-pathway-filter.ts` — fixture g:Profiler-shaped DataFrame; covers R5 generic-parent drop heuristic, non-GO:BP combined cap, default-on / dialog-off paths
- [ ] `src/tests/group-mean-correlation.ts` — covers D-12 derived columns (ensureFreshFloat re-run safety), Pearson/Spearman computation (inline; `@datagrok-libraries/statistics` only exports Kendall), title annotation string
- [ ] Extend `src/tests/volcano.ts` with: D-04 direction strings + color hex (G1), live-counter recompute on `df.filter`/`df.selection`/property change (R2), Volcano Options dialog preload (G2), Color → Location progress + cache (G3)
- [ ] Extend or add `src/tests/uniprot-panel.ts` with per-group bar chart mount (D-11/R3)
- [ ] Extend `src/tests/spectronaut-candidates-parser.ts` (or sibling) with multi-contrast Filters scoping assertion (G4)
- [ ] Confirm test categories register in `src/package-test.ts` (`tests` array entry per category — see existing `Proteomics: Phase 13` pattern)

*Framework install:* none — `grok test` is host-side and already wired in this package (see existing `src/tests/*.ts` and `package.json` scripts).

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Volcano visual parity vs CK-omics figure | G1 | Pixel/typographic comparison against `~/Downloads/ck/DMD_vs_WT/volcano_plots/` — title position, axis label rewriting, default top-N labels, magenta/cyan/gray legend with counts, search box dock arrangement | Open Spectronaut Candidates sample, run pipeline, compare side-by-side with reference PNG; verify magenta=numerator, cyan=denominator, ≥15 top-hit labels visible |
| Live counter overlay placement and update fluidity | R2 | Subjective UX — float anchor, font sizing, perceived lag | Open volcano on a >5k-row DataFrame, toggle filter categories on/off, confirm counters update within one animation frame and don't shift layout |
| Per-group quantity bar-chart readability | R3/D-11 | Visual mean±SD chart fit inside panel column width | Click 5 different proteins in the volcano, confirm bars + SD whiskers render legibly without overflow |
| Smart-filter result quality on real client dataset | R5/D-13 | Domain-expert judgement (analyst confirms the dropped parents are indeed generic) | Run enrichment on the BP DMD-vs-WT sample twice (default-on, then dialog-off), confirm default-on result is cleaner without losing analyst-relevant terms |
| Ensembl resolution on a real ENSRNOG/ENSMUSG mix | R1 | Network roundtrip + label readability | Import a rat or mouse Spectronaut Candidates file with ≥50 ENS prefixed IDs, verify Display Name column shows resolved names with `*`/`†` markers and Source ID column retains the raw ID for hover |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references (5 new/extended test files above)
- [ ] No watch-mode flags
- [ ] Feedback latency < 60s (quick run scoped to single plan category)
- [ ] `nyquist_compliant: true` set in frontmatter once planner pins task IDs and the seams table closes

**Approval:** pending
</content>
</invoke>
