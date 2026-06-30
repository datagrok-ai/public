---
phase: 13
slug: ck-omics-volcano-and-enrichment-parity
status: draft
nyquist_compliant: true
wave_0_complete: false
created: 2026-05-16
---

# Phase 13 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | `@datagrok-libraries/test` (grok test) |
| **Config file** | none — tests live under `src/tests/`, registered via `src/package-test.ts` |
| **Quick run command** | `grok test --category "<phase category>"` |
| **Full suite command** | `grok test` (requires running Datagrok instance) |
| **Estimated runtime** | ~120 seconds (full suite incl. webpack rebuild) |

---

## Sampling Rate

- **After every task commit:** Run `grok test --category "<phase category>"` (rebuild first — never `--skip-build` after touching test files, memory `feedback_grok_test_skipbuild_stale`)
- **After every plan wave:** Run `grok test` (full suite)
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** ~120 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 13-01-01 | 01 | 1 | R1 | T-13-01-02 | detectors.js literal matches SEMTYPE verbatim (grep gate) | unit (grep) | `grep -q "SUBCELLULAR_LOCATION: 'Proteomics-SubcellularLocation'" src/utils/proteomics-types.ts && grep -q "Proteomics-SubcellularLocation" detectors.js && echo OK` | ✅ | ⬜ pending |
| 13-01-02 | 01 | 1 | R3 | T-13-01-01 | Wave-0 findings file present with A2/A3/A4 headings (no secrets) | unit (file gate) | `F=.planning/phases/13-ck-omics-volcano-and-enrichment-parity/13-WAVE0-FINDINGS.md; test -f "$F" && grep -q A2 "$F" && grep -q A3 "$F" && grep -q A4 "$F" && echo OK` | ❌ W0 | ⬜ pending |
| 13-01-03 | 01 | 1 | R1, R6 | — | test scaffolds discoverable; no behavior surface | unit (grep) | `grep -q "category('SubcellularLocation'" src/tests/subcellular-location.ts && grep -q "category('Volcano'" src/tests/volcano.ts && grep -q "tests/subcellular-location" src/package-test.ts && grep -q "tests/volcano" src/package-test.ts && echo OK` | ❌ W0 | ⬜ pending |
| 13-02-01 | 02 | 1 | R2, R5 | T-13-02-01 / T-13-02-02 | g:Profiler body via JSON.stringify (no concat); empty-direction call skipped | unit + integration | `cd /Users/edjaeger/datagrok/src/public/.claude/worktrees/proteomics/packages/Proteomics && npm run build && grok test --category Enrichment` | ⚠️ extend | ⬜ pending |
| 13-02-02 | 02 | 1 | R2 | T-13-02-03 | malformed g:Profiler JSON tolerated (intersections optional) | unit + integration | `cd /Users/edjaeger/datagrok/src/public/.claude/worktrees/proteomics/packages/Proteomics && npm run build && grok test --category Enrichment` | ⚠️ extend | ⬜ pending |
| 13-03-01 | 03 | 2 | R3 | T-13-03-01 / T-13-03-02 | unparseable comparison → no flip (never invert on ambiguity); absent qty cols guarded | unit (pure) | `cd /Users/edjaeger/datagrok/src/public/.claude/worktrees/proteomics/packages/Proteomics && npm run build && grok test --category SpectronautCandidates` | ⚠️ extend | ⬜ pending |
| 13-03-02 | 03 | 2 | R4 | T-13-03-03 | >1 distinct Comparison docks Filter viewer; single does not | integration (shell) | `cd /Users/edjaeger/datagrok/src/public/.claude/worktrees/proteomics/packages/Proteomics && npm run build && grok test --category "SpectronautCandidates E2E"` | ⚠️ extend | ⬜ pending |
| 13-04-01 | 04 | 2 | R1 | T-13-04-01 | positional parse length-guarded; missing field → 'Unknown'; per-batch catch | unit (pure) | `cd /Users/edjaeger/datagrok/src/public/.claude/worktrees/proteomics/packages/Proteomics && npm run build && grok test --category SubcellularLocation` | ❌ W0 | ⬜ pending |
| 13-04-02 | 04 | 2 | R1 | T-13-04-02 / T-13-04-03 / T-13-04-04 | mandatory fetchProxy (no raw fetch); cache validated vs 12-name enum + __v; chunked OR-query under URL limit | unit (pure) + integration | `cd /Users/edjaeger/datagrok/src/public/.claude/worktrees/proteomics/packages/Proteomics && npm run build && grok test --category SubcellularLocation && grep -q "subcellular-location" src/panels/uniprot-panel.ts && echo OK` | ❌ W0 | ⬜ pending |
| 13-05-01 | 05 | 2 | R3 (D-09) | T-13-05-02 / T-13-05-03 | DE default = declared contrast (not alphabetical); edit isolated to differential-expression.ts | unit/integration | `cd /Users/edjaeger/datagrok/src/public/.claude/worktrees/proteomics/packages/Proteomics && npm run build && grok test --category Analysis` | ⚠️ extend | ⬜ pending |
| 13-06-01 | 06 | 3 | R1, R6 | T-13-06-01 / T-13-06-02 / T-13-06-03 | metric re-init in place (stable binding); p-value-absent guarded; threshold lines replace not stack | unit (pure) + integration | `cd /Users/edjaeger/datagrok/src/public/.claude/worktrees/proteomics/packages/Proteomics && npm run build && grok test --category Volcano` | ❌ W0 | ⬜ pending |
| 13-06-02 | 06 | 3 | R6 | T-13-06-02 | Volcano Options toggle drives recomputeVolcano; p-value choice guarded when column absent | integration (shell) | `cd /Users/edjaeger/datagrok/src/public/.claude/worktrees/proteomics/packages/Proteomics && npm run build && grok test --category Volcano` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

> Derived from `13-RESEARCH.md` § "Validation Architecture" (per-requirement test mapping
> R1–R6 + Wave-0 verification gaps A2/A3) and transcribed task-by-task from
> `13-01-PLAN.md`…`13-06-PLAN.md` `<verify><automated>` blocks. The gsd-planner
> populated this map from the generated PLAN.md tasks; the Nyquist auditor reconciles it
> during execute/verify. "File Exists" legend: ✅ = test/gate ready now · ⚠️ extend =
> existing test file extended by the plan · ❌ W0 = test scaffold created in Wave-0
> (13-01) then filled by the dependent plan.

---

## Wave 0 Requirements

- [ ] Verify real BP DMD/WT Candidates `AVG Log2 Ratio` sign vs declared comparison string (Assumption A2 — gates R3 flip logic; unconditional inversion would re-introduce the documented mirror defect)
- [ ] Trace exact file/function where report-DE group order is set (Assumption A3 — D-09 fix site; resolved to the `differential-expression.ts` `comparisonInput` default per RESEARCH/PATTERNS, Wave-0 confirms the exact line)
- [ ] Confirm `grok.dapi.userDataStorage` capacity for an ~8329-entry accession→location JSON map (Assumption A4 — has IndexedDB fallback)
- [ ] Test stubs for R1–R6 under `src/tests/` if not covered by existing infrastructure (`src/tests/subcellular-location.ts`, `src/tests/volcano.ts` created in 13-01 Task 3)

*Wave-0 work lives in 13-01 Task 2 (A2/A3/A4 → `13-WAVE0-FINDINGS.md`) and 13-01 Task 3 (test scaffolds). `wave_0_complete` flips to `true` after 13-01 executes.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Volcano figure visual parity vs CK-omics static HTML (colour palette, axes, threshold lines) | R1, R6 | Pixel/visual parity against an external client deliverable is not assertable in unit tests | Import BP DMD/WT Candidates file, toggle Subcellular Location colour + Q/P metric, compare against `~/Downloads/ck/DMD_vs_WT/volcano_plots/` static HTML |
| Up/Down enrichment scientific result parity vs CK-omics directional deliverables | R2, R5 | Requires domain judgement on enrichment term concordance | Run enrichment on the same gene sets, compare top terms (incl. WikiPathways) against CK-omics `*_up_*` / `*_down_*` outputs |

*If none: "All phase behaviors have automated verification."*

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify or Wave 0 dependencies
- [x] Sampling continuity: no 3 consecutive tasks without automated verify
- [x] Wave 0 covers all MISSING references (A2, A3, A4)
- [x] No watch-mode flags
- [x] Feedback latency < ~120s
- [x] `nyquist_compliant: true` set in frontmatter

**Approval:** map populated; `wave_0_complete` remains false until 13-01 executes.
