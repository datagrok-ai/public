---
phase: 16
slug: sample-level-spc-tracking
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-06-08
---

# Phase 16 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | grok test (Datagrok in-platform; tests under `src/tests/`) |
| **Config file** | `package.json` test scripts + `src/package-test.ts` |
| **Quick run command** | `grok test --category "SPC" --skip-build` (after first build) |
| **Full suite command** | `npm run build && grok test --category "SPC"` |
| **Estimated runtime** | ~60 seconds full suite (math is pure TS; AppData I/O is local) |

> NOTE: `feedback_grok_test_skipbuild_stale` — `--skip-build` reuses the stale test bundle when test files change. After adding new tests to `src/tests/spc.ts`, run a full `npm run build` before `grok test` to avoid null/0 results.

---

## Sampling Rate

- **After every task commit:** Run `grok test --category "SPC" --skip-build` (assumes prior build)
- **After every plan wave:** Run `npm run build && grok test --category "SPC"`
- **Before `/gsd:verify-work`:** Full suite green + manual dashboard click-through with the locked-baseline fixture
- **Max feedback latency:** ~60 seconds

---

## Per-Task Verification Map

This is the source RESEARCH.md derives — see `16-RESEARCH.md` § "Validation Architecture" for the
full test inventory mapped to SPC-01..SPC-08 + Pitfalls 5/6/7/8. The planner fills the per-task
rows once PLAN.md task IDs are minted.

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 16-XX-YY | XX | N | SPC-XX | T-16-XX / — | {expected behavior or "N/A"} | unit / integration / manual | `grok test --test "{TestName}"` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

> Planner: derive rows directly from RESEARCH.md § "Validation Architecture" — ~25 test ids mapped one-to-one to SPC-01..SPC-08 plus idempotency + belt-and-braces survival. Filter to one row per executable test in `src/tests/spc.ts` (or `src/tests/spc-storage.ts` / `src/tests/spc-dashboard.ts` if the planner splits the test file).

---

## Wave 0 Requirements

- [ ] `src/tests/spc.ts` — RED stubs for SPC-01..SPC-08 metric computations (median, missing %, control-corr, protein count) + Nelson 1+5 rule evaluations + baseline iterative outlier removal (cap=2)
- [ ] `src/tests/spc-storage.ts` (or merged into `spc.ts`) — RED stubs for AppData runs.csv append/read + baseline-<instrument_id>.json round-trip
- [ ] Spectronaut PG-report fixture under `src/tests/fixtures/` — minimal report with `R.RunDate` + `R.InstrumentMethod` columns to drive the parser-seed assertion (SPC-06)
- [ ] In-control + out-of-control synthetic numeric series (small inline arrays in test file) — drive Nelson 1+5 + iterative outlier removal RED stubs

*If existing infrastructure (`src/tests/` + `src/package-test.ts` registry) does not need new framework install, omit a framework-install task.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| `DG.Viewer.lineChart` formula lines (UCL/CL/LCL) render with dashed-vs-solid style on a datetime X-axis | SPC-02 | Renderer behavior is platform-side; can't be asserted from a Datagrok package unit test reliably | Open Visualize → SPC Dashboard..., pick instrument with ≥7 baseline runs, confirm UCL/LCL render as dashed horizontal lines and CL renders as solid (manual visual check + screenshot for VERIFICATION.md) |
| Drill-down opens source project on flagged-point click | SPC-07 | Requires a published project state that's user-environment-specific | With a published run in your AppData + a baseline that flags it, click the flagged point on the I-chart; confirm the project opens via grok.dapi.projects.find(id).open(). If source_project_id is null, confirm the fallback toast surfaces. |
| Empty-state banner + Define-baseline modal copy is biologist-readable | SPC-05 | Banned-words audit + biologist-in-the-room readability check carries Phase 15 audience pin; no automation possible | Pre-demo dress-rehearsal with biologist user; confirm zero `DataFrame` / `tag` / `semType` / `viewer factory` jargon in dashboard, modal, drill-down toast, and Pareto labels |
| Pareto bar order (P2 SPC-08) | SPC-08 | DG.Viewer.barChart sort behavior is platform-side | Run Compute SPC Status on ≥3 runs with multiple rule trips; open Pareto panel; confirm bars sort descending by trip count |

---

## Validation Sign-Off

- [ ] All planner-emitted tasks have either `<automated>` verify command or a documented Wave 0 dependency
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers SPC-01..SPC-08 RED stubs
- [ ] No watch-mode flags (`--watch`, `--gui` without `--skip-build`)
- [ ] Feedback latency < 60s for quick run, < 90s for full suite
- [ ] Manual verifications limited to renderer behavior + drill-down + audience-readability + Pareto sort (4 items above)
- [ ] `nyquist_compliant: true` set in frontmatter once planner rows land

**Approval:** pending
