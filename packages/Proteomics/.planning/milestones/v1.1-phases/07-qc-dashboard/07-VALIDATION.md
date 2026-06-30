---
phase: 7
slug: qc-dashboard
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-06
---

# Phase 7 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Datagrok test framework (@datagrok-libraries/test) |
| **Config file** | webpack.config.js (test entry: package-test.ts) |
| **Quick run command** | `grok test --test "QC" --host localhost` |
| **Full suite command** | `grok test --host localhost` |
| **Estimated runtime** | ~30 seconds |

---

## Sampling Rate

- **After every task commit:** Run `grok test --test "QC" --host localhost`
- **After every plan wave:** Run `grok test --host localhost`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 30 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 07-01-01 | 01 | 1 | QC-01 | integration | `grok test --test "QC Dashboard" --host localhost` | ❌ W0 | ⬜ pending |
| 07-01-02 | 01 | 1 | QC-02 | unit | `grok test --test "MA computation" --host localhost` | ❌ W0 | ⬜ pending |
| 07-01-03 | 01 | 1 | QC-03 | unit | `grok test --test "Missing values" --host localhost` | ❌ W0 | ⬜ pending |
| 07-01-04 | 01 | 1 | QC-04 | integration | `grok test --test "Correlation" --host localhost` | ❌ W0 | ⬜ pending |
| 07-01-05 | 01 | 1 | QC-05 | unit | `grok test --test "Box plot" --host localhost` | ❌ W0 | ⬜ pending |
| 07-01-06 | 01 | 1 | QC-06 | unit | `grok test --test "CV computation" --host localhost` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `src/tests/qc-dashboard.ts` — test stubs for QC-01 through QC-06
- [ ] Register test file in `src/package-test.ts`

*Existing test infrastructure (package-test.ts, webpack test config) covers framework setup.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Dashboard layout proportions look correct | QC-01 | Visual assessment of dock panel sizes | Open QC dashboard, verify all 5 viewers visible and reasonably sized |
| Missingness heatmap color scheme readable | QC-03 | Visual quality assessment | Open dashboard, verify white=missing/dark=present is clear |
| Cross-viewer selection highlights correctly | QC-01 | Requires visual confirmation of highlight propagation | Select proteins in MA plot, verify CV plot highlights same rows |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 30s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
