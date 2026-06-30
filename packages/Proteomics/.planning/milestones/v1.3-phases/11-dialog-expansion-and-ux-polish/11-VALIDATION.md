---
phase: 11
slug: dialog-expansion-and-ux-polish
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-07
---

# Phase 11 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Datagrok Puppeteer-based test framework |
| **Config file** | packages/Proteomics/webpack.config.js |
| **Quick run command** | `grok test --test "Dialog" --host localhost` |
| **Full suite command** | `grok test --host localhost` |
| **Estimated runtime** | ~60 seconds |

---

## Sampling Rate

- **After every task commit:** Run `grok test --test "Dialog" --host localhost`
- **After every plan wave:** Run `grok test --host localhost`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 60 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 11-01-01 | 01 | 1 | NORM-03 | manual-only | N/A — requires visual dialog interaction | N/A | ⬜ pending |
| 11-01-02 | 01 | 1 | NORM-04 | unit | Test `proteomics.preNormalized` tag detection | Wave 0 | ⬜ pending |
| 11-02-01 | 02 | 1 | IMP-03 | manual-only | N/A — requires DOM interaction | N/A | ⬜ pending |
| 11-02-02 | 02 | 1 | IMP-04 | unit | Test filter counting logic with mock DataFrame | Wave 0 | ⬜ pending |
| 11-02-03 | 02 | 1 | DE-01 | unit | Test pair generation from GroupAssignment | Wave 0 | ⬜ pending |
| 11-02-04 | 02 | 1 | DE-02 | manual-only | N/A — requires DOM interaction | N/A | ⬜ pending |
| 11-03-01 | 03 | 2 | UX-01 | smoke | Verify viewer creation includes `title` option | Wave 0 | ⬜ pending |
| 11-03-02 | 03 | 2 | UX-02 | unit | Test all parser output df.name values | Wave 0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] Unit test for `preNormalized` tag detection logic (NORM-04)
- [ ] Unit test for valid-values filter counting logic (IMP-04)
- [ ] Unit test for comparison pair generation from GroupAssignment (DE-01)
- [ ] Smoke test for viewer title option passthrough (UX-01)
- [ ] Unit test for df.name derivation from filename across all parsers (UX-02)

*Existing infrastructure covers dialog interaction patterns (NORM-03, IMP-03, DE-02) which are manual-only.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Normalization dialog shows method selector + box plot preview | NORM-03 | Requires DOM rendering in running Datagrok | Open dialog -> see method dropdown (Median/Quantile/VSN) -> box plot shows distributions -> change method -> box plot updates |
| Imputation dialog conditional params per method | IMP-03 | Requires DOM interaction | Open dialog -> select MinProb -> see downshift/width -> select kNN -> see k param, MinProb hides -> select Zero -> no extra params |
| DE dialog method-specific params shown/hidden | DE-02 | Requires DOM interaction | Open dialog -> select DEqMS -> peptide column appears -> select limma -> peptide column hides |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 60s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
