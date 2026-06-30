---
phase: 10
slug: spectronaut-parser-and-core-algorithms
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-07
---

# Phase 10 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | @datagrok-libraries/test (Puppeteer-based) |
| **Config file** | webpack.config.js (package-test.ts entry point) |
| **Quick run command** | `grok test --host localhost --category "Spectronaut"` |
| **Full suite command** | `grok test --host localhost` |
| **Estimated runtime** | ~30 seconds |

---

## Sampling Rate

- **After every task commit:** Run `grok test --host localhost --category "Spectronaut"`
- **After every plan wave:** Run `grok test --host localhost`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 30 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 10-01-01 | 01 | 1 | SPEC-01 | unit | `grok test --host localhost --test "Spectronaut: pivot"` | ❌ W0 | ⬜ pending |
| 10-01-02 | 01 | 1 | SPEC-02 | unit | `grok test --host localhost --test "Spectronaut: filter"` | ❌ W0 | ⬜ pending |
| 10-01-03 | 01 | 1 | SPEC-03 | unit | `grok test --host localhost --test "Spectronaut: groups"` | ❌ W0 | ⬜ pending |
| 10-01-04 | 01 | 1 | SPEC-04 | unit | `grok test --host localhost --test "Spectronaut: preNormalized"` | ❌ W0 | ⬜ pending |
| 10-02-01 | 02 | 1 | NORM-01 | unit | `grok test --host localhost --test "quantile"` | ❌ W0 | ⬜ pending |
| 10-02-02 | 02 | 1 | NORM-02 | integration | `grok test --host localhost --test "VSN"` | ❌ W0 | ⬜ pending |
| 10-02-03 | 02 | 1 | IMP-01 | unit | `grok test --host localhost --test "kNN"` | ❌ W0 | ⬜ pending |
| 10-02-04 | 02 | 1 | IMP-02 | unit | `grok test --host localhost --test "impute"` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `src/tests/spectronaut-parser.ts` — stubs for SPEC-01 through SPEC-04 with synthetic long-format data
- [ ] Tests for NORM-01, NORM-02, IMP-01, IMP-02 in `src/tests/analysis.ts` (extend existing file)
- [ ] Test helper: `makeLongFormatTsv()` utility for building Spectronaut-format test data
- [ ] Register new test file in `src/package-test.ts`

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| VSN R fallback | NORM-02 | Requires R runtime unavailable | Verify console warns and falls back to quantile |
| kNN progress indicator | IMP-01 | Visual UI feedback | Observe progress bar during imputation |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 30s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
