---
phase: 6
slug: generic-matrix-parser
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-05
---

# Phase 6 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | @datagrok-libraries/test (Puppeteer-based) |
| **Config file** | packages/Proteomics/src/package-test.ts |
| **Quick run command** | `grok test --category "Generic Parser" --host localhost` |
| **Full suite command** | `grok test --host localhost` |
| **Estimated runtime** | ~30 seconds |

---

## Sampling Rate

- **After every task commit:** Run `grok test --category "Generic Parser" --host localhost`
- **After every plan wave:** Run `grok test --host localhost`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 30 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 06-01-01 | 01 | 1 | IMPORT-01 | unit | `grok test --test "Generic: parses" --host localhost` | ❌ W0 | ⬜ pending |
| 06-01-02 | 01 | 1 | IMPORT-02 | unit | `grok test --test "Generic: assigns" --host localhost` | ❌ W0 | ⬜ pending |
| 06-01-03 | 01 | 1 | IMPORT-03 | unit | `grok test --test "Generic: auto" --host localhost` | ❌ W0 | ⬜ pending |
| 06-01-04 | 01 | 1 | IMPORT-04 | manual | Manual UI verification | N/A | ⬜ pending |
| 06-01-05 | 01 | 1 | IMPORT-05 | unit | `grok test --test "Generic: log2" --host localhost` | ❌ W0 | ⬜ pending |
| 06-01-06 | 01 | 1 | IMPORT-06 | integration | `grok test --test "Generic: downstream" --host localhost` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `src/tests/generic-parser.ts` — unit tests for generic parser logic (IMPORT-01 through IMPORT-06)
- [ ] Import `./tests/generic-parser` in `package-test.ts`
- [ ] Test helper: `makeGenericCsv()` function for building test CSV strings with various column layouts

*Existing test infrastructure (package-test.ts, parsers.ts, analysis.ts) covers framework setup.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Preview grid updates live as columns change | IMPORT-04 | Dialog interaction required | 1. Open Generic Matrix import 2. Load a CSV 3. Change column selections 4. Verify preview updates |
| Dialog layout: all fields visible at once | UX | Visual verification | 1. Open dialog 2. Verify all inputs visible without scrolling |
| Dataset input widget works correctly | IMPORT-01 | File picker requires OS interaction | 1. Click file picker 2. Browse and select CSV 3. Verify columns populate |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 30s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
