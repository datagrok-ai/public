---
phase: 9
slug: enrichment-visualization-volcano-integration
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-06
---

# Phase 9 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | @datagrok-libraries/test (Puppeteer-based) |
| **Config file** | packages/Proteomics/src/package-test.ts |
| **Quick run command** | `grok test --host localhost --category "Enrichment Visualization"` |
| **Full suite command** | `grok test --host localhost` |
| **Estimated runtime** | ~30 seconds |

---

## Sampling Rate

- **After every task commit:** Run `grok test --host localhost --category "Enrichment Visualization"`
- **After every plan wave:** Run `grok test --host localhost`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 30 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 09-01-01 | 01 | 1 | VIZ-02, VIZ-03 | unit | `grok test --host localhost --test "top-N filter"` | No - W0 | pending |
| 09-01-02 | 01 | 1 | VIZ-02 | unit | `grok test --host localhost --test "dot plot"` | No - W0 | pending |
| 09-01-03 | 01 | 1 | VIZ-03 | unit | `grok test --host localhost --test "bar chart"` | No - W0 | pending |
| 09-01-04 | 01 | 1 | VIZ-02, VIZ-03 | unit | `grok test --host localhost --test "enrichment dashboard"` | No - W0 | pending |
| 09-02-01 | 02 | 2 | ENRICH-04 | unit | `grok test --host localhost --test "cross-DF selection"` | No - W0 | pending |
| 09-02-02 | 02 | 2 | ENRICH-04 | unit | `grok test --host localhost --test "volcano highlight"` | No - W0 | pending |

*Status: pending / green / red / flaky*

---

## Wave 0 Requirements

- [ ] `src/tests/enrichment-visualization.ts` — test stubs for VIZ-02, VIZ-03, ENRICH-04, top-N filtering
- [ ] Add `import './tests/enrichment-visualization';` to `package-test.ts`
- [ ] No new framework install needed — @datagrok-libraries/test already configured

*Existing infrastructure covers framework requirements.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Dot plot visual layout (bubble sizes, color gradient) | VIZ-02 | Visual rendering quality | Open enrichment viz, verify bubbles scale with Gene Count, color gradient maps to significance |
| Bar chart sort order and orientation | VIZ-03 | Visual rendering quality | Open enrichment viz, verify bars sorted by significance descending, horizontal layout |
| Volcano highlight visual feedback | ENRICH-04 | Cross-view UI interaction | Click enrichment row, verify blue selection dots appear on volcano with gene name labels |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 30s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
