---
phase: 8
slug: gene-id-mapping-enrichment-analysis
status: draft
nyquist_compliant: true
wave_0_complete: false
created: 2026-03-06
---

# Phase 8 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | @datagrok-libraries/test + tsc (Datagrok package test pattern) |
| **Config file** | packages/Proteomics/tsconfig.json |
| **Quick run command** | `cd packages/Proteomics && npx tsc --noEmit --skipLibCheck` |
| **Full suite command** | `cd packages/Proteomics && npm run build` |
| **Estimated runtime** | ~15 seconds |

---

## Sampling Rate

- **After every task commit:** Run `npx tsc --noEmit --skipLibCheck`
- **After every plan wave:** Run `npm run build`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 15 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 08-01-T1 | 01 | 1 | MAP-01, MAP-02, ENRICH-01, ENRICH-02, ENRICH-03, VIZ-01 | unit | `npx tsc --noEmit src/analysis/enrichment.ts --skipLibCheck` | ❌ W0 | ⬜ pending |
| 08-01-T2 | 01 | 1 | MAP-01, MAP-02, ENRICH-01, ENRICH-02, ENRICH-03 | unit | `npx tsc --noEmit src/tests/enrichment.ts --skipLibCheck` | ❌ W0 | ⬜ pending |
| 08-02-T1 | 02 | 2 | MAP-01, MAP-02, ENRICH-01, ENRICH-02, ENRICH-03, VIZ-01 | integration | `npm run build` | ✅ | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `src/tests/enrichment.ts` — test stubs for enrichment pipeline (MAP-01, MAP-02, ENRICH-01, ENRICH-02, ENRICH-03)
- [ ] Mock fixtures for g:Profiler API responses (convert + GOSt)

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Interactive enrichment table sorts/filters | VIZ-01 | Requires Datagrok grid UI | Open DE result -> Run enrichment -> Verify table columns, sort by FDR, filter by source |
| Mapping statistics display | MAP-01 | UI balloon/notification | Run mapping -> Verify success rate shown |

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify or Wave 0 dependencies
- [x] Sampling continuity: no 3 consecutive tasks without automated verify
- [x] Wave 0 covers all MISSING references
- [x] No watch-mode flags
- [x] Feedback latency < 15s
- [x] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
