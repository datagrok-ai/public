---
gsd_state_version: 1.0
milestone: v1.4
milestone_name: Cross-Team Review
status: ready_to_plan
stopped_at: Phase 16 complete (7/7) — ready to discuss Phase 17
last_updated: 2026-06-08T18:50:24.684Z
last_activity: 2026-06-08 -- Phase 16 execution started
progress:
  total_phases: 6
  completed_phases: 1
  total_plans: 16
  completed_plans: 16
  percent: 17
---

# Project State

## Project Reference

See: .planning/PROJECT.md (last updated 2026-06-06 — v1.4 Cross-Team Review milestone started)

**Core value:** Scientists can import proteomics data and go from raw protein quantification to differential expression with biological interpretation — all within Datagrok, with no tool-switching or file exports.
**Current focus:** Phase 17 — campaign data model

## Current Position

Phase: 17
Plan: Not started
Status: Ready to plan
Last activity: 2026-06-08

## Performance Metrics

**Velocity (v1.3 baseline carried forward for trend reference):**

- Total plans completed (lifetime): 15+
- v1.3 average duration: ~2.7 min/plan (small/inline scope)
- v1.3 total execution time: ~16 min

**By Phase (v1.4):**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| Phase 15 | TBD | — | — |
| Phase 16 | TBD | — | — |
| Phase 17 | TBD | — | — |
| Phase 18 | TBD | — | — |
| 15 | 9 | - | - |
| 16 | 7 | - | - |

**Recent Trend:**

- v1.3 closed 2026-06-05 at PASSED verdict (17/17 reqs, 14/14 UAT)
- v1.4 starts 2026-06-06; per-phase velocity unknown until first plan ships

## Accumulated Context

### Roadmap Evolution

- v1.4 milestone started 2026-06-06; four phases (15–18) appended to ROADMAP.md
- Phase numbering continues from v1.3 (last shipped: Phase 14 + hotfixes 999.1 / 999.4)
- Suggested 4-phase structure: Publishing → SPC → Campaign Data Model → Cross-Run Comparison
- All four phases derived from REQUIREMENTS.md (37 v1.4 reqs across PUB / SPC / CAMP categories)
- Phase 17 leads with Task 0 = data model BEFORE viewer code (resolves FEATURES↔ARCHITECTURE divergence per research SUMMARY)
- Phase 15 carries the round-trip gate (`assertPublishedShape`) — Pitfall 3 (tag survival across DG.Project save) is the load-bearing unknown for the whole milestone
- Phase 999.3 (publish read-only analysis) marked SUPERSEDED BY v1.4 in BACKLOG — its strategic spine became PUB-01..13
- Phase 999.2 (banner wording fix) explicitly NOT pulled in — user chose "keep v1.4 tight"

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
v1.4-specific decisions captured in research synthesis (`.planning/research/SUMMARY.md`):

- Zero new npm dependencies — all three v1.4 capabilities on v1.3 baseline + new in-package modules
- Run identity = `(instrument_id, acquisition_datetime)` (Pitfall 7 mitigation; SPC-06)
- SPC baseline LOCKED at definition with iterative outlier removal — never rolling-recomputed (Pitfall 6; SPC-05)
- `CampaignSelectionBus` with per-viewer subscription auto-eviction — NOT a copy of v1.2's `activeSubscriptions[]` module-level array (Pitfall 13; CAMP-12)
- Compound canonicalization at save time via explicit ID or Chem SMILES — fuzzy name matching rejected (Pitfall 11; CAMP-04)
- AppData CSV storage for campaigns (HitTriage precedent); Postgres schema deferred to v1.5+ if O(1000) trigger hit (STACK.md)
- Belt-and-braces tag-AND-column encoding for critical published-analysis metadata (Pitfall 3; PUB-11)

### Pending Todos

21 pending — `/gsd:check-todos` to review. Most are March 2026 dialog-polish enhancement ideas deferred at v1.3 close; not v1.4-scoped.

### Blockers/Concerns

None active for v1.4 entry.

One research-flagged plan-time unknown (Phase 15): **which `proteomics.*` tags + `Proteomics-*` semTypes + `df.name` survive `DG.Project` save/reopen** — Phase 13 round-3 (commit `e527d07ba1`) already proved layout-config tags are partially stripped. Phase 15 plan-time work includes a throw-away round-trip script to enumerate survivors and decide whether belt-and-braces column-encoding is universal or selective.

## Session Continuity

Last session: 2026-06-08T15:27:36.349Z
Stopped at: Phase 16 UI-SPEC approved
Resume file: .planning/phases/16-sample-level-spc-tracking/16-UI-SPEC.md

## Deferred Items

Items acknowledged and deferred at v1.3 milestone close on 2026-06-05 (still standing):

| Category | Item | Status | Note |
|----------|------|--------|------|
| debug | enrichment-cross-link-not-firing | diagnosed | CONFIRMED perception/layout — addressed by Phase 13's 13-09 Option-A dock fix |
| debug | volcano-options-first-fetch-no-progress | diagnosed | CONFIRMED — addressed by Phase 13's 13-08 + 13-10 in-volcano busy overlay |
| debug | limma-de-slow | diagnosed | Performance debt from 2026-03 — carry forward to a future perf milestone |
| debug | showheatmap-hangs | diagnosed | Performance debt from 2026-03 — carry forward to a future perf milestone |
| uat | Phase 12 12-UAT.md | diagnosed | 0 open scenarios — false positive in audit-open scanner |
| uat | Phase 13 13-UAT.md | diagnosed | 0 open scenarios — same false positive shape; round-3 live UAT passed 2026-06-04 |
| verification | Phase 13 13-VERIFICATION.md | human_needed | Stale frontmatter — cosmetic, no action needed |
| todos | 20 pending enhancement todos | open | Mostly 2026-03 enhancement ideas — backlog for a future milestone, not v1.4-scoped |
| backlog | Phase 999.2 banner wording | open | Explicitly kept out of v1.4 (user direction: "keep v1.4 tight") |
| backlog | Parser expansion (FragPipe MaxLFQ, DIA-NN) | open | Not v1.4-scoped; future milestone |
| backlog | CK-omics cluster C (consolidated HTML export, custom contrast labels, organism selection) | open | Not v1.4-scoped |

None of these items block v1.4 progress.

## Operator Next Steps

- Run `/gsd:plan-phase 15` to begin the Read-Only Publishing Foundation phase
- Phase 15 plan should include the `assertPublishedShape` round-trip research task as a leading deliverable
