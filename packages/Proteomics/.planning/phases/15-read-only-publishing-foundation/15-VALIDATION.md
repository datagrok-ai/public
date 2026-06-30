---
phase: 15
slug: read-only-publishing-foundation
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-06-07
---

# Phase 15 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
> Detailed test design lives in `15-RESEARCH.md` §"Validation Architecture".

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | `@datagrok-libraries/test ^1.1.0` (already installed) |
| **Config file** | `src/package-test.ts` (test entry — exports tests array for `grok test`) |
| **Quick run command** | `grok test --category Publishing` |
| **Full suite command** | `grok test` |
| **Estimated runtime** | ~30–60 seconds for Publishing category once Wave 0 spike is locked |

Tests are registered in `src/tests/` and exported via `src/package-test.ts`. The Proteomics package already runs `grok test` against a live Datagrok server (no mocks for round-trip semantics — Pitfall 3 demands real save/reopen).

---

## Sampling Rate

- **After every task commit:** Run `grok test --category Publishing --test "<specific test name>"` for the tests directly exercised by the task.
- **After every plan wave:** Run `grok test --category Publishing` (full Publishing-category suite).
- **Before `/gsd:verify-work`:** Full suite green: `grok test --category Publishing` AND `grok test --category Pipeline` (regression guard — no upstream pipeline tag/SEMTYPE invariants regressed by the new publishing primitives).
- **Max feedback latency:** ~60 seconds per task commit (single-test invocation); ~3 minutes per wave (category suite).

---

## Per-Task Verification Map

> This map is populated by the planner once PLAN.md files exist. Each `<task>` should land one row; `automated` block on each task points at the row.

| Task ID | Plan | Wave | Requirement | Threat Ref | Secure Behavior | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|------------|-----------------|-----------|-------------------|-------------|--------|
| 15-00-01 | 00 | 0 | PUB-01..PUB-13 (spike) | — | publish-spike.ts prints actual save/reopen survival for tags, semTypes, df.name, Project.options, viewer position against the live server | spike | `grok test --test "Publishing: spike"` | ❌ W0 | ⬜ pending |
| 15-XX-YY | … | … | PUB-XX | … | … | unit/integration | `grok test --test "Publishing: <name>"` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

**Load-bearing tests (per RESEARCH.md §Validation Architecture):**

| Test | Covers | Source criterion | Pitfall |
|------|--------|------------------|---------|
| `Publishing: assertPublishedShape round-trip` | Success criterion 3 (tags + columns + df.name + viewer dock position survive save/reopen) | SC-3 | Pitfall 3 |
| `Publishing: source mutation after publish leaves clone unchanged` | Success criterion 4 (frozen clone) | SC-4 | Pitfall 1 |
| `Publishing: verify-and-rollback rejects Edit-inherited grant` | Success criterion 5 (view-only ACL enforced) | SC-5 | Pitfall 2 |
| `Publishing: no-enrichment source publishes cleanly` | D-05 enrichment carry is opportunistic | — | — |
| `Publishing: with-enrichment source carries both DFs + cross-DF wiring` | D-05 enrichment-DF carry + protein-highlight subscription re-establishes on reopen | — | — |
| `Publishing: republish creates v2 + superseded_by pointer + bidirectional supersedes` | Success criterion 5 (republish chain) | SC-5 | Pitfall 4 |
| `Publishing: Spectronaut Candidates fixture publishes via de_complete shortcut` | Coverage for both source-path shapes | — | — |
| `Publishing: published_target slug sanitization handles freeform input` | D-01 slug rules | — | — |

---

## Wave 0 Requirements

Wave 0 is a one-shot spike against the live Datagrok server to resolve assumptions A1, A2, A3, A4, A6, A8 from RESEARCH.md before the trim contract gets locked. This is the load-bearing technical unknown for the phase.

- [ ] `src/tests/publish-spike.ts` — enumerates which `proteomics.*` tags, `Proteomics-*` semTypes, `df.name`, `Project.options`, and viewer dock positions survive `DG.Project.save` → `find(id).open()` in a fresh session. Output drives the trim contract and the belt-and-braces metadata column column-list.
- [ ] `src/tests/fixtures/` — synthetic demo fixture (small, deterministic — reuse v1.0 test fixture) + a Spectronaut Candidates micro-fixture (covers the `proteomics.de_complete` shortcut path where DE is pre-computed by the parser).
- [ ] Register `Publishing` test category in `src/package-test.ts` so `grok test --category Publishing` resolves.

*Framework install: not required — `@datagrok-libraries/test ^1.1.0` and the `grok test` pipeline are already in place; the Proteomics package has existing tests under `src/tests/`.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Dialog UX flow (target input + reviewer group ChoiceInput + note + republish banner pre-fill) | PUB-04, PUB-08, PUB-10 | UI flow inside `ui.dialog(...)` is not exercisable inside `grok test`; needs interactive validation | Open a DE-complete table, `Proteomics → Share → Share Analysis for Review...`, exercise target free-text + group picker + note + (after first publish) re-open dialog and assert prefill + republish banner appears |
| Audit context panel auto-docks on Project open and is biologist-readable (Pitfall 14 jargon audit) | PUB-07 | Panel rendering + biologist-comprehension is a human-judgment surface | Open published Project as a reviewer impersonation (`grok s` user-switch), confirm panel auto-docks on a PROTEIN_ID column row-current, verify wording avoids `DataFrame`/`tag`/`semType`/`ACL`/`viewer factory` |
| Mailto launch (PUB-13) | PUB-13 | Browser `mailto:` URL opens the user's mail client — outside test runner | Click the mailto button in the audit panel; assert OS launches mail client with subject + body pre-filled |
| Pre-demo dress rehearsal with biologist-class user | Pitfall 14 recovery | Dress-rehearsal is by definition a human review pass | Schedule with one biologist before any Cytokinetics demo touch surfaces; capture jargon audit findings |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies (populated by planner)
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify (planner enforces task ordering)
- [ ] Wave 0 covers all MISSING references (spike script lands first)
- [ ] No watch-mode flags
- [ ] Feedback latency < 60 seconds (single-test grok test invocation)
- [ ] `nyquist_compliant: true` set in frontmatter (after planner populates per-task map and Wave 0 spike resolves)

**Approval:** pending
