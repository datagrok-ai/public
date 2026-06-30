---
phase: 12
slug: spectronaut-input-coverage
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-05-15
---

# Phase 12 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.
> Derived from `12-RESEARCH.md` § Validation Architecture.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | `@datagrok-libraries/test` (`category`/`test`/`expect`), bundled, runs via `grok test` in Puppeteer/Chrome against a running Datagrok instance |
| **Config file** | none — entry `src/package-test.ts` side-effect-imports `src/tests/*` |
| **Quick run command** | `grok test --category "Spectronaut" --host localhost` |
| **Full suite command** | `grok test --host localhost` |
| **Estimated runtime** | ~30–90 s (Spectronaut category); full suite minutes |

---

## Sampling Rate

- **After every task commit:** `grok test --category "Spectronaut" --host localhost` (existing 19 + new Spectronaut tests)
- **After every plan wave:** `grok test --host localhost` (full package suite — guards downstream pipeline parity)
- **Before `/gsd-verify-work`:** Full suite green + manual 2.6 GB E2E gate (R4)
- **Max feedback latency:** ~90 s (Spectronaut category quick run)

---

## Per-Task Verification Map

| Req | Behavior | Test Type | Automated Command | File Exists |
|-----|----------|-----------|-------------------|-------------|
| R1 | Streaming import yields a DataFrame from a precursor fixture (no whole-file buffering) | unit | `grok test --test "streams precursor fixture" --host localhost` | ❌ W0 (extend `src/tests/spectronaut-parser.ts`) |
| R2 | Streaming aggregation == duckdb golden (`PG.Quantity` ≤1e-3, `EG.Qvalue` equal, same group keys) | unit (golden) | `grok test --test "streaming output equals duckdb golden" --host localhost` | ❌ W0 (new fixture + golden `.tsv`) |
| R2 | `EG.Qvalue ≤ 0.01` numeric-drop / non-numeric-pass / null-pass; `CON__`/`REV__`/null-PG drop | unit | `grok test --test "streaming filter parity" --host localhost` | ❌ W0 |
| R3 | Streaming output shape == `parseSpectronautText` output (rows, sample cols, tags, 2-cond groups) on shared fixture | unit | `grok test --test "stream path matches text path" --host localhost` | ❌ W0 |
| R3 | Tags: `proteomics.source=spectronaut`, `proteomics.preNormalized=true`, `proteomics.groups` 2 conditions | unit (tag) | covered by the equivalence test above | ❌ W0 |
| R4 | Progress + responsiveness on 2.6 GB file (monotonic bar, no Page-Unresponsive, tab switch works) | manual | release-gate, manual on `~/Downloads/2026-05-13 BP DMD WT.tsv` | n/a (manual gate per spec) |
| R5 | `tools/spectronaut-aggregate.{sql,sh}` committed; README references it; failure path hints it | unit/static | assert files exist + string test on catch-path hint; README check manual | ❌ W0 (file-presence test automatable) |
| R6 | A `grok test` for the streaming path passes against the committed small fixture | unit | `grok test --category "Spectronaut" --host localhost` | ❌ W0 |
| — | Existing 19 Spectronaut tests still green (regression) | unit | `grok test --category "Spectronaut" --host localhost` | ✅ exists (`src/tests/spectronaut-parser.ts`) |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] Extend `makeLongFormatTsv` in `src/tests/spectronaut-parser.ts` to emit precursor rows: ≥2 precursor/fragment rows per protein×sample (varying `EG.ModifiedPeptide`/`FG.Charge`/`FG.Id`), `CON__`/`REV__` decoy rows, sub-threshold + non-numeric + empty `EG.Qvalue` rows, conditions kept non-`DMD`/`WT`, carry-along cols group-constant — covers every R2 branch
- [ ] New committed synthetic precursor fixture `.tsv` under `files/demo/`
- [ ] New committed duckdb golden `.tsv` beside it (D-04) + documented one-liner regen command (mirror `node tools/generate-spectronaut-candidates-fixture.mjs` convention; here a duckdb invocation via `tools/spectronaut-aggregate.sh`)
- [ ] `tools/spectronaut-aggregate.sql` + `tools/spectronaut-aggregate.sh` committed from `/tmp` (D-05), with the `R.Condition` flip documented as reference-file-only
- [ ] New tests: streaming-vs-text equivalence, streaming-vs-golden equivalence, filter-branch parity, tag set, tools-file presence
- [ ] `files/demo/README.md` rows for the new fixture + golden, with regen one-liner and the flip caveat
- [ ] No framework install needed — existing test infra covers it

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| 2.6 GB precursor report imports with monotonic progress, no Page-Unresponsive, tab switch works mid-import, then full pipeline runs | R4 | Reference file `~/Downloads/2026-05-13 BP DMD WT.tsv` (2.6 GB) is not committable; CI uses the synthetic fixture | Import via `Proteomics \| Import \| Spectronaut Report`; observe TaskBar progress monotonic; switch browser tabs during import (must remain responsive); confirm result ≈ 8,328 proteins × 24 samples with correct tags; run Annotate→Normalize→Impute→DE→Volcano end-to-end |
| `files/demo/README.md` references the `tools/` duckdb fallback and regen one-liner | R5 | Doc-content check is qualitative | Confirm README lists the new fixture/golden rows, the regen command, and the `R.Condition` flip caveat |

---

## Validation Sign-Off

- [ ] All tasks have automated verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references (fixture, golden, tools files, new tests)
- [ ] No watch-mode flags
- [ ] Feedback latency < 90s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
