---
phase: 12
slug: spectronaut-input-coverage
status: secured
threats_open: 0
threats_total: 14
threats_closed: 14
asvs_level: 1
created: 2026-05-15
---

# SECURITY — Phase 12: Spectronaut Input Coverage

**Audited:** 2026-05-15
**Auditor:** gsd-security-auditor (verification, not discovery)
**Scope:** Threat register T-12-01..T-12-14 (PLAN 12-01 / 12-02 / 12-03 `<threat_model>` blocks)
**ASVS Level:** default (no auth, no network endpoints, no PII; client-side parsing of a user-selected file into a `DG.DataFrame`)
**Result:** SECURED — 14/14 threats resolved (10 mitigations verified present in code, 4 accepted risks documented and confirmed against code)

This audit verifies each declared mitigation is actually present in the implemented code (not in
documentation or intent). Implementation files were treated as read-only.

---

## Threat Verification

| Threat ID | Category | Disposition | Status | Evidence (file:line) |
|-----------|----------|-------------|--------|----------------------|
| T-12-01 | Tampering | mitigate | CLOSED | `tools/spectronaut-aggregate.sql:30-42` loud `REFERENCE-FILE-ONLY` warning block around the `CASE "R.Condition" WHEN 'DMD'…WHEN 'WT'` flip (L43-47); `tools/spectronaut-aggregate.sql:20-28` `DIVERGENCE FROM /tmp` comment, no `PG.Genes`/`PG.ProteinAccessions` token anywhere in the SQL; `src/parsers/spectronaut-parser.ts` confirmed to contain NO `WHEN 'DMD'`/`WHEN 'WT'`/`'DMD'`/`'WT'` (grep clean); fixture uses `CondA`/`CondB` (`tools/generate-spectronaut-precursor-fixture.mjs:36`) so flip is a structural no-op; golden produced (`files/demo/spectronaut-hye-precursor-golden.tsv` 9.4 KB, committed). |
| T-12-02 | Tampering | mitigate | CLOSED | `files/demo/README.md:148-166` "Flip caveat & regeneration order" documents the mandatory 3-step chain in order; `tools/derive-precursor-golden-sidecar.mjs:7-12` does pure transcription, no aggregation logic (reads golden `.tsv`, parses columns, emits JSON — no filter/group/max/min); Plan-03 golden test `streaming output equals duckdb golden` (`src/tests/spectronaut-parser.ts:472-522`) loads the committed sidecar and asserts exact key-set + per-key quantity ≤1e-3, failing loudly on drift. |
| T-12-03 | Information disclosure | accept | CLOSED (accepted) | Fixture fully synthetic — `tools/generate-spectronaut-precursor-fixture.mjs` builds rows from in-script constants only (no client/external data source); `files/demo/README.md:175` License: "Synthetic data, no restrictions." See Accepted Risks below. |
| T-12-04 | Denial of service | accept | CLOSED (accepted) | `files/demo/spectronaut-hye-precursor.tsv` = 41,597 B (~42 KB), golden = 9,374 B, sidecar = 15,867 B — all tens of KB, comparable to the existing ~20 KB candidates fixture. Bounded by `N_REGULAR = 36` + 5 edge proteins × 2 conditions × 3 replicates (`generate-spectronaut-precursor-fixture.mjs:83,37`). See Accepted Risks below. |
| T-12-05 | Denial of service | mitigate | CLOSED | `src/parsers/spectronaut-parser.ts:282-432` `parseSpectronautStream`: only retained state is the bounded carry-over `buffer` string (L292, sliced after each `\n` L383-384) and the `agg` Map keyed by `protein\x1fcondition\x1freplicate` (L291,324) — bounded by #proteins×#samples, not input rows. Grep confirms NO raw-line array push, NO `lines.push`, NO whole-text `.split('\n')`, NO `file.text()`/`readAsText` in the streaming path. Doc block L262-265 states the invariant. |
| T-12-06 | Denial of service | mitigate | PARTIAL → CLOSED with documented residual | `src/package.ts:93-94` `maxBytes = 1 MB` pre-newline sanity guard caps the sniff read; main loop `src/parsers/spectronaut-parser.ts:375-409` streams chunk-by-chunk with no per-line size cap; trailing partial line flushed after `done` (L411-420). Residual: REVIEW WR-05 — a precursor header line >1 MB truncates the sniffed header and silently misroutes the file to the non-streaming `file.text()` path. See Residual Risk note below; threat's *declared* mitigation (1 MB guard + chunked loop + trailing flush) is present and verified, so CLOSED, with WR-05 logged as a tracked residual for the routing-robustness improvement. |
| T-12-07 | Tampering | mitigate | CLOSED | `src/parsers/spectronaut-parser.ts:287-289` `file.stream().pipeThrough(new TextDecoderStream('utf-8')).getReader()` — stateful WHATWG decoder carries multibyte across chunk boundaries; malformed rows skipped+counted via `handleFields` returning false → `skipped++` (L307-308,320,392-393), surfaced via `grok.shell.info` (L425-426) not fatal — duckdb `ignore_errors=true` parity. |
| T-12-08 | Tampering/Spoofing | accept | CLOSED (accepted) | `src/package.ts:87-106` `sniffIsPrecursor` is header-only and pure (reads to first `\n`, `await reader.cancel()` in `finally` L103-105, no side effects, returns boolean). Misroute degrades safely: a spurious signature column → streaming path (still correct via shared `finalizeSpectronaut`); a precursor file lacking all three → `file.text()` fast-fail + the D-05 fallback hint. No secrets/auth. See Accepted Risks below. |
| T-12-09 | Denial of service | mitigate | CLOSED | `src/parsers/spectronaut-parser.ts:399-408` explicit `await new Promise<void>((r) => setTimeout(r, 0))` macrotask yield on a `performance.now()`-tracked ~16 ms cadence (L302,400,407), paired with `pi.update(...)` bytes-read progress (L401-405); `DG.TaskBarProgressIndicator` created L285, closed in `finally` L429-431. |
| T-12-10 | Tampering | mitigate | PARTIAL → CLOSED with documented residual | `src/tests/spectronaut-parser.ts:472-522` golden test asserts streamed (protein×sample) key-set EQUALS the JSON sidecar exactly + per-key quantity ≤1e-3; `streaming filter parity` (L430-458) and per-cell `stream path matches text path` (L379-428) lock edge/value behavior; parser grep confirms no ported flip. Residual: REVIEW WR-01/WR-02/WR-03 are real text-vs-stream-vs-duckdb divergences on *non-fixture* input the synthetic fixture cannot exercise. The declared mitigation (golden key-set + per-cell + filter-parity tests) is present and would catch a ported flip / aggregator key drift, so CLOSED; the WR-01/02/03 correctness gaps are logged as tracked residuals (they degrade real-input correctness but do not defeat the anti-flip regression net). |
| T-12-11 | Tampering | mitigate | CLOSED | `src/tests/spectronaut-parser.ts:477-478` golden test loads the committed sidecar via `_package.files` (no in-test re-aggregation, no conditional fallback — confirmed `! grep "if committed-file reads"`); regen-without-rederive fails the exact key-set/per-cell assertions L487-522 on the next `grok test` run, making drift loud. |
| T-12-12 | Repudiation | mitigate | CLOSED | `tools/generate-spectronaut-precursor-fixture.mjs:92-93,97-122` seeds `CON__`/`REV__`, a mixed >0.01 sub-threshold protein, a `Profiled` non-numeric, and an empty-string q-value protein; `src/tests/spectronaut-parser.ts:430-458` `streaming filter parity` asserts each branch (CON__/REV__/HIGHQ excluded; PROF/EMPTYQ included; rowCount==4); per-cell equivalence test (L379-428) rejects a no-op aggregator; `sniffIsPrecursor routes by header` (L355-377) pins routing. |
| T-12-13 | Denial of service | accept | CLOSED (accepted) | Same bound as T-12-04: fixture ~42 KB, sidecar ~16 KB — no CI stall surface; tests read via `_package.files.readAsText` (`src/tests/spectronaut-parser.ts:320-322`) of a tens-of-KB file. See Accepted Risks below. |
| T-12-14 | Repudiation | mitigate | CLOSED | All three plans' `<automated>` gates parse the runner's machine-readable failed-count via `node -e '... s.match(/[Ff]ailed[:\s]+(\d+)/) ... process.exit(m[1]==="0"?0:1)'` (PLAN 12-01 Task 2, 12-03 Task 1/2 `<verify>` blocks), not a substring grep of "pass"/"fail". SUMMARYs 01/02/03 document the documented clean-pass interpretation (`Tests passed.` + exit 0 + no `[1-9]\d* failed`) consistently — deterministic, not brittle substring. |

---

## Accepted Risks Log

The following threats were dispositioned `accept` at plan time. Each rationale was re-checked
against the implemented code and holds.

### T-12-03 — Information disclosure: client data in committed fixture
**Accepted.** `tools/generate-spectronaut-precursor-fixture.mjs` constructs every row from
hard-coded constants (`PRECURSORS`, `ORGANISMS`, synthetic `P00000…` ids, `CON__P99999`,
`REV__Q88888`). No file read, no network, no external/client data source. `files/demo/README.md`
License section explicitly marks it "Synthetic data, no restrictions." No information-disclosure
surface. Rationale holds.

### T-12-04 — DoS: fixture bloats the repo
**Accepted.** Fixture 41,597 B, golden 9,374 B, sidecar 15,867 B — all tens of KB, smaller than
or comparable to the existing 20 KB `spectronaut-hye-candidates.tsv`. Generation is bounded
(`N_REGULAR = 36` + 5 edge proteins × 2 conditions × 3 replicates × 2 precursors). Not a
repo-size DoS vector. Rationale holds.

### T-12-08 — Tampering/Spoofing: header column-name routing spoof
**Accepted.** `sniffIsPrecursor` (`src/package.ts:87-106`) is pure and header-only: reads only
to the first newline, `await reader.cancel()` in `finally`, returns a boolean from
`PRECURSOR_SIGNATURE_COLUMNS.some(c => header.includes(c))`. Both misroute directions degrade
safely to existing behavior (streaming path is a correct superset via the shared
`finalizeSpectronaut` tail; a non-precursor-classified precursor file fast-fails with the
existing error plus the D-05 fallback hint). No auth, no secrets, no privilege or
data-integrity escalation. Exposing the function widens no attack surface. Rationale holds.
*(Note: REVIEW IN-02 observes the `includes()` substring match is looser than a tab-split
exact-token check — an informational hardening suggestion, not a change to this accepted
disposition; the worst case remains a safe degrade.)*

### T-12-13 — DoS: test reads a huge committed file and stalls CI
**Accepted.** Same bound as T-12-04. Tests read the ~42 KB fixture / ~16 KB sidecar via
`_package.files.readAsText`. No CI-stall surface. Rationale holds.

---

## Residual Risks (tracked, non-blocking)

The independent code review (`12-REVIEW.md`, 0 critical / 5 warning) surfaced correctness
divergences that the synthetic fixture deliberately avoids. These do NOT defeat any declared
threat mitigation (each threat's declared mitigation was verified present), but they degrade
correctness/robustness on real Spectronaut input and are recorded here so they are not lost:

- **WR-05 (cross-ref T-12-06):** `sniffIsPrecursor`'s 1 MB pre-newline guard, if tripped by a
  >1 MB precursor header line, silently truncates the sniffed header and routes a multi-GB
  file to the non-streaming `file.text()` path — the exact OOM this phase exists to prevent.
  T-12-06's *declared* mitigation (1 MB guard + chunked main loop + trailing-line flush) is
  present and verified; the silent-misroute failure mode is a routing-robustness residual
  (suggested fix: treat "1 MB, no newline" as an explicit signal — route to streaming or
  throw, never treat the truncated buffer as the full header).
- **WR-01 (cross-ref T-12-10):** truncated-row handling uses `f.length < expectedFields`
  where `expectedFields` is the max *used* column index +1, not the full header width;
  diverges from duckdb `ignore_errors` (which rejects on header-width mismatch) for rows
  truncated after the last used column. Fixture emits only full-width rows so the golden
  test cannot expose it.
- **WR-02 (cross-ref T-12-10):** a group whose q-value passes but whose every quantity cast
  is null is dropped entirely by the streaming path (`aggToPivotResult` L236) but emitted
  all-null by the text/duckdb path — a structural divergence the constant-quantity fixture
  masks.
- **WR-03 (cross-ref T-12-10/T-12-12):** text path takes first-encountered quantity, stream
  path + duckdb take `max`; they coincide only because the fixture holds quantity constant
  per group. The `stream path matches text path` test proves equivalence only for
  constant-quantity input.

T-12-10's declared mitigation (golden key-set + per-cell + filter-parity regression net that
catches a ported flip or aggregator key drift) is present and verified — hence CLOSED — but
the WR-01/02/03 real-input correctness gaps are genuine and should be addressed or explicitly
documented in the parser's parity doc block (`src/parsers/spectronaut-parser.ts:267-281`).
Routing to a follow-up correctness work item is recommended; none is a security BLOCKER for
this phase (no crash, no data exfiltration, no privilege escalation — client-side parse only).

---

## Unregistered Flags

None. The Plan 02 and Plan 03 SUMMARY "Threat Surface Scan" sections both explicitly state
"No threat flags" and confirm no new security-relevant surface (no network endpoints, no auth,
no schema changes, no runtime file writes — only client-side parsing of a user-selected file
into a `DG.DataFrame`). Verified against the code: the only new entry points are
`parseSpectronautStream` and `sniffIsPrecursor`, both covered by the registered threats
T-12-05..T-12-09; no unmapped attack surface appeared during implementation.

---

## Verdict

All 14 threats resolved: 10 `mitigate` mitigations verified present in the implemented code,
4 `accept` risks documented with rationale re-checked against the code. No declared mitigation
is absent or contradicted. Phase 12 is **SECURED**. The WR-01/02/03/05 review findings are
real-input correctness/robustness residuals (logged above), not gaps in any declared threat
mitigation, and are non-blocking for the security disposition of this phase.
