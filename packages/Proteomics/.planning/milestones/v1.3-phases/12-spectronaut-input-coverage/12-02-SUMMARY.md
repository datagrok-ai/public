---
phase: 12-spectronaut-input-coverage
plan: 02
subsystem: api
tags: [spectronaut, streaming, proteomics, parser, importer, web-streams, textdecoderstream, duckdb-parity]

# Dependency graph
requires:
  - phase: 12-spectronaut-input-coverage
    provides: "12-01 — tools/spectronaut-aggregate.{sql,sh} duckdb parity contract, precursor fixture + duckdb golden + JSON sidecar, extended makeLongFormatTsv"
provides:
  - "parseSpectronautStream(file, qValueThreshold?) — Blob.stream() + TextDecoderStream single-pass aggregator with bounded memory, duckdb filter/aggregate parity, TaskBarProgressIndicator by bytes-read, explicit setTimeout(0) macrotask yield"
  - "finalizeSpectronaut(result: PivotResult) — shared post-pivot tail so streaming + text paths are byte-identical in shape/semTypes/log2/tags/groups"
  - "tryCastDouble — duckdb nullstr + TRY_CAST AS DOUBLE mirror used for the q-value filter and quantity/qvalue aggregates"
  - "aggToPivotResult — folds the streaming aggregate Map into the existing PivotResult shape"
  - "importSpectronaut header-sniff branch (sniffIsPrecursor) routing precursor reports to the stream path; PG-level path byte-identical; D-05 duckdb-fallback failure hint"
affects: ["12-03 (equivalence/golden/parity tests pin against this streaming output)"]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "WHATWG Streams precursor importer: file.stream().pipeThrough(TextDecoderStream).getReader() + while-loop (no for-await), carry-over line buffer, trailing-buffer flush, bounded aggregate Map keyed by control-char-joined (protein × condition × replicate) — no raw-row buffering"
    - "Shared post-pivot tail (finalizeSpectronaut) so two parse entrypoints converge to byte-identical DataFrames"
    - "Explicit setTimeout(0) macrotask yield on a ~16 ms wall-clock cadence alongside a bytes-read TaskBarProgressIndicator (the stream-read await alone does not yield on OS-buffered data)"

key-files:
  created: []
  modified:
    - packages/Proteomics/src/parsers/spectronaut-parser.ts
    - packages/Proteomics/src/package.ts

key-decisions:
  - "Place sniffIsPrecursor at module scope next to requireSampleLevelData (mirrors the existing module-helper placement) rather than as a static class method, so it is reusable and matches the parser-helper style"
  - "Use grok.shell.info for the malformed-skip surface (added `import * as grok` to spectronaut-parser.ts) — minimal new import, matches the project's user-facing-notice convention"
  - "PRECURSOR_SIGNATURE_COLUMNS = [EG.ModifiedPeptide, FG.Charge, PEP.StrippedSequence] as a named const (D-01); 1 MB pre-newline sanity guard caps the sniff read"

patterns-established:
  - "Streaming aggregator parity contract: tryCastDouble mirrors duckdb nullstr+TRY_CAST; null/non-numeric q-value passes; per-group max(quantity)/min(qvalue)/first-non-null FileName+Organisms; the reference-only DMD<->WT flip is deliberately NOT ported"
  - "AGG_KEY_SEP = 0x1F (ASCII Unit Separator) for the aggregate map key — not TSV-legal, cannot collide with data"

requirements-completed: [R1, R2, R3, R4]

# Metrics
duration: 7min
completed: 2026-05-15
---

# Phase 12 Plan 02: Streaming Spectronaut Precursor Importer Summary

**`parseSpectronautStream` streams a precursor-level Spectronaut TSV via `Blob.stream()` + `TextDecoderStream`, single-pass aggregating to the same wide protein×sample DataFrame with duckdb-parity filters/aggregates, bounded memory, a bytes-read progress bar and an explicit macrotask yield — wired into `importSpectronaut` behind a header sniff while the proven PG-level `file.text()` path stays byte-identical.**

## Performance

- **Duration:** ~7 min
- **Started:** 2026-05-15T14:36:00Z
- **Completed:** 2026-05-15T14:43:00Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Extracted the post-pivot tail of `parseSpectronautText` verbatim into a module-private `finalizeSpectronaut(result: PivotResult)`; the legacy text path now ends with `return finalizeSpectronaut(result);` and is byte-identical (all 36 Spectronaut tests green before and after).
- Added `parseSpectronautStream(file, qValueThreshold=0.01)`: WHATWG `file.stream().pipeThrough(new TextDecoderStream('utf-8')).getReader()` with an explicit `while` loop (no `for await`), carry-over line buffer, header parse + required-column / quantity-column validation matching the text path, trailing-buffer flush, a bounded `Map<AggRow>` keyed by `protein\x1fcondition\x1freplicate` (never buffering raw rows), a `TaskBarProgressIndicator` advanced by bytes-read in a `try/finally`, and an explicit `await new Promise(r => setTimeout(r, 0))` macrotask yield on a ~16 ms wall-clock cadence.
- Added `tryCastDouble` (duckdb `nullstr=['','NaN','NA']` + `TRY_CAST AS DOUBLE` mirror) used for both the q-value filter and the quantity/qvalue aggregates, with a code comment documenting the intentional strictness divergence from `pivotSpectronaut`'s loose `Number()`; added `aggToPivotResult` folding the aggregate into the existing `PivotResult` shape with `sampleKey = \`${condition}_${replicate}\``; both paths share the exact tail.
- Did NOT port the reference-file-only DMD↔WT `R.Condition` flip (RESEARCH Pitfall 1 — highest-risk parity trap); verified by grep guard.
- Added module-scope `sniffIsPrecursor(file)` (reads only to the first newline with a 1 MB guard, `await reader.cancel()` in `finally`) and branched `importSpectronaut` on it: precursor reports → `parseSpectronautStream(file)`, PG-level → `parseSpectronautText(await file.text())` (byte-for-byte the original). Added a D-05 failure-path `grok.shell.warning` pointing at `tools/spectronaut-aggregate.sh` + `files/demo/README.md`, keeping the original `grok.shell.error` intact.

## Task Commits

Each task was committed atomically:

1. **Task 1: Extract finalizeSpectronaut shared tail; add tryCastDouble + parseSpectronautStream + aggToPivotResult** - `7ef16d32aa` (feat)
2. **Task 2: Add D-01 header-sniff branch and D-05 failure-path hint to importSpectronaut** - `d130f207c4` (feat)

**Plan metadata:** committed with this SUMMARY (docs)

_Plan frontmatter declared `type: execute`; Task 1's `tdd="true"` behavior tests are explicitly owned by Plan 03 (Wave 3) per the `<behavior>` block — this plan is the production code, verified by static guards + tsc + the existing 36-test regression gate, not a RED/GREEN cycle here._

## Files Created/Modified
- `packages/Proteomics/src/parsers/spectronaut-parser.ts` — added `import * as grok`; extracted `finalizeSpectronaut`; added `tryCastDouble`, `AggRow`, `AGG_KEY_SEP`, `aggToPivotResult`, exported `parseSpectronautStream`; `parseSpectronautText` now returns `finalizeSpectronaut(result)` (legacy behavior unchanged).
- `packages/Proteomics/src/package.ts` — `parseSpectronautStream` added to the parser import; module-scope `PRECURSOR_SIGNATURE_COLUMNS` + `sniffIsPrecursor`; `importSpectronaut` open callback branches on the sniff; catch block adds the D-05 duckdb-fallback `grok.shell.warning` after the unchanged `grok.shell.error`.

## Decisions Made
- `sniffIsPrecursor` placed at module scope beside `requireSampleLevelData` (consistent with the existing module-helper convention) rather than a class method.
- Added `import * as grok` to `spectronaut-parser.ts` so the streaming path can surface the malformed-line skip count via `grok.shell.info` (the plan's "surface via grok.shell.info/the progress message"); minimal, matches package convention.
- Yield cadence implemented as ~16 ms wall-clock (`performance.now()` delta) per RESEARCH Q3 resolution / D-02 (the simpler fixed-line-count fallback was the alternative; the time-based form is more responsive across variable row widths).
- `AGG_KEY_SEP` = `\x1f` (ASCII Unit Separator) for the aggregate-map key — not a TSV-legal character, so it cannot collide with any protein/condition/replicate value.

## Deviations from Plan

None - plan executed exactly as written.

The plan's `<verify>` greps and acceptance criteria all passed as specified. The Plan 01 SUMMARY had already noted that the "19 existing tests" figure was approximate (the Spectronaut category actually contains 36 tests across 3 suites); all 36 pass before and after the tail extraction, satisfying the plan's "legacy path byte-identical / existing tests stay green" requirement. No code or scope changes were needed.

## Issues Encountered
- This freshly-created git worktree had no `node_modules` (same situation Plan 01 documented). Resolved by symlinking the sibling `proteomics` worktree's `node_modules` (verified same base commit `d1d817fec8`; this plan changes zero dependencies) to run `grok build` / `grok test`. The symlink was never `git add`ed; it was dereferenced by `grok`'s build operations into an untracked, gitignored `node_modules` directory and is absent from `git status` and `git ls-files` — no repository pollution.

## Threat Surface Scan
The plan's `<threat_model>` covers the streaming/sniff surface (T-12-05..T-12-09). No new security-relevant surface beyond it: no network endpoints, no auth, no schema changes, no file writes — only client-side parsing of a user-selected file into a `DG.DataFrame`. Mitigations were implemented as specified (bounded aggregation Map + no raw-row buffering; 1 MB sniff guard + trailing-line flush; `TextDecoderStream` multibyte carry-over; malformed-row skip+count; explicit `setTimeout(0)` yield). No threat flags.

## Next Phase Readiness
- Plan 03 (Wave 3) can now pin its per-cell streaming-vs-text equivalence test (within 1e-3) and its streaming-vs-duckdb-golden test to `parseSpectronautStream` against Plan 01's `files/demo/spectronaut-hye-precursor.tsv` + golden `.tsv`/`.json` sidecar.
- `grok build` succeeds (function metadata + webpack bundle); `npx tsc --noEmit` is fully clean (0 errors); `grok test --category "Spectronaut" --host localhost` is green (17 + 18 + 1 = 36).
- The manual 2.6 GB responsiveness/memory release-gate (R1/R4) is a manual gate per the SPEC and is out of scope for automated verification here.
- No blockers.

## Self-Check: PASSED

Both modified files exist on disk (`src/parsers/spectronaut-parser.ts`, `src/package.ts`); both task commits present in git history (`7ef16d32aa`, `d130f207c4`); `npx tsc --noEmit` 0 errors; `grok build` succeeded; `grok test --category "Spectronaut"` 36/36 green.

---
*Phase: 12-spectronaut-input-coverage*
*Completed: 2026-05-15*
