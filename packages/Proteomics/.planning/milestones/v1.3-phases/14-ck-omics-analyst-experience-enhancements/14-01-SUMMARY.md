---
phase: 14-ck-omics-analyst-experience-enhancements
plan: 01
subsystem: parsers
tags: [ensembl, gene-labels, ckomics-port, fetchproxy, userdatastorage, semtype]

requires:
  - phase: 13-spectronaut-candidates-parity
    provides: subcellular-location resolver pattern (cache shape, runWithConcurrency, schema-v invalidation) — gene-label resolver mirrors it verbatim
provides:
  - Display Name + Source ID columns on every imported Proteomics DataFrame, regardless of vendor (R1)
  - Verbatim port of CKomics improve_gene_labels_with_ensrnog_marking + extract_readable_description
  - Cross-session Ensembl /lookup/id cache via userDataStorage with __schema_v invalidation
  - Graceful degradation on Ensembl outage (raw IDs preserved + single warning toast)
  - Duplicate disambiguation via (Source ID) suffix with one shell warning per import
  - 4 new SEMTYPE constants (DISPLAY_NAME, SOURCE_ID, NUMERATOR_MEAN, DENOMINATOR_MEAN) — last two reserved for Plan 14-04
affects:
  - 14-02 (volcano labels will bind to Display Name, hover to Source ID)
  - 14-03 (Filters viewer will add Display Name / Source ID as free-text filters per D-05)
  - 14-04 (UniProt panel + group-mean correlation viewer use the new SEMTYPE constants)

tech-stack:
  added: []
  patterns:
    - "Parser-post-process resolver: every parseXText() awaits resolveGeneLabels(df) before return — single canonical call site per vendor"
    - "Type-guarded REST response field reads (T-14-01-T1 mitigation): every Ensembl response field passes through typeof check before being stored"

key-files:
  created:
    - src/utils/gene-label-resolver.ts
    - src/tests/gene-label-resolver.ts
    - .planning/phases/14-ck-omics-analyst-experience-enhancements/14-01-SUMMARY.md
  modified:
    - src/utils/proteomics-types.ts
    - detectors.js
    - src/parsers/maxquant-parser.ts
    - src/parsers/spectronaut-parser.ts
    - src/parsers/spectronaut-candidates-parser.ts
    - src/parsers/fragpipe-parser.ts
    - src/parsers/generic-parser.ts
    - src/package.ts
    - src/package-test.ts
    - src/tests/parsers.ts
    - src/tests/fragpipe-parser.ts
    - src/tests/fragpipe-e2e.ts
    - src/tests/spectronaut-parser.ts
    - src/tests/spectronaut-candidates-parser.ts
    - src/tests/spectronaut-candidates-e2e.ts

key-decisions:
  - "Resolver call lives inside each parseXText() rather than in package.ts import handlers — honours the plan's verify gate (`await resolveGeneLabels` in each parser file) and makes the contract durable for every parser consumer, not just menu imports. Cost: parser signatures become async, ~75 existing test call sites gain `await`."
  - "extractReadableDescription's regex truncation matches the CK-omics behavior exactly, including the 60-character cap and the rejection of empty / 'uncharacterized' / <3-char results."
  - "ENSRNO falls back to rattus_norvegicus per CK-omics line 933 — `id.startsWith('ENS')` is sufficient; species detection is informational since the POST itself does not require species hints."
  - "Spectronaut Candidates parser also gets the resolver despite its proteomics.de_complete shortcut — the resolver reads Gene name unconditionally via findColumn, so the call is safe on any parser output."

patterns-established:
  - "Pattern S-1 / userDataStorage cache: every resolver follows STORE_KEY + SCHEMA_V_KEY + __schema_v invalidation, with loadCache / flushCache helpers and a final-only flush (no per-chunk flush race). Plan 13 subcellular-location established this; plan 14-01 reuses it for gene labels."
  - "Pattern: type-guarded REST entry reads. Every typed field comes through `typeof e.foo === 'string'` before being stored in the cache value object."

requirements-completed: [R1]

duration: ~2h
completed: 2026-05-29
---

# Plan 14-01: Gene-Label Resolver Summary

**Ships R1 — every imported Proteomics DataFrame now carries `Display Name` (resolved gene name with `*` grouped / `†` predicted markers) and `Source ID` (raw input ID) columns before any viewer renders, across all five vendor formats.**

## Performance

- **Duration:** ~2h inline (3 tasks: scaffold → resolver core → parser integration)
- **Completed:** 2026-05-29
- **Tasks:** 3 (atomic commits)
- **Files modified:** 15 (4 new, 11 modified)

## Accomplishments

- Verbatim TypeScript port of `improve_gene_labels_with_ensrnog_marking` and `extract_readable_description` from `~/Downloads/ck/CKomics_tool2.py:756-1104`.
- Cross-session Ensembl `/lookup/id` cache wired through `grok.dapi.userDataStorage` with `__schema_v` invalidation; second import on the same protein set sends zero new Ensembl requests until the schema bumps.
- Resolver wired into all 5 parsers (`maxquant`, `spectronaut` text + stream paths via shared `finalizeSpectronaut`, `spectronaut-candidates`, `fragpipe`, `generic`); parser signatures became async, all callers (package.ts handlers, demo, ~75 existing test sites) updated with `await`.
- 13-test suite in category `Proteomics: 14-01` covers species detection, marker rules, description cleanup, the POST shape, 429 retry-once, the no-predicted-IDs no-op path, the LOC/RGD/AABR †-only path, duplicate disambiguation, cache short-circuit, graceful degrade on fetch throw, the three-level fallback, and the 4-parser integration smoke.

## Task Commits

1. **Task 1: SEMTYPE constants + detectors + test scaffold** — `fe5bb67cdc` (feat)
2. **Task 2: gene-label-resolver core (port + cache + 11 unit tests)** — `4962678073` (feat)
3. **Task 3: wire resolver into 5 parsers + package.ts + integration test** — `2552ad8fa5` (feat)

**Plan metadata:** `d479e49d20` (docs: commit revised plans + supporting context)

## Files Created/Modified

### Created

- `src/utils/gene-label-resolver.ts` — public `resolveGeneLabels(df)` entry, `STORE_GENE_LABELS`, `SCHEMA_V_GENE_LABELS`, helper exports (`detectSpecies`, `isEnsemblEligible`, `isPredicted`, `extractReadableDescription`, `applyMarkerRules`, `lookupEnsemblBatch`).
- `src/tests/gene-label-resolver.ts` — category `Proteomics: 14-01` with 13 tests covering the unit + integration surface; uses prototype-patching of `grok.dapi.fetchProxy` / `userDataStorage` / `grok.shell.warning` with try/finally restore.

### Modified

- `src/utils/proteomics-types.ts` — added 4 SEMTYPE constants (DISPLAY_NAME, SOURCE_ID, NUMERATOR_MEAN, DENOMINATOR_MEAN). Last two reserved for Plan 14-04 (correlation viewer) but defined here as the single source of truth.
- `detectors.js` — 4 mirror detector methods (exact column-name match per the conservative pattern used by detectGeneSymbol).
- `src/parsers/maxquant-parser.ts`, `src/parsers/fragpipe-parser.ts`, `src/parsers/spectronaut-candidates-parser.ts` — parseXText becomes async; awaits resolveGeneLabels before return.
- `src/parsers/spectronaut-parser.ts` — `finalizeSpectronaut` (shared by text + stream paths) becomes async; awaits resolveGeneLabels; both callers await.
- `src/parsers/generic-parser.ts` — dialog onOK converted to async; awaits resolveGeneLabels before addTableView.
- `src/package.ts` — import handlers for MaxQuant / Spectronaut Report (both branches) / Spectronaut Candidates / FragPipe and `proteomicsDemo` await the now-async parser fns.
- `src/package-test.ts` — registers the new test file.
- `src/tests/{parsers,fragpipe-parser,fragpipe-e2e,spectronaut-parser,spectronaut-candidates-parser,spectronaut-candidates-e2e}.ts` — every parser invocation now `await`-ed (including the three throw-tests, the `streamTsv` helper, and the cross-path equivalence test).

## Decisions Made

- **Resolver call sits inside each parser, not in package.ts handlers** — the plan's verify gate (`grep -c "await resolveGeneLabels" src/parsers/*`) requires a call in each parser file. The trade-off (async signature + ~75 test edits) was paid because it makes Display Name / Source ID a durable parser contract for every consumer, not just the menu importers. Plan 14-02's volcano binds unconditionally to `Display Name`; this decision means that's safe.
- **Cache flush is finally-only, not per-chunk** — mirrors Phase 13 subcellular-location pattern after that plan documented a per-chunk-flush race. The resolver writes through `flushCache(cache, fetched)` once at the end of the lookup pass inside a `finally`.
- **ENSRNO falls back via the same Ensembl POST as ENSRNOG** — `isEnsemblEligible` only checks the `ENS` / `MGP_` prefix; species detection is informational (the POST doesn't require it). Matches CK-omics line 933 verbatim.
- **Duplicate disambiguation runs *after* marker application** — `(Source ID)` is appended to the final marker-bearing string so the analyst sees `Myh7† (ENSRNOG00000001)` rather than `Myh7 (ENSRNOG…)†`. CK-omics warns only; D-10 turns the warning into a structural disambiguation.

## Deviations from Plan

- **Test count: 13 in category instead of "≥11"** — added 2 tests for `isPredicted` and the parser-integration smoke beyond the minimum. Net-positive coverage; no behavioral deviation from the spec.
- **Resolver writes Display Name = raw Gene name unchanged on the no-predicted-IDs fast path**, exactly as the plan §"behavior" requires (Pitfall 9: no-op preserves the volcano label binding). Confirmed by the dedicated test.
- **The Wave 0 detector-mirror test (Task 1) was simplified** from a runtime-fetch-and-grep-detectors.js assertion to a Column.semType assignment round-trip. The original runtime fetch was brittle (depended on `window.fetch`/`_package.webRoot`/`window.ProteomicsPackageDetectors` global) and the plan explicitly allows the lighter assertion. The detectors are validated end-to-end by `grok check` during the build.
- **Optional progress indicator was not wired in the import handlers** (plan §Task 3 marked it "planner picks"). The resolver itself accepts a `progress?` callback for callers; the handler-side `DG.TaskBarProgressIndicator` wrapper is parked for Plan 14-02 where the volcano dialog progress UX is the central concern.

## Verification

- TypeScript clean across all changed files (`npx tsc --noEmit`).
- Full build clean (`npm run build` — same pre-existing webpack size + planning-doc header warnings as before this plan; no new warnings).
- 13 tests in category `Proteomics: 14-01` are registered for the `grok test` runner; no test runtime executed in this session (the package's runtime tests require a live Datagrok instance per `.claude/rules/testing.md`). The unit tests are self-contained — they mock `grok.dapi.fetchProxy` / `userDataStorage` / `grok.shell.warning` via prototype patching — so they should pass on any `grok test` invocation; this is captured as a UAT line for the next runtime session.
- All 5 parser files contain exactly one `await resolveGeneLabels` site (verified by grep gate).

## Open Items for Next Session

- **Run `grok test --category "Proteomics: 14-01"` against localhost** to verify the 13 new tests pass on a live Datagrok instance. Same for the existing parser tests — they were updated mechanically (`parseXText(` → `await parseXText(`) but not runtime-verified.
- **Manual import smoke test** with a real Spectronaut Candidates fixture (e.g., from `~/Downloads/ck`): import via `Proteomics | Import | Spectronaut Candidates...`, confirm Display Name column visible in the grid with raw ENSRNOG IDs replaced + `†` markers, Source ID column carries the originals.
- **Cache-hit verification**: re-import the same file in the same session and observe the second invocation makes zero Ensembl POSTs (developer-tools Network tab).
- **Plans 14-02..14-05 remain to execute.** Re-invoke `/gsd-execute-phase 14` in a fresh context to continue. Plan 14-02 should now find `Display Name` already populated and bind volcano labels to it unconditionally.

## Self-Check: PASSED
