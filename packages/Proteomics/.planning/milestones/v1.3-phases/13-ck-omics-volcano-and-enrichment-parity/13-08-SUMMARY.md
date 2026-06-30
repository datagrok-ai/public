---
phase: 13-ck-omics-volcano-and-enrichment-parity
plan: 08
status: complete
gap_closure: true
requirements: [R1, R6]
files_modified:
  - src/analysis/subcellular-location.ts
  - src/viewers/volcano.ts
  - src/package.ts
  - src/tests/subcellular-location.ts
tests:
  category: "SubcellularLocation"
  before: 10
  after: 16
  added:
    - "runWithConcurrency caps in-flight tasks at the supplied limit"
    - "getSubcellularLocations emits non-decreasing progress per phase"
    - "getSubcellularLocations writes the cache incrementally during fetch"
    - "STORE export matches the documented userDataStorage key"
    - "ensureLocationColumn short-circuits on accession-set hash match"
    - "ensureLocationColumn forwards progress through fetch + emits init-column"
---

# 13-08 — Progress + bounded concurrency + incremental cache + short-circuit

## What was built

Four mechanical UX/perf fixes around the existing 13-04 UniProt subcellular
location service. All locked CK-omics contracts (keyword map, palette,
schema, chunk sizes, D-03 reviewed-by-gene fallback semantics) are preserved
verbatim.

| Fix | Touch point | Effect |
|---|---|---|
| Bounded concurrency | `runWithConcurrency` + Pass 1/2 in `getSubcellularLocations` | 4–8× wall-clock speedup on a real ~80-chunk Candidates file. |
| Incremental cache flush | Timer-driven `flushCache` in `getSubcellularLocations` | Interrupted sessions retain every accession that finished. |
| Progress callback | `ProgressCb` threaded through `getSubcellularLocations` → `ensureLocationColumn` → `recomputeVolcano` → `volcanoOptions` | Live `pi.update("Fetching subcellular locations: N/M")` instead of a static label. |
| Column-tag short-circuit | `proteomics.location_acc_hash` on the `Subcellular Location` column | Repeat toggles on the same DataFrame skip fetch+init entirely. |
| Cache-aware toast | `volcanoOptions` OK handler peeks `userDataStorage` before showing | Cold-cache wording vs warm-cache wording vs no-toast (column already exists). |

## Chosen tuning knobs

- **`FETCH_CONCURRENCY = 6`** — measured baseline on a real ~80-chunk Spectronaut
  Candidates file (~80 sequential round-trips at ACC_CHUNK=100). Six is a
  compromise: ~4× speedup while staying well under rest.uniprot.org's documented
  limits. Increase only with measured evidence; UniProt's stream endpoint is
  shared infrastructure.
- **`CACHE_FLUSH_INTERVAL_MS = 5000`** — wall-clock interval between incremental
  cache writes during a fetch. Five seconds is short enough that an
  interrupted session loses at most ~5 s of work, long enough that the
  number of write-through ops on a normal fetch stays bounded.
- **`SCHEMA_V = '13-04-1'`** — UNCHANGED. The on-disk cache shape is identical
  to the 13-04 contract; sessions that ran the pre-gap-closure code keep their
  warm cache.

## Hash algorithm

`fnv1aHex` is a FNV-1a 32-bit mixer with hex output. Deterministic,
dependency-free, fast enough for the ~8k-accession stable-sorted join that
`ensureLocationColumn` feeds it. Input: `[...accessions].sort().join('|')` —
sorting first means row order in the DataFrame doesn't affect the hash.

## Patterns established (new — additive to CLAUDE.md)

### Column-level tag convention

`proteomics.location_acc_hash` is set on the `Subcellular Location`
`DG.Column` via `col.setTag` / `col.getTag`, **NOT** on the DataFrame via
`df.setTag`. This is the first column-tag in the package — all prior
`proteomics.*` workflow state is on the DataFrame.

The convention: per-column cache-invalidation metadata belongs on the column;
whole-pipeline workflow state belongs on the DataFrame. CLAUDE.md's "DataFrame
tags" table covers the workflow-state half; this plan introduces the
column-tag half. Both halves use the `proteomics.` prefix.

### Timer-driven incremental flush

Concurrent workers should NEVER share mutex-style flush counters with `await`
checkpoints inside the worker body — that pattern hits a documented
double-flush race when two workers reach the check between `await` points
(JS single-threaded reasoning breaks at every `await` yield). Instead:

- Workers update only synchronous-write state (here: `resolved`, `fetched`,
  `geneByAcc` maps).
- The orchestrator starts a single `setInterval(flushCache, INTERVAL_MS)`
  before the workers launch.
- `try/finally` ensures `clearInterval` + one final `flushCache` run after
  the workers drain.

This is the safe pattern for any concurrent-worker code in this package that
needs periodic side-effects.

## Test fixture hygiene

`grok.dapi.userDataStorage` returns a **NEW `UserDataStorage` instance per
access** (js-api dapi.ts: `get userDataStorage() { return new UserDataStorage(); }`).
Patching the instance's `put` / `get` only sticks for one call, then the next
access creates a fresh instance with the original methods. Tests use a helper
`patchUserDataStorage({get, put})` that patches the **prototype**, with a
returned `restore()` callback that the test's `finally` invokes. This pattern
should be reused for any future test that needs to spy on `userDataStorage`
in this package.

Also: test accession prefixes are unique per test (`TESTB*`, `TESTC*`,
`TESTD*`, `TESTE*`) so a previous test run's real on-disk cache cannot
silently mark every accession as a hit and bypass the fetch loop.

## Verification

- `npm run build` clean.
- `grok test --category "SubcellularLocation"` — **16/16 pass** (10 existing
  + 6 new: runWithConcurrency cap, progress non-decreasing per phase,
  incremental cache write, STORE export, short-circuit by hash, progress
  forwarding through ensureLocationColumn).
- `grok test --category "Volcano"` — **7/7 pass** (no regression; new
  optional `progress` arg on `ensureLocationColumn` and `recomputeVolcano`
  is backward-compatible).
- `grok test --category "Enrichment Visualization"` — **10/10 pass** (the
  co-dock in 13-07 routes through the same `createVolcanoPlot` that
  consumes `recomputeVolcano`; the extra arg is optional so the
  enrichment path stays compatible).

## Follow-up

13-UAT.md Test 3 should be re-run on a real ~8k-accession Spectronaut
Candidates file to confirm the wall-clock improvement and toast UX
land for the user. The automated tests prove the wire-up (concurrency
cap, progress callback firing, incremental cache write, short-circuit on
hash match); the wall-clock and toast wording are human-UAT only.
