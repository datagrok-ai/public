---
phase: 13-ck-omics-volcano-and-enrichment-parity
plan: 10
status: complete
gap_closure: true
requirements: [R1, R6]
files_modified:
  - src/viewers/volcano.ts
  - src/package.ts
  - src/tests/volcano.ts
tests:
  category: "Volcano"
  before: 7
  after: 9
  added:
    - "showVolcanoBusy attaches overlay, updateVolcanoBusy mutates, hideVolcanoBusy detaches"
    - "ensureLocationColumn short-circuits warm-cache path within one tick"
regression:
  category: "SubcellularLocation"
  status: pass
---

# 13-10 — In-volcano busy overlay for subcellular-location fetch

## What was built

A presentation-only fix layered on top of 13-08's mechanics.

`src/viewers/volcano.ts` gains an exported helper trio —
`showVolcanoBusy(sp, label)`, `updateVolcanoBusy(sp, label, detail?)`,
and `hideVolcanoBusy(sp)` — that attaches/updates/detaches a centered
white card on `sp.root` marked `data-volcano-busy`. The card sits at
zIndex 6 (above the counter overlay at 5, axis labels at 4) and is
non-interactive (`pointerEvents: 'none'`). `disposeVolcanoAttachments`
also sweeps `[data-volcano-busy]` so `createVolcanoPlot` re-entry tears
down any stranded overlay.

`src/package.ts` `volcanoOptions` OK handler attaches the overlay when
`colorDim === 'location'`, threads ProgressCb ticks into both `pi.update`
and `updateVolcanoBusy`, and detaches in the finally block alongside
`pi.close()`.

## Why 13-08 alone was insufficient

13-08 shipped bounded concurrency, incremental cache, accession-hash
short-circuit, and ProgressCb threading into `pi.update`. But
`DG.TaskBarProgressIndicator` lives at the bottom of the platform shell,
and the pre-OK `grok.shell.info` toast is a one-shot transient. Test 3
of the post-13-08 UAT round confirmed the user still perceives a hang
because both surfaces are off the surface they are staring at (the
volcano viewer). The fix puts the progress signal directly on that
viewer.

## Short-circuit suppression strategy

Optimistic always-attach. When `colorDim === 'location'` the OK handler
always attaches the overlay — no pre-check for whether
`ensureLocationColumn` will actually fetch. The 13-08 short-circuit
returns in microseconds when the accession-hash matches, so the overlay
shows for at most one render frame before `hideVolcanoBusy` runs in the
`finally`. The new `ensureLocationColumn short-circuits warm-cache path
within one tick` test asserts that short-circuit completes in <50ms with
exactly one `init-column:1/1` tick, guarding against a future regression
that would re-walk the accession set.

The cleaner option — recompute the accession-set hash in the handler to
decide attachment — duplicates `volcano.ts` internals (the hash function
is module-local) and creates a cross-module coupling. The optimistic
form is simpler and observably indistinguishable to the user on warm
toggles.

## What is preserved verbatim

All 13-08 mechanics in `src/analysis/subcellular-location.ts`:
`FETCH_CONCURRENCY = 6`, `runWithConcurrency`, `CACHE_FLUSH_INTERVAL_MS
= 5000`, `flushCache`, `STORE`, `SCHEMA_V`, `ProgressCb`, the
`LOCATION_HASH_TAG` short-circuit in `volcano.ts`, and the pre-OK
`shell.info` toast in `package.ts`. Verified via the grep gates in Task
3 plus a clean `SubcellularLocation` regression rerun (15/15 pass).

## Test details

- **`showVolcanoBusy attaches overlay…`** — synthetic DataFrame +
  `createVolcanoPlot`, then asserts the attach → update text → detach
  lifecycle on `sp.root`. No network, fully synchronous DOM assertions.
- **`ensureLocationColumn short-circuits warm-cache path within one tick`**
  — patches `userDataStorage` prototype so the first call resolves all
  accessions from the (faked) cache. Asserts the second call returns in
  <50ms with one `init-column:1/1` tick. `patchUserDataStorage` is
  duplicated inline in `src/tests/volcano.ts` because the source helper
  in `src/tests/subcellular-location.ts` has no `export` keyword.

## Verification

- `npm run build` succeeds.
- `grok test --category "Volcano"` — 9/9 pass (7 baseline + 2 new).
- `grok test --category "SubcellularLocation"` — 15/15 pass; no
  regression to 13-08 mechanics.
- Grep gates from Task 3 confirm `FETCH_CONCURRENCY`,
  `runWithConcurrency`, `CACHE_FLUSH_INTERVAL_MS`, `flushCache`,
  `setInterval`, `LOCATION_HASH_TAG`, `ProgressCb`, `pi.update` still
  present.

## Human-UAT follow-up

13-UAT Test 3 should be re-run on a real BP DMD/WT Spectronaut
Candidates session:
- **First-time fetch:** clicking OK on Volcano Options with Color by =
  Subcellular Location surfaces a centered white card on the volcano
  showing "Classifying subcellular locations…" + a live "N/80 (P%)"
  detail line that updates per chunk completion.
- **Warm-cache toggle:** subsequent toggles (same accession set) either
  do not show the overlay or flash it for less than one render frame
  — no strobe.
