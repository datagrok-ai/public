---
status: diagnosed
trigger: "Switching the Volcano Options Color by from significance to Subcellular Location for the first time in a session feels like a hang — no progress UI during the chunked UniProt fetch"
created: 2026-05-22T00:00:00Z
updated: 2026-05-22T00:00:00Z
---

## Current Focus

hypothesis: CONFIRMED — three compounding root causes:
  (1) UniProt chunks are awaited STRICTLY SEQUENTIALLY (~80 chunks for 8k accessions × per-chunk RTT)
  (2) The TaskBarProgressIndicator that IS created in the OK handler is never `.update(...)`d — it shows an indeterminate "Updating volcano..." label only, with no percent, no chunk count, no phase
  (3) ensureLocationColumn re-runs the full fetch+parse on every toggle (no short-circuit when the LOCATION_COL is already populated for the same accession set); only the userDataStorage cache (sub-call) saves the round-trip
test: Inspected the three named code paths and confirmed all three conditions by reading.
expecting: (n/a — diagnosed)
next_action: Return ROOT CAUSE FOUND

## Symptoms

expected: Toggling Volcano Options metric or Color by returns control with timely user feedback. First-time Color by location should show some progress signal during chunked UniProt fetch.
actual: "changing options kicks off some calculations that never return" → later amended to "eventually returned". Spec functionally met but no progress UI during multi-second to multi-minute operation.
errors: None.
reproduction: Phase 13 UAT test 3. With Spectronaut Candidates file (~8k unique accessions) and volcano open, run Proteomics | Visualize | Volcano Options… first time in session (cache empty), switch Color by significance→Subcellular Location, click OK.
started: 2026-05-21 during 13-UAT testing; 13-04 implementation completed 2026-05-17.

## Eliminated

(none yet)

## Evidence

- timestamp: 2026-05-22T00:00:00Z
  checked: src/analysis/subcellular-location.ts:227-247 (Pass 1 accession-chunk loop) and 260-286 (Pass 2 gene-fallback loop)
  found: Both loops are `for (const group of chunk(...))` with `await grok.dapi.fetchProxy(url)` INSIDE the body. No Promise.all anywhere in the file. Each batch waits for the previous batch's network round-trip to fully complete before the next URL is even constructed.
  implication: For 8k unique accessions at ACC_CHUNK=100 → 80 sequential UniProt round-trips. Even at a fast 1–2 s/chunk, total wall time is 80–160 s of opaque waiting. UniProt-side variability can push individual chunks to 5–10 s, easily yielding multi-minute first-time fetches. This alone produces "calculations that never return" behaviour from the user's POV.

- timestamp: 2026-05-22T00:00:00Z
  checked: src/package.ts:317-324 (volcanoOptions OK handler)
  found: `const pi = DG.TaskBarProgressIndicator.create('Updating volcano...');` is created, and `pi.close()` runs in `finally`. There is NO `pi.update(...)` call anywhere in package.ts, volcano.ts, or subcellular-location.ts (grep across all three confirmed). The fetch loop in subcellular-location.ts has no progress callback, no `onProgress` parameter, no shell.info between chunks. Only `console.warn` on per-batch failures (silent to the user).
  implication: The user sees a static, indeterminate task-bar entry labelled "Updating volcano..." for the entire 80–160+ second wall-clock duration with zero indication of (a) which phase is running (fetch vs DataFrame init), (b) progress through the 80 chunks, (c) whether anything is happening at all. UI is indistinguishable from a hang. This is the dominant user-facing symptom.

- timestamp: 2026-05-22T00:00:00Z
  checked: src/viewers/volcano.ts:87-116 (ensureLocationColumn) — searched for an early-exit on `df.col(LOCATION_COL)` having values
  found: The function ALWAYS walks the id column, ALWAYS calls `getSubcellularLocations([...accessions])`, ALWAYS bulk-inits the column. The line `const col = df.col(LOCATION_COL) ?? df.columns.addNewString(LOCATION_COL);` reuses the column object (Pitfall 5 — preserves scatter binding) but does not skip the fetch/parse work when the column already holds valid data.
  implication: After a Q↔P metric toggle that keeps `colorDim='location'`, the entire location pipeline re-runs. The `userDataStorage` cache (subcellular-location.ts:189-198) saves the network round-trip on cache hit, but every call still pays: id-column scan, accession dedup, cache load (one userDataStorage.get round-trip), Map→object conversion, full re-init of N row values. The user-reported "second toggle eventually returned" suggests the cache is doing its job on the warm path, but the work is not negligible on large tables and is wholly invisible.

- timestamp: 2026-05-22T00:00:00Z
  checked: src/analysis/subcellular-location.ts:296-303 (cache write-through)
  found: `grok.dapi.userDataStorage.put(STORE, {...cache, ...fetched, [SCHEMA_KEY]: SCHEMA_V})` runs ONCE at the very end of getSubcellularLocations, AFTER both passes complete. If the user navigates away or refreshes mid-fetch (very plausible during a multi-minute opaque wait), every fetched accession is lost and the next session restarts from zero.
  implication: Atomic write-through compounds the perceived-hang problem: an impatient user who reloads gets no partial credit, recreating the multi-minute wait next session. Not a root cause of THIS bug, but a load-bearing constraint for any fix.

- timestamp: 2026-05-22T00:00:00Z
  checked: src/viewers/volcano.ts:179-197 (recomputeVolcano) and src/package.ts:311-326 (Volcano Options dialog lifecycle)
  found: `ui.dialog('Volcano Options').onOK(async () => { ... }).show()` — the dialog closes immediately on OK click (datagrok dialog default behaviour); the awaited `recomputeVolcano` runs AFTER the dialog has dismissed. So the dialog itself doesn't block. But there is no replacement modal/spinner — the user sees the volcano in its old state (significance colouring) until the fetch completes, then it snaps to the new state.
  implication: Combined with (1) and (2), the user has zero visual feedback between OK-click and final render. No spinner, no "fetching subcellular locations…" toast, no progress bar tick. Pure dead time.

## Resolution

root_cause: |
  The Phase 13-04 implementation lacks any user-facing progress reporting during what is intrinsically a slow operation. Three compounding causes:
  (R1) [PRIMARY] No progress UI signal — `volcanoOptions` creates a TaskBarProgressIndicator but never calls `pi.update(percent, description)`; `getSubcellularLocations` exposes no progress callback. The static indeterminate label "Updating volcano..." is the only signal across 80+ chunks of UniProt traffic.
  (R2) [PRIMARY] Sequential chunk awaiting in `getSubcellularLocations` (subcellular-location.ts:227 and :260) — both fetch loops are `for…of` with `await` in the body. For 8k accessions at ACC_CHUNK=100 = 80 strictly sequential round-trips, plus a second pass over unresolved-with-gene that adds another GENE_CHUNK=20 sized batch series. Each round-trip is to rest.uniprot.org (1–10 s variable). Wall-clock easily reaches 80–300 s on first run.
  (R3) [SECONDARY] No short-circuit in `ensureLocationColumn` — every Volcano Options OK click that selects `colorDim='location'` re-walks the id column, re-calls getSubcellularLocations (which itself does a userDataStorage round-trip + Map allocation), and re-inits the column, regardless of whether the column already holds the same values. The cache prevents the network re-fetch on the warm path but the orchestration overhead remains visible on large tables.
fix: (left to downstream planner — diagnose-only mode)
verification: (n/a — diagnose-only mode)
files_changed: []
