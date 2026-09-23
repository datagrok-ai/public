# Crux filter demo

Compiled `@datagrok/crux-js` artifacts, including the Rust WASM binary. Exact source
revisions and the upstream license are recorded in `provenance.json`. Source lives
in the private `datagrok-ai/crux-js` and `datagrok-ai/crux-core` repositories.

To refresh, build `crux-js` with `npm run build`, then run from Chem:

```sh
node scripts/update-crux.mjs /path/to/crux/crux-js
grok build
grok publish localhost --skip-build --skip-docker-rebuild
```

Rspack bundles the ESM worker entry and emits the WASM file into Chem's `dist/`.
No separate Crux server or Vite server is needed. Reload Datagrok after publishing.

The molecule filter's **Contains** mode uses Crux. Other modes, radical molecule
queries (whose radical constraints are lost in SMARTS serialization), and the public
search APIs still use RDKit. RDKit also supplies query normalization, molblock
conversion, sketching, rendering and highlighting. Each active filter lazily owns
up to eight workers (at most half the available cores), indexes its column in 25,000-row shards, caches the index
between sketches, rebuilds it after column edits, and disposes its workers on
detach. Query cancellation takes effect between shards.

Result updates use progressively spaced milestones (1, 5, 10, 20, 40 and 70%),
emitting the actual completed percentage when a milestone is crossed, and a
single final 100%. Intermediate updates are at least 150 ms apart; milestones
crossed during that interval are combined. There are at most seven updates per
search, usually fewer for fast searches. Every hit is still accumulated.

Chem snapshots the column's category codes instead of expanding it through
`toList()`, normalizes each distinct string once, and loads one 25,000-molecule
chunk per shard. Opening the filter no longer eagerly starts RDKit workers;
RDKit-only modes and molblock conversion initialize them when needed.

To try it at `http://localhost:8889`, open a molecule table, open Filters, and
sketch a benzene ring in the molecule filter with **Contains** selected. Chrome's
worker inspector shows the Crux worker entry (`crux-pool`).

Validation:

```sh
grok test --host localhost --category 'Crux substructure' --no-retry --csv --skip-playwright
```

This is a local integration demo, not a claim of complete RDKit query parity or
a measured speedup. First use includes parsing/index construction. Compare warm
queries separately when benchmarking.

Crux skips unparseable/unsupported target molecules while preserving their row
positions. In Chem's malformed SMILES fixture, it rejects the non-kekulizable
seven-member aromatic ring at row 13, which RDKit's permissive fallback accepts:
benzene therefore gives 35 hits instead of 36. This difference is recorded in the
filter regression test; it is not an indexing offset error.

## Local validation (2026-09-23)

* `pnpm run build` succeeds. The new service and test file pass ESLint.
* `Crux substructure`: 7/7 browser tests pass, including consecutive cell edits.
* `substructure filters`: 12/12 existing tests pass with the malformed-input
  expectation described above.
* Live filter at localhost:8889: `tests/smi10K.csv` gives 8,970 benzene hits
  from 10,001 rows. `tests/smiles_50K.csv` gives 45,199 benzene hits and 7,511
  pyridine hits from 50,000 rows. Four Crux workers start, two shards load once,
  and the second query reuses those shards without loading again.
* The 10K check used the ribbon's Filters button and sketcher dialog. On the 50K
  table the panel opened empty and the add-column control did not add a card;
  `FilterGroup.updateOrAdd` added it before both queries were entered through
  the sketcher UI.
* Closing the filter viewer terminates all four Crux workers and restores all
  50,000 rows. The check used `FilterGroup.close()` after a ribbon menu intercepted
  clicks on the close button.
* Package-wide type checking still reports seven existing buffer-type errors
  in `docker/api.ts`, MMP and substructure search tests, and the filter's existing
  `BitSet.fromBytes` calls. Existing files also have unrelated lint errors.

## Progress and startup tuning (2026-09-23)

The million-row regression uses the 50,000-row fixture repeated twenty times.
On the local 18-core machine, the original first search took 8.30 s and emitted
41 updates (including duplicate completion). With the Chem-side preparation
changes and eight workers it took 4.50 s and emitted two updates, at 2.5% and
100%. A subsequent cached query took 19 ms. Results stayed at 903,980 benzene
hits and 150,220 pyridine hits. These are local observations, not speed guarantees
for other datasets or machines. Fingerprint prescreening remains enabled.

After tuning, all eight `Crux substructure` tests and all twelve
`substructure filters` tests pass. The build and lint checks for the new service
and tests pass; the existing package-wide type errors listed above remain.
