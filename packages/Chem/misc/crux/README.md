# Crux substructure engine — working notes

Status as of 2026-09-25. How the engine is wired into Chem is in the package `CLAUDE.md`, section
"Crux substructure engine"; this file keeps what is needed to continue: setup, verification, results and the
open items with their analysis.

## Where things are

- `src/crux/` — the integration: `crux-searches.ts` (engine switch, the crux counterpart of
  `chemSubstructureSearchLibrary`), `crux-service.ts` + `crux.worker.ts` (segmented column indexes over a worker
  pool), `crux-smarts.ts` (RDKit query mol → SMARTS crux matches the same way).
- `src/crux/crux_wasm*` — the crux build, three files from one `wasm-pack` run, replaced together:
  - `crux_wasm_bg.wasm` — the engine itself; the only one that lands in `dist/` as a file;
  - `crux_wasm.js` — the JS glue `wasm-bindgen` generates for exactly this binary (bundled into `package.js` and
    the crux worker chunk);
  - `crux_wasm.d.ts` — the glue's TypeScript declarations.
- `src/tests/crux-substructure-search-tests.ts` — package test category `crux substructure search`.
- This folder — the notes, the verification scripts and [`PARITY.md`](PARITY.md), the RDKit-parity ledger (not bundled,
  not published).

The vendored build is crux-core `078ca7c` — branch `spike/rdkit-parity-gaps`, rebased on crux-core main a7e2f26 and
not merged yet —, wasm-pack 0.13.1, Rust 1.95.0. The Chem worker calls `setLenient(true)` on it (`crux.worker.ts`).

## Setting up crux on a machine

- Crux lives in its own GitHub repos. Clone the umbrella `datagrok-ai/crux` and run `bash bootstrap.sh` in it
  (needs GitHub auth: `gh auth login`, `gh auth setup-git`); it clones `crux-core`, `crux-js`, `crux-bench`, … as
  siblings inside it.
- Work in `crux-core` there. It is its own clone (own `.git` and remote), not a submodule — the umbrella ignores
  it — so its `git status` / `diff` show your changes and nothing leaves the machine until you push. Keep the
  sibling layout: crux-js builds from `../crux-core`, crux-core's tools read `../crux-bench/datasets`.
- Open the Claude Code session in `crux-core` itself: its `CLAUDE.md` defines the workflow (`/spike <slug>` →
  branch `spike/<slug>` in a worktree beside the umbrella, `crux-spikes/<slug>`; `/accept` merges it into the local
  `main` and never pushes) and the rule that chemistry fixes are ports of RDKit's C++ (vendored in the repo), not
  rewrites. A spike worktree has no RDKit venv of its own: point `CRUX_RDKIT_PYTHON` at `crux-core/.venv/bin/python`.
- Rust: `rustup` with the toolchain `crux-core/rust-toolchain.toml` pins (1.95.0) and the `wasm32-unknown-unknown`
  target: `rustup-init -y --default-toolchain 1.95.0 --profile minimal -t wasm32-unknown-unknown`. On Windows it
  needs the MSVC build tools (Visual Studio Build Tools).
- Build and vendor (about 1.5 minutes cold). wasm-pack writes a `*` `.gitignore` into its output folder, so never
  point `--out-dir` at Chem:

  ```bash
  cd crux/crux-js
  npx --yes wasm-pack@0.13.1 build ../crux-core/crates/crux-wasm --target web --release --out-name crux_wasm --out-dir ../../../crux-js/src/wasm
  cp src/wasm/crux_wasm.js src/wasm/crux_wasm.d.ts src/wasm/crux_wasm_bg.wasm <reddata>/public/packages/Chem/src/crux/
  ```

  For a spike worktree, build its crate — wasm-pack resolves `--out-dir` against the crate, so the path changes too:

  ```bash
  npx --yes wasm-pack@0.13.1 build ../../crux-spikes/<slug>/crates/crux-wasm --target web --release --out-name crux_wasm --out-dir ../../../../crux/crux-js/src/wasm
  ```

  To check a build before vendoring it, point the scripts below at its output folder with `CRUX_WASM`.

## Verification

Run from this folder, with Node 20 and the workspace installed (`grok setup`); they load Chem's own RDKit build and
transpile Chem's query code (`utils/mol-creation_rdkit.ts`, `crux/crux-smarts.ts`) on the fly, so they check the
code the package runs. Datasets: `<reddata>/data/demo/chem/chembl/chembl-100k.csv` and
`<reddata>/data/demo/chem/smiles_1M.zip` (unzip first); SMARTS sets from `crux-bench/queries/bdb-*`.

- `node parity.mjs [molecules] [sets]` — Chem's RDKit search pipeline against Chem's crux engine on one thread each,
  hit set by hit set. Sets: `smiles`, `molblock`, `molblockH`, `self`, or any crux-bench query set by name (`bdb-3488`,
  `filtercatalog`, `zinc-50`, …, with `CRUX_BENCH=<crux-bench checkout>`). `DATASET=<csv|smi>` to change the data,
  `CRUX_WASM=<wasm-pack out dir>` to test another build, `PURE=1` to drop the RDKit verdict Chem gives to molecules
  crux cannot parse (shows crux's own results), `VERBOSE=1` to print every query. Exits non-zero on any difference.
- `node census.mjs <csv|smi> [molecules]` — molecule by molecule, where crux and Chem's RDKit disagree: A (RDKit reads
  it, crux cannot), B (the reverse), C (crux misses the molecule by RDKit's own canonical SMILES). `OUT=<tsv>` lists
  them; worker threads, ~1 minute for the 1M set.
- `node mem.mjs <csv> [molecules] [segment size]` — memory per indexed molecule, build speed, molecules crux cannot
  parse. The bytes per molecule include the build's high-water mark (wasm memory never shrinks), so measure on a
  big file: on the 1M set it reads 378 bytes.
- `node bench.mjs <csv> [api,filter] [RDKit,Crux]` — in the browser on a running stand with Chem published:
  `grok.chem.searchSubstructure` timings and hit-set hashes per engine, then the substructure filter (time to the
  first result, time to the end, filter updates). `STAND_API` / `STAND_WEB` / `DEV_KEY` default to the local stand
  (`http://localhost:8082`, `http://localhost:8888`, `admin`; on the Mac dev stack both are `http://localhost:8889`, the
  API under `/api`); `STAND_WEB` must serve the client on the API's origin — package web workers cannot start
  cross-origin, so not the pub serve port. `CHROME=<path>` runs an installed Chrome when Playwright's own browser is
  not installed. It leaves the engine property on RDKit.
- Package tests: `grok test --host localhost --category "crux substructure search"`. To run the existing suites
  (`substructure search`, `substructure filters`, `scaffold tree`) on crux, initialize `engineOverride` in
  `crux-searches.ts` to `SubstructureSearchEngine.Crux` for the run — and put it back.

## Results

Parity with Chem's RDKit — molecule by molecule over 42 datasets (1.2M rows) and hit set by hit set over 5,405
queries — is in [`PARITY.md`](PARITY.md): with this build no dataset molecule is unreadable to crux, 7 rows are read
differently (5 of them RDKit 2024.09 vs 2026.03, see open item 1), and every query set gives identical hits.

Package tests on the local stand with this build: 407, all passing but `clone and layout tests` /
`25_clone_layout_scenario`, which fails the same way without these changes. With Crux as the default engine
(`engineOverride`), `substructure search` (all four categories), `substructure filters`, `scaffold tree`,
`chem exported` and `clone and layout tests` (but for the same scenario) pass.

In the browser (`bench.mjs` with this build; 2026-09-25, 18-core arm64 Mac, local stand): all 14
`searchSubstructure` queries give identical hit sets on ChEMBL 100k and on the 1M set, and the filter ends on the same
rows for its 3 queries. (With the previous build, 43f333e, the 1M set gave 13/14: a pyridine N-oxide.)

| | RDKit | Crux |
|---|---|---|
| 100k, first search (fingerprints / index) | 9.8 s | 0.7 s |
| 100k, later searches | 70–890 ms | 3–8 ms |
| 100k filter, first query | 6.4 s, 31 updates | 0.5 s, 3 updates |
| 100k filter, later queries | 0.2 s, 31 updates | 8–11 ms, 1 update |
| 1M, first search / later searches | 94 s / 1.0–16 s | 4.0 s / 6–27 ms |
| 1M filter, first query (RDKit workers already up) | 92 s, 41 updates | 4.3 s, 5 updates |
| 1M filter, later queries | 1.1–1.5 s, 40 updates | 21–24 ms, 1 update |

Crux itself, one thread (`mem.mjs`, 1M set, arm64 Mac): index build ~24 µs per molecule, 375 bytes per molecule;
starting the crux workers takes ~0.2 s (the wasm is compiled once on the main thread and shared).

Done in this round (details in `PARITY.md`):

- crux-core `rdkit-parity-gaps`: `MolOps::cleanUp` (nitro / N-oxides, perchlorates and other halogen oxides, azides),
  elements through Og, RDKit's SMARTS element tokens and whitespace handling, Kekulize aromatic flags, valence lists,
  and a lenient read — RDKit's sanitization without Kekulize, as Chem's `kekulize: false` retry — that the Chem worker
  turns on with `setLenient(true)`. The census's 67 unreadable molecules and 270 of 277 misread ones are gone;
  `rdkitMatch` in `crux-searches.ts` now matters only for rows neither engine reads.
- crux-wasm memory: 981 → ~375 bytes per molecule landed in crux-core main; `MAX_INDEXED_ROWS` stays 3M (~1.1 GB
  over the workers) — not raised, for browsers with less memory.
- Typed radical queries (`[OH]`, `[CH3]`, `[F]`, `[N+]`, `c1cc[n+]cc1`): both engines read them as SMARTS
  (`getQueryMolSafe`); the SMILES reading, a radical, found ~nothing.
- `r` (ring size) goes to RDKit: crux read `[#6;r16,r17,r18]` differently (1,212 hits vs 28).
- A crux search that fails — a worker out of memory, say — restarts crux (freeing its memory) and runs the same search
  on RDKit into the same result (`cruxSubstructureSearch`'s `rdkitSearch`).
- RDKit start-up: in Crux mode the substructure filter no longer starts the RDKit worker pool (it starts at the first
  thing that needs it), and the RDKit wasm is fetched and compiled once on the main thread and handed to every worker
  (`initRdKitFrom`) instead of one download and compile per worker.

## Open items

### 1. RDKit version (Chem)

Chem bundles RDKit MinimalLib 2024.09.1pre; crux ports 2026.03.1. They put Kekulé double bonds differently in some
fused systems — a phosphazene, a fluoranthene, biphenylenes: 5 census rows, and `C=C` / `CCC` / biphenyl queries on
them differ. crux agrees with 2026.03.1 in each (checked with the RDKit Python build). Updating Chem's RDKit closes
them; it is its own change, with the whole Chem test run.

### 2. crux-core follow-ups (no Chem dataset needs them)

- Charged-atom valence: crux lacks RDKit's strict valence check and accepts `C[O-2]C`, `C[N-2]C`,
  `C[S+2](C)(C)(C)C`, which RDKit rejects (census: 0 such molecules). Porting RDKit's valence model
  (`Atom::calcExplicitValence` with the effective atomic number) would also close most of the aromaticity fuzz's
  ~48k soft mismatches — crux-core `research/rdkit-parity-gaps/README.md`, "The fuzz shift".
- SMARTS grammar: 15 FilterCatalog queries run on RDKit because crux cannot parse their translation (a `$()` holding
  branches inside a bracket or another `$()`; `;` in a bond expression).
- Oligonucleotides: crux misses a ~320-atom oligonucleotide by its own SMILES; it finds its fragments up to 300 atoms
  (1–4 s each) and misses at 310 — the matcher on large repetitive queries. Only whole-molecule queries that big hit it.
- Dative bonds (`cleanUpOrganometallics`): 4 Cu / Zn complexes of the 1M set read differently in crux-core's own
  parse parity; not visible in Chem's census.

### 3. Typed SMARTS with explicit hydrogens (a Chem semantics question, both engines)

A typed SMARTS written with hydrogen atoms (`c-[#8]-[#1]`, PAINS style — 190 of the 1,585 FilterCatalog SMARTS) finds
nothing on ordinary columns in either engine: RDKit matches `[#1]` only against explicit hydrogens, and Chem reads
typed SMARTS as they are (molblocks with H atoms get `mergeQueryHs`). RDKit's FilterCatalog merges them into H counts
(`c-[#8]-[#1]` → phenol). Chem's MinimalLib `get_qmol` takes no options, so merging would go through a molblock or a
newer RDKit. Leaving it as is keeps both engines in agreement.

## Gotchas met on the way

- `grok build` regenerates `src/package.g.ts`, `src/package-api.ts` and `src/generated/db.ts`, picking up unrelated
  upstream changes; revert them unless they are yours.
- The workspace build can fail in `@datagrok-libraries/compute-utils` (`json-logic-js` missing) — unrelated; Chem
  still bundles.
- `col.toList()` leaves holes for empty cells, and `Array.prototype.map` skips holes: use `Array.from(list, fn)`.
- A debug `grok publish` lands as the publishing user's version; the stand's other users keep the released one.
- Chem's ESLint is 1TBS (`} catch {` on one line).
- `grok test` runs a Playwright pass after the Puppeteer one and reports "Tests failed" when Playwright's browser is
  not installed: add `--skip-playwright`. `--category` matches by prefix (`substructure search` runs its three
  sub-categories too).
