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
- This folder — the notes and the verification scripts (not bundled, not published).

The vendored build is crux-core `43f333e` (main), wasm-pack 0.13.1, Rust 1.95.0.

## Setting up crux on a machine

- Crux lives in its own GitHub repos. Clone the umbrella `datagrok-ai/crux` and run `bash bootstrap.sh` in it
  (needs GitHub auth: `gh auth login`, `gh auth setup-git`); it clones `crux-core`, `crux-js`, `crux-bench`, … as
  siblings inside it.
- Work in `crux-core` there. It is its own clone (own `.git` and remote), not a submodule — the umbrella ignores
  it — so its `git status` / `diff` show your changes and nothing leaves the machine until you push. Keep the
  sibling layout: crux-js builds from `../crux-core`, crux-core's tools read `../crux-bench/datasets`.
- Open the Claude Code session in `crux-core` itself: its `CLAUDE.md` defines the workflow (`/spike <slug>` →
  branch `spike/<slug>` in a worktree under `crux/crux-spikes/`; `/accept` merges it into the local `main` and
  never pushes) and the rule that chemistry fixes are ports of RDKit's C++ (vendored in the repo), not rewrites.
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

  For a spike worktree, pass its `crates/crux-wasm` path instead of `../crux-core/crates/crux-wasm`. To check a
  build before vendoring it, point the scripts below at its output folder with `CRUX_WASM`.

## Verification

Run from this folder, with Node 20 and the workspace installed (`grok setup`); they load Chem's own RDKit build and
transpile Chem's query code (`utils/mol-creation_rdkit.ts`, `crux/crux-smarts.ts`) on the fly, so they check the
code the package runs. Datasets: `<reddata>/data/demo/chem/chembl/chembl-100k.csv` and `smiles_1M.zip` next to it
(unzip first); SMARTS sets from `crux-bench/queries/bdb-*`.

- `node parity.mjs [molecules] [sets]` — Chem's RDKit search pipeline against Chem's crux engine on one thread each,
  hit set by hit set. Sets: `smiles`, `molblock`, `molblockH`, `self`, `bdb-100|bdb-1k|bdb-3488` (with
  `CRUX_BENCH=<crux-bench checkout>`). `DATASET=<csv>` to change the data, `CRUX_WASM=<wasm-pack out dir>` to test
  another build, `PURE=1` to drop the RDKit verdict Chem gives to molecules crux cannot parse (shows crux's own
  results — use it for open items 1 and 2). Exits non-zero on any difference.
- `node mem.mjs <csv> [molecules] [segment size]` — memory per indexed molecule, build speed, molecules crux cannot
  parse. The bytes per molecule include the build's high-water mark (wasm memory never shrinks), so measure on a
  big file: on 20k molecules in one segment it reads ~1.6 KB, on the 1M set ~1 KB.
- `node bench.mjs <csv> [api,filter] [RDKit,Crux]` — in the browser on a running stand with Chem published:
  `grok.chem.searchSubstructure` timings and hit-set hashes per engine, then the substructure filter (time to the
  first result, time to the end, filter updates). `STAND_API` / `STAND_WEB` / `DEV_KEY` default to the local stand
  (`http://localhost:8082`, `http://localhost:8888`, `admin`); `STAND_WEB` must serve the client on the API's origin —
  package web workers cannot start cross-origin, so not the pub serve port. It leaves the engine property on RDKit.
- Package tests: `grok test --host localhost --category "crux substructure search"`. To run the existing suites
  (`substructure search`, `substructure filters`, `scaffold tree`) on crux, initialize `engineOverride` in
  `crux-searches.ts` to `SubstructureSearchEngine.Crux` for the run — and put it back.

## Results (2026-09-25, 32-thread machine, local stand)

Parity: `parity.mjs` over 3,730 queries (all 3,488 bdb SMARTS plus SMILES, molblock, explicit-H and whole-molecule
queries) on ChEMBL 20k–100k: no difference; 130 of them run on RDKit, all radicals (item 3) but one typed `Se`,
which is no molecule. In the browser, ChEMBL
100k: 14/14 queries identical; the 1M set: 13/14 (item 1). Package tests pass, and the existing substructure
suites pass with either engine as the default.

| | RDKit | Crux |
|---|---|---|
| 100k, first search (fingerprints / index) | 24.8 s | 0.9 s |
| 100k, later searches | 150–590 ms | 6–11 ms |
| 100k filter, first query | 4.8–21 s, 30 updates | 0.6–1.0 s, 3–4 updates |
| 100k filter, later queries | 0.3–0.6 s, 30 updates | 10–16 ms, 1 update |
| 1M, first search / later searches | 97 s / 1.2–14 s | 7.5 s / 13–50 ms |
| 1M filter, first query (RDKit workers already up) | 88 s, 40 updates | 9.5 s, 5 updates |

Crux itself, one thread: index build ~70 µs per molecule, ~981 bytes per molecule (item 4); starting the crux
workers takes ~0.2 s (the wasm is compiled once on the main thread and shared).

## Open items

### 1. crux-core: pyridine N-oxide written `n(=O)` is not aromatic

crux parses these molecules but does not perceive the ring as aromatic; RDKit does, so `c1ccncc1` finds them in
RDKit only. Four rows of the 1M set:

```
Oc1c2c(CCc3ccccc23)nn1c4ccccn4=O
CC(C)(C)c1ccc(cc1)S(=O)(=O)Nc2ccc(Cl)cc2C(=O)c3ccn(=O)cc3
[Na+].CC(C)(C)c1ccc(cc1)S(=O)(=O)[N-]c2ccc(Cl)cc2C(=O)c3ccn(=O)cc3
COc1ccc(cc1)S(=O)(=O)N(C(=O)C)c2ccccc2C=Cc3ccn(=O)cc3
```

A crux-core parity bug (aromaticity perception); Chem has no workaround. `PURE=1 node parity.mjs` with these as the
dataset shows it.

### 2. crux-core: molecules crux cannot parse

25 rows of the 1M set; Chem currently gives them RDKit's verdict (`rdkitMatch` in `crux-searches.ts`, on the main
thread) and logs them once per column ("crux could not parse N molecules"). Without that, crux misses 17 benzene,
22 halogen and 6 amine hits on the set.

- 23 phosphazenes with aromatic phosphorus, e.g.
  `N1c2ccccc2p3(c4ccccc14)c5ccccc5nc6ccccc36`, `FC1(F)COp2(OCC1(F)F)np(Cl)(Cl)np(Cl)(Cl)n2`,
  `S=P1(NNp2(NN1)np(np(n2)(N3CC3)N4CC4)(N5CC5)N6CC6)Oc7ccccc7`. RDKit's default SMILES parser rejects them too
  (kekulization fails); Chem accepts them only through `getMolSafe`'s retry with `kekulize: false`. The same holds
  for `COc1cc(C)ccc1OC(=O)C(C)c2cccc(CC(C)C)cc2` (a 7-atom aromatic ring) in
  `files/tests/Test_smiles_with_empty_and_malformed.csv`. So in crux this is an opt-in lenient parse, not an RDKit
  parity fix.
- 2 radium salts, `[Cl-].[Cl-].[Ra+2]` and `[Cl-].[Cl-].[223Ra+2]`: RDKit's default parser accepts them, crux does
  not — a plain crux gap.

### 3. Radical queries (a Chem semantics decision, both engines)

What happens today, with the RDKit engine:

- For a typed query that is not a molblock, `getQueryMolSafe` (`utils/mol-creation_rdkit.ts`) decides SMILES or
  SMARTS with `_isSmarts`, a regex looking for `#<digit>`, `$`, `&`, `;`, `,` or `!`. `[OH]`, `[NH2]`, `[NH]`,
  `[CH3]`, `[CH2]`, `[F]`, `[Cl]`, `[N+]`, `[O-]` have none, so they are parsed as SMILES and that molecule is the
  query.
- In SMILES a bracket atom's H count is exact and free valence becomes unpaired electrons: `[OH]` is a hydroxyl
  radical, `[CH3]` a methyl radical, `[F]` a fluorine atom, `[N+]` a radical N+.
- RDKit matches a query atom's unpaired electrons exactly when it has any, so only radicals match: 0 hits on
  ordinary data (ChEMBL 100k: `[OH]` finds 26, real radicals such as nitroxides). Read as SMARTS, as almost
  certainly meant, a 5k ChEMBL sample gives 1,995 hits for `[OH]`, 965 for `[F]`, 116 for `[N+]`.
- Scope: typed queries only (API, filter operators, typed SMILES) — 121 of the 3,488 bdb queries. A sketcher
  molblock carries a radical only when the user sets one (`M  RAD`), and then the radical meaning is right.
- crux cannot reproduce it: SMARTS has no radical primitive and crux stores no unpaired electrons, so
  `getCruxSmarts` returns null for these and they run on RDKit.

Options:

1. Keep: identical results, these queries run on RDKit and find ~0.
2. Recommended: read such typed queries as SMARTS in both engines. In `getQueryMolSafe`'s non-molblock branch, when
   the SMILES reading has unpaired electrons and a SMARTS reading exists, keep the SMARTS one (about three lines;
   leave `_isSmarts` alone — it also decides which context panels show). Both engines then return the SMARTS
   result and crux takes these queries (`[OH]` → `[O&H1]`). It changes RDKit results for them (~0 → real hits) and
   applies wherever `getQueryMolSafe` is used: the RDKit search, the crux translation, the scaffold tree, R-group
   core parsing. Molblock queries keep radical semantics.
3. Drop radicals only in the crux translation: crux would then disagree with RDKit on exactly these queries.
4. Radicals in crux: no SMARTS syntax for it — not feasible.

### 4. crux-wasm memory (crux-core)

- In `crates/crux-wasm/src/lib.rs`, `CollectionBuilder` keeps every molecule's parsed graph (`MolProps`) in
  `MolStore.props`. `finish()` encodes them into the CXMOL blob the matcher reads, then moves the store — graphs
  included — into the `Collection`.
- An indexed collection never reads them again: search reads molecules from the CXMOL index, screening from its
  fingerprint sidecar, similarity from the ECFP4 sidecar. Only `len()` and the non-indexed "direct" mode use them;
  Chem builds indexed collections only.
- The change, tested on a copy of crux-core (~5 lines; direct mode unchanged):

  ```diff
       let (blob, _report) = build_blob_from_props(&store.props);
       let idx = RkIndex::open(blob);
  +    let mut store = store;
       if build_index {
  +        // the indexed path never reads the parsed graphs again
  +        store.props = Vec::new();
           let rarity = build_rarity_table(screen_fps.iter().map(|v| v.as_slice()));
  ...
       pub fn len(&self) -> usize {
  -        self.store.props.len()
  +        self.store.orig_input_ids.len()
       }
  ```

- Measured on the 1M set (`mem.mjs`, 33k-molecule segments): 981 → 377 bytes per molecule, identical hits, build
  74.6 → 66.2 µs per molecule, searches unchanged.
- What remains: CXMOL ~144 bytes, the 1024-bit screening fingerprint ~130, the 256-bit ECFP4 ~35, input ids 4, and
  the build high-water mark. Smaller follow-ups: a builder option to skip the ECFP4 sidecar (Chem does not use crux
  similarity: ~35 bytes per molecule and some build time), and a flat buffer for the screening fingerprints held
  until `finish()` (lower build peak).
- Chem after it lands: vendor the new build and raise `MAX_INDEXED_ROWS` in `crux-service.ts` from 3M
  (≈3 GB at 1 KB per molecule) to about 8M.

### 5. RDKit worker pool start-up

- Measured on the local stand: opening the substructure filter starts 31 downloads of the 5.7 MB RDKit wasm (main
  thread plus one per worker, 30 workers). The last ended at 18.3 s, and crux's small worker script loaded only
  after them (18.4–18.6 s): the dev nginx speaks HTTP/1.1, so the browser opens 6 connections to it and 171 MB
  squeeze through them; each worker also compiles its own copy.
- The only eager start is the filter constructor, `initRdKitService(); // No await`
  (`widgets/chem-substructure-filter.ts`). Everything else calls `getRdKitService()`, which starts the pool on
  first use, once.
- Planned fix: skip that line when the engine is Crux. The pool then starts at the first thing that needs it:
  - a search crux hands to RDKit (Included in, Exact, Stereo agnostic, Similar, or a query crux cannot express);
  - a crux search on a molblock column — the workers convert molblocks to SMILES — so there the first crux search
    still waits for the pool;
  - features on the workers: R-groups, MMP, MCS, SAR matrix, BitBIRCH, reactions, the structural alerts /
    pharmacophore / InChI panels, SDF export, recalculate coordinates, flatten / beautify, notation conversion.

  Grid rendering and crux's query translation use the main-thread RDKit, which is loaded anyway.
- Trade-off: whatever needs RDKit first pays the pool start then (~17 s on the local stand) and, for a search,
  RDKit's first fingerprint pass on the column (~25 s on 100k), which RDKit-mode users pay today. After switching
  the property back to RDKit, the pool starts at the first search, as now.
- Complementary fix, also for RDKit mode: compile the RDKit wasm once on the main thread and give the module to
  the workers — what cut crux's worker start from 1.7 s to 0.2 s. RDKit's build supports it through its
  `instantiateWasm` option; it removes the 30 downloads and 30 compiles and touches `rdkit.worker.ts` (module
  init), `RdKitServiceWorkerClient.moduleInit` and the main-thread init in `utils/chem-common-rdkit.ts`.
- On production, HTTP/2 and caching should make the delay smaller (not measured); the 30 compiles remain.

## Gotchas met on the way

- `grok build` regenerates `src/package.g.ts`, `src/package-api.ts` and `src/generated/db.ts`, picking up unrelated
  upstream changes; revert them unless they are yours.
- The workspace build can fail in `@datagrok-libraries/compute-utils` (`json-logic-js` missing) — unrelated; Chem
  still bundles.
- `col.toList()` leaves holes for empty cells, and `Array.prototype.map` skips holes: use `Array.from(list, fn)`.
- A debug `grok publish` lands as the publishing user's version; the stand's other users keep the released one.
- Chem's ESLint is 1TBS (`} catch {` on one line).
