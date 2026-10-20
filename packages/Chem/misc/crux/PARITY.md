# Crux vs RDKit — parity ledger

What Chem's crux engine gives against Chem's RDKit engine, before and after crux-core's `rdkit-parity-gaps` spike and
the Chem changes that came with it (2026-09-25, arm64 Mac, Node 20). RDKit here is Chem's own: MinimalLib
2024.09.1pre, read the way Chem reads molecules (`getMolSafe`) and queries (`getQueryMolSafe`).

| | crux build | Chem |
|---|---|---|
| before | crux-core main afc192f (the census; the vendored 43f333e gives identical results) | public 2562ecd |
| after | crux-core `spike/rdkit-parity-gaps` c3230f3, `setLenient(true)`; the vendored 078ca7c is the same branch rebased on crux-core main a7e2f26 and gives the same census, row for row, and the same search results | this change |

crux-core's own gates for the same steps (RDKit 2026.03.1 oracle, parse parity, rk-parity, fuzz) are in crux-core
`research/rdkit-parity-gaps/README.md` on that branch.

## Molecules — `census.mjs`

42 datasets, 1,222,093 rows: Chem's `files/` csvs, `data/demo/chem/*.csv`, ChEMBL 100k, the 1M set. Classes:

- A — RDKit reads the molecule, crux cannot parse it;
- B — crux parses it, RDKit cannot;
- C — both read it, but crux does not find the molecule with RDKit's canonical SMILES as the query.

| | before | after |
|---|---|---|
| A | 67 | 0 |
| B | 0 | 0 |
| C | 277 | 7 |

| dataset | rows | A | C |
|---|---|---|---|
| smiles_1M | 1,000,000 | 25 → 0 | 217 → 6 |
| smiles_50K | 50,000 | | 16 → 1 |
| smiles_20K | 20,000 | | 6 → 0 |
| smi10K, smiles_10K_with_activities | 10,000 each | | 1 → 0 each |
| mmp_demo | 20,267 | | 2 → 0 |
| solubility_train | 1,025 | 4 → 0 | 34 → 0 |
| pmpo | 663 | 35 → 0 | |
| stereochemistry, Test_smiles_malformed, Test_smiles_with_empty_and_malformed | | 1 → 0 each | |

What they were and what closed them (crux-core commits):

| before | what | closed by |
|---|---|---|
| 64 A | aromatic rings RDKit reads only without Kekulize (Chem's `kekulize: false` retry): a missing `[nH]` (`c4ncnc4C`), aromatic P phosphazenes, a 7-atom aromatic ring | `parse_smiles_lenient` + `setLenient(true)` |
| 2 A | radium salts `[Ra+2]` | elements through Og |
| 1 A | leading spaces | RDKit's whitespace handling |
| 227 C | halogen and P oxides (mostly perchlorate salts `[O-]Cl(=O)(=O)=O`, `OCl=O`), an azide, ring N-oxides `n4=O` | `MolOps::cleanUp` |
| 34 C | nitro / N-oxides written `N(=O)=O` | `MolOps::cleanUp` |
| 4 C | pyridine N-oxides written `n(=O)` | `MolOps::cleanUp` |
| 5 C | aromatic P that RDKit reads only without Kekulize | Kekulize aromatic-bond flags, then the lenient read |

Counts are rows; a molecule in two datasets counts twice (the 227 are 204 distinct molecules).

The 7 left:

- 5 — a phosphazene (in the 1M and 50K sets), a fluoranthene, two biphenylenes: RDKit 2026.03.1, the version crux
  ports, reads them as crux does; Chem's 2024.09 puts the Kekulé double bonds elsewhere
  (`Cn1nc(N)c2c(np(nc12)…)c5ccccc5` is `…N=P(…)N=C2…` in 2026.03.1 and `…=NP(…)=NC=2…` in 2024.09). They go when
  Chem's RDKit is updated.
- 2 — oligonucleotides of ~320 atoms. crux finds their fragments up to 300 atoms (1–4 s per query) and misses at
  310: the matcher on a large repetitive query (crux-core).

## Searches — `parity.mjs`

Chem's RDKit pipeline against Chem's crux pipeline, hit set by hit set (`queries` / `differ` / `RD` = run on RDKit):

| set | molecules | queries | before: differ, RD | after: differ, RD |
|---|---|---|---|---|
| Chem's sets (SMILES, molblock, explicit H, whole molecules) | ChEMBL 100k | 282 | 0, 9 | 0, 4 |
| bdb-3488 (crux-bench) | ChEMBL 20k | 3,488 | 0, 121 | 0, 0 |
| FilterCatalog (crux-bench) | ChEMBL 20k | 1,585 | 3, 292 | 0, 85 |
| zinc-50 (crux-bench) | ChEMBL 100k | 50 | 0, 0 | 0, 0 |
| census rows, SMILES + whole molecules, `PURE=1` | the 319 census molecules | 122 | 84, 7 | 3, 2 |
| the same with Chem's RDKit verdict for unparsed molecules | the 319 census molecules | 122 | 40, 7 | 3, 2 |

In the browser, on the local stand (`bench.mjs`: `grok.chem.searchSubstructure` and the substructure filter): 14 of 14
queries give identical hit sets on ChEMBL 100k and on the 1M set (13 of 14 on the 1M set before: the pyridine
N-oxide), and the filter ends on the same rows for its 3 queries.

- FilterCatalog before: `[#6;r16,r17,r18]…` (twice; crux 1,212 hits, RDKit 28 — crux reads `r` differently, now sent
  to RDKit) and a metal list crux misread through `[se]` (RDKit's SMARTS element tokens now).
- Still on RDKit, by design: `r`, `Rn`, `v`, `x`, `h`, `z`, `^`, isotopes, radicals in molblocks (`M  RAD`); and 15
  FilterCatalog SMARTS whose translation crux cannot parse — 13 with a `$()` holding branches inside a bracket or
  another `$()`, 2 with `;` in a bond expression (`-,=;!@`): crux-core SMARTS grammar. (Raw, 82 of the 1,585 do not
  parse in crux; the translation writes most of them, e.g. unbracketed `a` as `[a]`, in a form crux reads.)
- The 3 census-row differences are the RDKit-version molecules above: `CCC` in the phosphazene, `C=C` in it and the
  biphenylenes, biphenyl in the biphenylenes — crux gives RDKit 2026.03.1's answer in each.

### Typed radical queries

A typed query whose SMILES reading has unpaired electrons (`[OH]`, `[CH3]`, `[F]`, `[N+]`, `c1cc[n+]cc1`) is now read
as SMARTS in both engines (`getQueryMolSafe`). Measured on 20k ChEMBL rows, RDKit's hits change only where the reading
does: `[OH]` and the like already read as SMARTS whenever the SMARTS matched its own molecule; `c1cc[n+]cc1`, which the
radical made a non-aromatic ring, goes 0 → 72 hits; 190 FilterCatalog SMARTS written with explicit `[#1]` keep their 0
hits (RDKit matches `[#1]` only against explicit hydrogens, in either reading) and now run on crux.

## Memory — `mem.mjs`, the 1M set, 33k-molecule segments

| crux build | bytes per molecule | molecules not parsed | index build |
|---|---|---|---|
| main afc192f | 372 | 25 | 23.3 µs/mol |
| c3230f3, lenient | 378 | 0 | 23.6 µs/mol |
| 078ca7c (rebased), lenient | 375 | 0 | — (measured under load) |

Benzene hits 837,395 → 837,412 and `[#7]` 935,527 → 935,550 on the set: the molecules crux could not parse before.
`MAX_INDEXED_ROWS` stays 3M (~1.1 GB over the workers).
