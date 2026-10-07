# HMM

Profile hidden Markov models in the browser, powered by an exact
Rust/WebAssembly port of [HMMER 3.4](http://hmmer.org). Nothing leaves the
browser: sequences are compared in Web Workers, and results are identical to
the C programs (verified against the canonical Linux build, see "Verification").

## Antibody and T-cell receptor numbering (ANARCI)

`Bio | Annotate | Apply Numbering Scheme...` gains the **ANARCI (HMMER)**
engine next to Bio's own. It implements
[ANARCI](https://github.com/oxpig/ANARCI) (Dunbar & Deane, 2016) on the
HMMER engine:

* schemes: IMGT, Kabat, Chothia, Martin, AHo, Wolfguy;
* chains: heavy, kappa, lambda and TCR alpha, beta, gamma, delta (TCRs in IMGT
  and AHo, as ANARCI);
* the aligned column, FR/CDR regions and per-row annotations, as with any
  numbering engine;
* scFv and other multi-domain chains: the N-terminal domain is numbered.

`Bio | Annotate | Germlines and Species (ANARCI)...` adds the chain type,
species, closest V and J germline genes with their identities, and the HMMER
E-value and bit score. Species default to ANARCI's preference for human and
mouse; pick "all" or one species to change it.

Results are identical to ANARCI 2024.05.21 run on the same germline data
(numbering, germlines, E-values), including for long CDR3s and N/C-terminal
extensions. The germline models and tables are rebuilt from current IMGT® data
(October 2026) with ANARCI's own build pipeline.

## Profile HMM search

`Bio | Search | Profile HMM Search...` searches a sequence column with one
profile HMM (`hmmsearch`): a family of the built-in Pfam library, any Pfam
family by accession (fetched from InterPro), or an HMM file (HMMER2/3 text or
`hmmpress` `.h3m`, optionally gzipped). It adds score, E-value, domain count,
significance (inclusion threshold) and domain region columns, and annotates the
domains on the sequences. E-values treat the column as the database (Z = the
number of non-empty rows). Pfam gathering thresholds (`--cut_ga`) are
available for models that have them.

## Protein domains

`Bio | Annotate | Find Domains (HMMER)...` annotates domains (`hmmscan`) from
the built-in Pfam library of antibody- and biologics-relevant families, Pfam
accessions or an HMM file, with Pfam gathering thresholds or E-values. Domains
are drawn as colored regions on the sequences, with a summary column.

## Scripting

* `HMM:anarciNumbering(df, seqCol, scheme)` — the numbering engine (Bio's
  five-column contract plus chain, species, germline and E-value columns).
* `HMM:searchWithHmm(table, sequence, model, evalue, gathering)` — model is
  HMMER text, Pfam accession(s) or `builtin:<index>`.
* `HMM:findDomains(table, sequence, library, evalue, gathering)` — library is
  `builtin`, Pfam accessions or HMMER text; returns one row per domain.

## Performance

The engine is a 360 KB WebAssembly module (145 KB gzipped) compiled with SIMD,
loaded on first use and shared by up to eight workers. On a laptop, numbering
1,000 antibody chains against all 29 ANARCI germline models takes about two
seconds. Size and speed budgets are enforced by the engine's test harness.

## Verification

The engine and the ANARCI port live in the Rusty-HMMER repository
(`crates/hmmer-web`, `web/`), which verifies them against the C HMMER 3.4
oracle (hmmsearch/hmmscan tabular outputs, including searches split across
workers) and against unmodified ANARCI run in Docker (numbering, germlines and
E-values for 2,500 sequences × 6 schemes × 2 species settings). The package
tests repeat a sample of those comparisons in the browser.

## Credits

The engine is a Rust translation of HMMER 3.4 and Easel (BSD-3-Clause,
BSD-2-Clause) with math routines from Arm Optimized Routines (MIT); the
numbering is a TypeScript translation of ANARCI (BSD-3-Clause). The germline
models and tables are derived from IMGT® data (CC BY 4.0); the domain library
is Pfam 37.0 (CC0). See [CREDITS.md](CREDITS.md) for the full attribution and
`licenses/` for the license texts.
