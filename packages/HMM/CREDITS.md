# HMM — Third-Party Software and Data

The `@datagrok/hmm` package is distributed under the MIT license that covers the
rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)). It
incorporates the software and data listed below; this file reproduces the
attribution and notices their licenses require. Full license texts are in
[`licenses/`](licenses/).

Everything bundled is under a permissive license (BSD-3-Clause, BSD-2-Clause,
MIT, Unicode-3.0) or an open data license (CC BY 4.0, CC0). No copyleft
(GPL/LGPL/MPL) component is bundled. The names of the authors and institutions
below are used for attribution only; they do not endorse this package.

---

## 1. Software bundled in the published package

### HMMER 3.4 — translated to Rust, compiled to WebAssembly

`wasm/hmmer_web.wasm` (copied to `dist/`) is compiled from a Rust translation of
HMMER 3.4: profile configuration, the MSV, bias, Viterbi, Forward and Backward
algorithms, domain definition, null2 correction, alignments, E-value statistics,
hit ranking and thresholds, and the HMM file readers. The translation was made
in the Rusty-HMMER project and verified to give the same results as the
original C programs. It is a modified (translated) version of HMMER.

- Upstream: http://hmmer.org — https://github.com/EddyRivasLab/hmmer
- License: **BSD-3-Clause**, full text in [`licenses/HMMER.txt`](licenses/HMMER.txt)
- Copyright (C) 1992–2023 Sean R. Eddy; (C) 2015–2023 President and Fellows of
  Harvard College; (C) 2000–2023 Howard Hughes Medical Institute; (C) 1995–2006
  Washington University School of Medicine; (C) 1992–1995 MRC Laboratory of
  Molecular Biology.
- Reference: Eddy SR. Accelerated profile HMM searches. PLoS Comput Biol
  7(10):e1002195 (2011).

### Easel 0.49 — translated to Rust, in the same engine

The parts of Easel that HMMER depends on (alphabets and digitization, sequence
and model input, the random number generator, statistical distributions, score
matrices, vector routines).

- Upstream: https://github.com/EddyRivasLab/easel
- License: **BSD-2-Clause**, full text in [`licenses/Easel.txt`](licenses/Easel.txt)
- Copyright (C) 1990–2023 Sean R. Eddy; (C) 2015–2023 President and Fellows of
  Harvard College; (C) 2000–2023 Howard Hughes Medical Institute; (C) 1995–2006
  Washington University School of Medicine; (C) 1992–1995 MRC Laboratory of
  Molecular Biology, UK.

### Arm Optimized Routines — math functions in the engine

`expf`, `logf`, `log2f`, `exp`, `log` and `log2`, translated to Rust from
commit `fe09d4e7ed62aa24c81b5791c810bc9de70006e9`, so that results match the
reference Linux build bit for bit.

- Upstream: https://github.com/ARM-software/optimized-routines
- License: MIT OR Apache-2.0 WITH LLVM-exception; used under **MIT**. Full
  text in [`licenses/ARM-optimized-routines.txt`](licenses/ARM-optimized-routines.txt)
- Copyright (c) 2017–2025 Arm Limited.

### Rust Standard Library and runtime crates

The engine statically links parts of the Rust Standard Library (`core`,
`alloc`, `std`, with `compiler-builtins`) and the runtime crates `dlmalloc`,
`hashbrown` and `cfg-if`.

- License: MIT OR Apache-2.0, used under **MIT**; Unicode tables in `core`
  under **Unicode-3.0**. Notices in
  [`licenses/Rust-standard-library.txt`](licenses/Rust-standard-library.txt)
- Copyright (c) The Rust Project Developers; dlmalloc and cfg-if: Alex Crichton;
  hashbrown: Amanieu d'Antras; Unicode data: (c) 1991–2024 Unicode, Inc.

### ANARCI 2024.05.21 — translated to TypeScript

`src/hmmer/anarci/*.ts` translate ANARCI's numbering method to TypeScript:
the parsing of hmmscan hits into IMGT state vectors and `check_for_j` from
`anarci.py`, the numbering schemes of `schemes.py` (IMGT, Kabat, Chothia,
Martin, AHo, Wolfguy) and germline assignment. Modifications: translated from
Python; a numbering error fails one sequence instead of the whole batch.

- Upstream: https://github.com/oxpig/ANARCI (commit
  `79f6c575056dedef86cb8f405ebb039197923eec`)
- License: **BSD-3-Clause**, full text in [`licenses/ANARCI.txt`](licenses/ANARCI.txt)
- Copyright 2019 Charlotte Deane, James Dunbar, Alexsandr Kovaltsuk, Claire
  Marks; `schemes.py` Copyright (C) 2016 Oxford Protein Informatics Group.
- Reference: Dunbar J, Deane CM. ANARCI: antigen receptor numbering and receptor
  classification. Bioinformatics 32(2):298–300 (2016).

---

## 2. Data bundled in the published package

### IMGT® germline gene sequences — ANARCI germline models and tables

`wasm/anarci/ALL.hmm.h3m.gz` (29 profile HMMs: heavy, kappa and lambda chains
of eight species, TCR alpha, beta, gamma and delta chains of human and mouse)
and `wasm/anarci/germlines-*.json`, `species.json` and `hmm-lengths.json` are
derived from IMGT/GENE-DB germline V and J gene sequences, downloaded from
https://www.imgt.org/genedb/ on 7 October 2026 (15:26 UTC).

Modifications: ANARCI's build pipeline (commit above) extracted the gapped
V-gene and J-gene amino acid sequences, aligned the J genes with MUSCLE,
assembled IMGT-numbered alignments per species and chain, built profile HMMs
with HMMER 3.4 `hmmbuild --hand`, pressed them with `hmmpress`, and stored the
aligned sequences as lookup tables (generated by Rusty-HMMER
`web/anarci/generate-data.py`).

- Source: IMGT®, the international ImMunoGeneTics information system®,
  https://www.imgt.org
- License: **CC BY 4.0** (https://creativecommons.org/licenses/by/4.0/), as
  stated in IMGT's terms of use (https://www.imgt.org/about/termsofuse.php),
  in effect since 1 July 2026. Notice in [`licenses/IMGT-data.txt`](licenses/IMGT-data.txt)
- Reference: Giudicelli V, Chaume D, Lefranc M-P. IMGT/GENE-DB: a comprehensive
  database for human and mouse immunoglobulin and T cell receptor genes.
  Nucleic Acids Res. 33(suppl_1):D256–D261 (2005).
- IMGT® is a registered trademark of CNRS; it is used here only to attribute
  the data.

### Pfam 37.0 — built-in domain library

`wasm/pfam/pfam-biologics.h3m.gz`: 61 Pfam-A families (listed in
`wasm/pfam/pfam-biologics.json`) extracted from the Pfam 37.0 release
(https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/) and pressed
with HMMER 3.4 `hmmpress`.

- License: **CC0 1.0** (public domain dedication); no attribution required.
- Reference: Mistry J, et al. Pfam: The protein families database in 2021.
  Nucleic Acids Res. 49(D1):D412–D419 (2021).

Pfam families requested by accession in the dialogs are fetched at run time
from the InterPro API of EMBL-EBI (https://www.ebi.ac.uk/interpro/) and are not
bundled.

---

## 3. Test fixtures (package tests only)

- `files/tests/search-*` and `files/tests/scan-*`: Pfam 37.0 models (CC0);
  synthetic variants of Pfam family sequences (CC0); five unmodified
  UniProtKB/Swiss-Prot entries (P01009, P02144, P02647, P68871, P69905) from
  The UniProt Consortium under **CC BY 4.0**
  (https://www.uniprot.org/help/license). Expected values were produced by the
  C HMMER 3.4 programs.
- `src/tests/anarci-fixture.json`: ANARCI results for chains from Datagrok's Bio
  sample data and for synthetic chains recombined from IMGT germline genes
  (IMGT data, see above).

---

## 4. Tools used only to build the bundled data (not distributed)

- MUSCLE (public domain; Edgar RC, Nucleic Acids Res. 32(5):1792–1797, 2004),
  run by ANARCI's build pipeline to align J genes.
- HMMER 3.4 `hmmbuild` and `hmmpress` (BSD-3-Clause), which built the HMMs.
- Unmodified HMMER 3.4 and ANARCI 2024.05.21 (with Biopython) were run in Docker
  as reference oracles to verify the package's results.
