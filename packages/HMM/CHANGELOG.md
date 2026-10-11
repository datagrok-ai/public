# HMM changelog

## 0.1.0 (unreleased)

* ANARCI numbering engine (IMGT, Kabat, Chothia, Martin, AHo, Wolfguy) for Bio's
  Apply Numbering Scheme dialog, with chain, species and germline assignment.
* Germlines and Species (ANARCI) for antibody and TCR columns.
* Profile HMM Search over a sequence column (hmmsearch, split across workers).
* Find Domains (hmmscan) with a built-in Pfam library, Pfam accessions or HMM files.
* HMMER 3.4 engine in WebAssembly (SIMD), verified against the C programs.
* Numbering options in Bio's dialog (species, bit score threshold, domain of
  multi-domain chains); separator-notation columns.
* Pfam clan competition in Find Domains, Profile HMM Search and the Pfam panel;
  InterPro families cached in the browser.
* Align to profile (hmmalign) as an MSA engine; Homology Search (phmmer);
  Build Profile HMM (hmmbuild) with search and logo; profile HMM logos drawn
  by Bio's WebLogo (its Profile property; hmmlogo heights from this package).
* Profile HMM Logo: choose the WebLogo's sequence column and the profile (Pfam
  family, accession, ANARCI germline or a profile built in the table); a column
  aligned to a profile proposes it.
* Domain annotations are written with Bio's shared annotation helpers; the
  ANARCI panel takes its IMGT CDRs from Bio's scheme tables.
