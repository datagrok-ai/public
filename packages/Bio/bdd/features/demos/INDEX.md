# Demo clips

The bioinformatics showcase, in the order its reel cuts it. Each clip is a feature in the `bdd/features/demos/`
of the package named, and a test of what it shows. `npx grok-bdd guide features/demos/<feature>` in that
package's `bdd/` films it into `guides/<feature slug>/<scenario slug>/`; `guide-render.py --reel` cuts the
clips together. The `/bdd-answer` skill, section 5, is the recipe.

| # | Clip | Package | Feature |
|---|---|---|---|
| 1 | Every notation opens as a Macromolecule column | Bio | `01-sequence-formats.feature` |
| 2 | Rings, cross-links and chains drawn as connections | Bio | `02-biln-connections.feature` |
| 3 | HELM to BILN and back to plain text | Bio | `03-notation-conversion.feature` |
| 4 | Any HELM drawn in the grid, edited in one editor | Helm | `03b-helm-showcase-and-editor.feature` |
| 5 | Custom notation to HELM and molecules | Bio | `03c-polytool-custom-notation.feature` |
| 6 | Align, then compare every sequence with the one you click | Bio | `04-alignment-and-differences.feature` |
| 7 | Rank the peptides by their similarity to a lead | Bio | `04b-rank-by-similarity.feature` |
| 8 | A sequence's structure in 2D and 3D | Bio | `05-sequence-structure.feature` |
| 9 | Sequences become molecules, monomer by monomer | Bio | `06-to-atomic-level.feature` |
| 10 | DNA and RNA, nucleotide by nucleotide | Bio | `07-nucleotides-to-atoms.feature` |
| 11 | siRNA duplexes, every modification in view | Bio | `08-sirna.feature` |
| 12 | One siRNA from notation to duplex | UsageAnalysis | `21-oligo-toolkit.feature` |
| 13 | Antibody numbering in the browser | Bio | `09-antibody-numbering.feature` |
| 14 | CDR3 extracted with its IMGT positions | Bio | `10-numbered-regions.feature` |
| 15 | Pfam domains, germlines and profile logos | Bio | `11-hmmer-domains.feature` |
| 16 | A profile HMM built from selected chains | Bio | `11b-profile-hmm.feature` |
| 17 | From antibodies to a CDR3 SAR | Peptides | `12-cdr3-sar.feature` |
| 18 | Mutation cliffs and the invariant map | Peptides | `13-mutation-cliffs.feature` |
| 19 | Clusters, residues and cliffs in more viewers | Peptides | `14-sar-viewers.feature` |
| 20 | Find peptides by motif and by similarity | Bio | `15-search.feature` |
| 21 | Map the sequence space and find the activity cliffs | Bio | `16-sequence-space-and-cliffs.feature` |
| 22 | Cluster the peptides into a tree and cut it | Bio | `18-clustering.feature` |
| 23 | Pick out the peptides with one residue at one position | Bio | `18b-composition.feature` |
| 24 | Browse the monomer libraries | Bio | `17-monomer-libraries.feature` |
| 25 | Enumerate a peptide library | Bio | `19-polytool-enumeration.feature` |
| 26 | Protein structures in 3D | UsageAnalysis | `20-protein-structures.feature` |
