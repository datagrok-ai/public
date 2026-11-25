# Demo clips

The bioinformatics showcase, in the order its reel cuts it: six chapters, from one column for every notation to
molecules made from sequences and Chem's tools on them, then oligonucleotides, peptide series, library design and
antibodies. Each clip is a feature in the `bdd/features/demos/` of the package named, and a test of what it shows;
its feature name is `<chapter> · <chapter title> · <topic>`, which the clip's title card shows above the scenario.
`npx grok-bdd guide features/demos/<feature>` in that package's `bdd/` films it into `guides/<feature slug>/<scenario
slug>/`; `guide-render.py --reel` cuts the clips together in this order. The `/bdd-answer` skill, section 5, is the
recipe.

| # | Clip | Package | Feature |
|---|---|---|---|
| **1** | **Any notation, one column** | | |
| 1.1 | Every notation opens as a Macromolecule column | Bio | `01-sequence-formats.feature` |
| 1.2 | Rings, cross-links and conjugates drawn as connections | Bio | `02-biln-connections.feature` |
| 1.3 | HELM to BILN and back to plain text | Bio | `03-notation-conversion.feature` |
| 1.4 | Any HELM drawn in the grid, edited in one editor | Helm | `03b-helm-showcase-and-editor.feature` |
| 1.5 | Custom notation to HELM and molecules | Bio | `03c-polytool-custom-notation.feature` |
| **2** | **Every sequence is a molecule** | | |
| 2.1 | A peptide drug in 2D and 3D | Bio | `05-sequence-structure.feature` |
| 2.2 | Sequences become molecules, monomer by monomer | Bio | `06-to-atomic-level.feature` |
| 2.3 | Peptide drugs measured like small molecules | Bio | `06b-chemical-properties.feature` |
| 2.4 | One chemical motif across every peptide | Bio | `06c-substructure-filter.feature` |
| 2.5 | Monomers are molecules too | Bio | `06d-monomers-as-molecules.feature` |
| 2.6 | RNA and DNA, nucleotide by nucleotide | Bio | `07-nucleotides-to-atoms.feature` |
| **3** | **Oligonucleotides** | | |
| 3.1 | siRNA duplexes, every modification in view | Bio | `08-sirna.feature` |
| 3.2 | One siRNA from notation to pattern | UsageAnalysis | `21-oligo-toolkit.feature` |
| **4** | **A peptide series, by sequence and by chemistry** | | |
| 4.1 | Align a peptide family, then compare it with one peptide | Bio | `04-alignment-and-differences.feature` |
| 4.2 | Rank the peptides by their similarity to a lead | Bio | `04b-rank-by-similarity.feature` |
| 4.3 | Find peptides by motif and by similarity | Bio | `15-search.feature` |
| 4.4 | Map the sequence space and find the activity cliffs | Bio | `16-sequence-space-and-cliffs.feature` |
| 4.5 | The same peptides in chemical space | Bio | `16b-chemical-space.feature` |
| 4.6 | Cluster the peptides into a tree and cut it | Bio | `18-clustering.feature` |
| 4.7 | Pick out the peptides with one residue at one position | Bio | `18b-composition.feature` |
| **5** | **Design a library** | | |
| 5.1 | Browse the monomer libraries | Bio | `17-monomer-libraries.feature` |
| 5.2 | Enumerate a peptide library | Bio | `19-polytool-enumeration.feature` |
| 5.3 | Rank the library by chemistry | Bio | `19b-library-chemistry.feature` |
| **6** | **Antibodies and proteins** | | |
| 6.1 | Antibody numbering in the browser | Bio | `09-antibody-numbering.feature` |
| 6.2 | CDR3 extracted with its IMGT positions | Bio | `10-numbered-regions.feature` |
| 6.3 | Pfam domains, germlines and profile logos | Bio | `11-hmmer-domains.feature` |
| 6.4 | A profile HMM built from selected chains | Bio | `11b-profile-hmm.feature` |
| 6.5 | From antibodies to a CDR3 SAR | Peptides | `12-cdr3-sar.feature` |
| 6.6 | Mutation cliffs and the invariant map | Peptides | `13-mutation-cliffs.feature` |
| 6.7 | Clusters, residues and cliffs in more viewers | Peptides | `14-sar-viewers.feature` |
| 6.8 | Protein structures in 3D | UsageAnalysis | `20-protein-structures.feature` |
