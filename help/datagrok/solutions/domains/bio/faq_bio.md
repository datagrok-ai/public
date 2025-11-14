---
title: "FAQ: Peptides"
sidebar_position: 1
keywords:
 - peptides
 - macromolecules
 - biologics
 - peptide SAR
---

To view FAQ on other topics, see [FAQ](../../../resources/faq.md).

## Data format and interoperability

### How does Datagrok access sequence data?

Datagrok can load sequences from CSV, XLSX, PDB, FASTA, or retrieve them from connected
databases and other [data sources](../../../../access/access.md#data-sources).

### What sequence formats does Datagrok support?

* BILN
* FASTA
* HELM
* Separator-based
* Connection-rules notations, and [other formats](../../../../access/files/supported-formats.md#bioinformatics)

Datagrok also natively supports SMILES and can convert sequences into molecular
form for further analysis.

### Does Datagrok support natural and unnatural amino acids, custom monomers, and peptide linkers?

**Yes**. Datagrok supports sequences composed of any building blocks defined in
the monomer library, including:
* Natural and non-natural amino acids
* Custom residues and linkers (represented via HELM components such as CHEM/BLOB)

Monomer definitions can be extended or modified through the [built-in monomer library management system](bio.md#manage-monomer-libraries). Bulk import with [Pistoia HELM](https://github.com/PistoiaHELM/HELMMonomerSets)
JSON. Once defined, monomers are recognized consistently across the platform for
rendering, SAR, enumeration, and other operations.

### Does Datagrok support different peptide topologies, including cyclic or modified structures?

**Yes**. Datagrok supports a wide range of peptide topologies through HELM,
BILN, and connection-rules notations, including:

* Linear, cyclic, and bicyclic peptides
* Branched and dendritic structures
* Lipidated peptides
* Custom cyclization patterns and noncanonical linkages

### Does Datagrok support peptide conjugates?

**Yes**. They are either represented in HELM completely, or separated in
different columns for analysis.

### Can Datagrok convert between molecular and sequence representations?

**Yes**. Datagrok supports sequence-to-sequence and sequence-to-molecule
conversion between all supported formats. When viewing a sequence, hovering over
a monomer highlights the corresponding fragment in the generated molecule.

Molecule-to-sequence conversion is supported for SMILES-to-HELM.

For more details, see [Format conversion](bio.md#format-conversion) 
or watch the [RDKit UGM](https://www.youtube.com/watch?v=la-kj52djeI) presentation.

### Can Datagrok integrate with registration systems and ingest associated assay results?

**Yes**. Datagrok integrates with:

* [Registration systems](../chem/chem.md#compound-registration-systems)
* [Assay and dose-response data](https://github.com/datagrok-ai/public/tree/master/packages/Curves#integration-with-assay-plates)
* [Other data sources](../../../../access/access.md#data-sources)

> _Developers_: Read about [integration options](../../teams/it/integration-story.md) and [JS development](../../../../develop/develop.md).

### Can Datagrok integrate with a proprietary monomer library and support periodic updates?

**Yes**. Datagrok can connect to proprietary monomer libraries through custom
library providers that work with file-based, database, or API-based sources.
Libraries can be accessed directly for immediate read/write access or
synchronized via ETL workflows.

Admin-controlled editing is supported through Datagrok's 
[permissions management system](../../../../govern/access-control/access-control.md#permissions).

## Visualization and UI

### Can Datagrok show sequences, molecules, properties, and assay data in a single table?

**Yes**. Tabular views are powered by the [grid viewer](../../../../visualize/viewers/grid.md), 
which displays sequences, molecular structures, calculated properties, and assay data side-by-side in a
single interactive table. You can add, rearrange, or 
[compute columns](../../../../transform/add-new-column.md) as needed.

Peptide properties and [descriptors](../chem/descriptors.md) can be calculated
from the **Top Menu** and/or dedicated [info panes](../chem/info-panels/info-panels.md). 
Custom calculations can be defined
through [user-defined functions](../../../../compute/scripting/scripting.mdx)
(JavaScript, Python, Julia, R, or MATLAB) or by integrating external services.

To learn about visualization and analytics capabilities for peptides, see
[Bioinformatics](bio.md). For chemical utilities, see
[Cheminformatics](../chem/chem.md). For assay data integration options, see
[Curves](https://github.com/datagrok-ai/public/tree/master/packages/Curves#integration-with-assay-plates).

## Sequence analysis

### Can I filter or search peptides by architecture or sequence similarity?

* Text-based filtering: available for FASTA using [text filter](../../../../visualize/viewers/filters.md#text-filter)
* Sequence similarity and diversity searches: use **Bio > Search > Similarity** and **Diversity** tools
* Substructure search: convert sequences to molecular form and apply the [substructure filter](../chem/chem.md#structure-search)

### Does Datagrok support sequence alignment?

**Yes**. Datagrok supports MSA using [K-Align](https://github.com/TimoLassmann/kalign) 
for canonical sequences and [PepSEA](https://github.com/Merck/PepSeA) 
for non-canonical sequences ([learn more](bio.md#multiple-sequence-alignment-msa)). 
For visual summaries of aligned sequences, use the [WebLogo](../../../../visualize/viewers/web-logo.md) viewer.

### Does Datagrok support peptide SAR analysis?

**Yes**. Datagrok provides an interactive environment for exploring
sequence-activity relationships. It combines visualization, clustering, and
statistical tools to help you identify key mutation sites, visualize sequence
variability, and compare peptide properties with assay results. [Learn more](peptides-sar.md).

### Does Datagrok support mutation-cliff analysis for peptides?

**Yes**. The [Mutation Cliffs viewer](peptides-sar.md#sequence-variability-map)
identifies sequence positions where single-monomer substitutions cause
significant activity changes. You can configure the minimum activity difference
threshold to focus on the most impactful mutations and filter by specific
positions.

### Does Datagrok support matched pairs analysis for peptides?

**Yes**. The **Mutation Cliffs** mode of the 
[Sequence Variability Map](peptides-sar.md#sequence-variability-map) systematically compares peptides
that differ by a single monomer substitution. For each position-monomer pair, it
displays the number of matched pairs (circle size), mean activity change
(color), and detailed statistics in the **Context Panel**.

### Does Datagrok support invariance mapping and positional filtering for peptide sequences?

**Yes**. Datagrok includes tools for visualizing and filtering sequence variability across positions:

* **Invariance mapping**: use [WebLogo](../../../../visualize/viewers/web-logo.md) to show monomer frequencies at each position, or [Sequence Variability Map](peptides-sar.md#sequence-variability-map) in **Invariant Map** mode to show position-specific frequencies and related statistics
* **Positional filtering**: click any cell in the **Sequence Variability Map** to
  select all sequences containing that monomer-position pair. Selected sequences
  are displayed in the **Selection** table and **Context Panel**. You can also use the monomer
  search filter to find specific residues across positions

### Does Datagrok support statistical analysis for peptides?

**Yes**. The sequence [Position Statistics](../../../../visualize/viewers/position-statistics-viewer.md) viewer displays box or violin plots of
selected properties across different motifs or positions, helping correlate
sequence patterns with assay results or molecular properties. You can also use
hypothesis testing tools, including
[MVA](../../../../explore/multivariate-analysis.md) and
[ANOVA](../../../../explore/anova.md).

<!-- add a wiki page and a link for the Sequence Position Statistics viewer-->