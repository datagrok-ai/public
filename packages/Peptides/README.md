# Peptides

Provides sequence-activity relationship analysis for peptides.

## Detection and usage

Once a dataframe is opened with Datagrok, the column(s) containing peptides are detected automatically. 
Each amino acid is rendered with a different color based on used monomer library, which helps to visually identify patterns. 
Clicking on the column header brings up the options to launch peptide-specific tools, such as [SAR](#sar) 
and others.

## Sar

Wiki: _Structure-Activity Relationship (SAR) is an approach designed to find relationships between chemical structure
(or structural-related properties) and biological activity (or target property) of studied compounds._

Peptides SAR is an interactive analysis tool to help identify point mutations and residues that cause major change in
activity. Analysis can be triggered from context menu of a peptide column or from top menu `Bio | Analyze | SAR`.
You can choose the peptides column, activity column, configure scaling, clustering and other options.
You can read more about the tool in our [help section](https://datagrok.ai/help/datagrok/solutions/domains/bio/peptides-sar).

[Datagrok community: Macromolecules](https://community.datagrok.ai/t/macromolecules-updates/661/1)