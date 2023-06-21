# Peptides

Provides advanced tools for analyzing collections of peptides.

Status: Beta
Owner: Leonid

## Detection and usage

Once a dataframe is opened with Datagrok, the column(s) containing peptides are detected automatically. Each amino acid is rendered with a different color, which helps to visually identify patterns. Clicking on the column header brings up the options to launch peptide-specific tools, such as [SAR](#sar) and [Peptide space](#peptides-space).

## Sar

Wiki: _Structure-Activity Relationship (SAR) is an approach designed to find relationships between chemical structure (or structural-related properties) and biological activity (or target property) of studied compounds._

Peptides SAR is an interactive analysis tool to help identify point mutations that cause large change in activity. Given a `peptide` column with aligned sequences (max length N, number of unique amino acids A) and `activity` column, it produces an N*A dataframe where (n, a) element contains a measure of the change of activity caused by mutating the amino acid at position n to a.

[Issue tracker](https://github.com/datagrok-ai/public/issues/93)

## Peptides space

Visualizes a collection of peptides in 2-dimensional space, using [t-SNE](https://en.wikipedia.org/wiki/T-distributed_stochastic_neighbor_embedding) with the [edit distance](https://en.wikipedia.org/wiki/Edit_distance) function.

[Issue tracker](https://github.com/datagrok-ai/public/issues/93)
