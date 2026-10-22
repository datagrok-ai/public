# Mol Track changelog

## v.next

* MolTrack EMS: Added a second app over the `moltrack_ems` domain schema, beside the container-backed one: compounds, batches, additions, assays, runs and results as lists and entity pages, with registration numbers shown as DG-000123 and DGB-000123
* MolTrack EMS: Registering a compound computes its parent structure, registration hashes, InChI and formula in the browser (Crux); a structure already registered under the compound matching rule (a package setting) is named under Structure as soon as it is drawn, with the registered compound shown, and the hashes are folded under Registration details
* MolTrack EMS: Every list has a Filters pane: conditions row by row, and on a structure column has substructure, is substructure of, is similar to (with a threshold) and is exactly
* MolTrack EMS: A new batch is added to a registered compound picked by its Grok ID or name; the batch page shows that compound, which cannot be changed once the batch is saved
* MolTrack EMS: Added bulk registration of compounds (Compounds › ⋯ › Import…): a mapping of every column, a preview of what each row becomes, and a per-row report with each row's Grok ID — a compound registered already is reported, not written — and a link to the compounds registered; a file is written whole or not at all by default
* MolTrack EMS: Batches are imported through Batches › ⋯ › Import…, a column of Grok IDs naming each batch's compound; a Grok ID no compound has fails on its own row
* MolTrack EMS: A Grok ID anywhere in the platform shows the registered compound in the context panel, without opening the app first
* MolTrack EMS: Lists open newest first, a compound clicked in a list shows in the context panel, and the tree reuses the open view

## 0.0.1 (2025-03-08)