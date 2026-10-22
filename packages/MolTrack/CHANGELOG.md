# Mol Track changelog

## v.next

* MolTrack EMS: Added a second app over the `moltrack_ems` domain schema, beside the container-backed one: compounds, batches, additions, assays, runs and results, with registration numbers shown as DG-000123 and DGB-000123
* MolTrack EMS: A compound's parent structure, hashes, InChI and formula are computed in the browser as it is drawn; a structure registered already is named under Structure and cannot be saved, and the registration identity is folded under Registration details
* MolTrack EMS: Lists filter by structure (substructure, superstructure, similarity with a threshold, exact match) and by any other column
* MolTrack EMS: A batch is added to a registered compound picked by its Grok ID or name, and cannot be moved to another compound once saved
* MolTrack EMS: Added bulk registration (Compounds › ⋯ › Import…): a column mapping, a preview and a per-row report with Grok IDs; by default a file is registered whole or not at all
* MolTrack EMS: Batches are imported with a column of Grok IDs naming their compounds
* MolTrack EMS: The context panel shows the compound of a selected compound or batch, of a batch being registered, and of a Grok ID anywhere in the platform

## 0.0.1 (2025-03-08)