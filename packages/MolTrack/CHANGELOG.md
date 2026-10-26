# Mol Track changelog

## v.next

* Compounds: Batches are a grid you type new batches into; the Grok Batch ID opens the batch
* Assay results: Showed concentrations and CLint values as formatted numbers instead of raw floats
* Compounds, batches and assay results show and check the same in child tabs as in their own lists, and in the `/domains` view
* Added MolTrack EMS: compounds, batches, additions, assays and results on a domain schema, each compound identified from its structure before it is saved
* Compounds: Search by structure (substructure, superstructure, similarity, exact) and by any column in a filter builder with nested conditions
* Compounds: Bulk registration previews every row (new, already registered, repeated, unreadable) and writes the Grok IDs back to the table
* Grok IDs open their compound from any table, from the context panel and from `?corporate_compound_id=` links
* Compounds: Every saved compound shows its stored identity (parent structure, formula, InChIKey, matching rule), not only a new one
* Batches: A batch page shows its compound (structure, Grok ID and name) at the top
* Batches: A Grok Batch ID (DGB-) shows its batch and compound in the context panel; the legacy "Databases | MolTrack" panel is gone
* Compounds: ⋯ > Open as table opens the listed compounds, structure filters included, in a table view
* Added `getEmsCompound`, `getEmsBatch`, `searchEms` and `registerEmsCompounds` for other packages and AI assistants
* Batches and assay results: Lines name the compound and line up in columns; cut text shows in full on hover
* Compounds: Card titles line up, and bulk registration previews draw larger structures

## 0.0.1 (2025-03-08)