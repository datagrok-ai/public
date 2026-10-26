# Mol Track changelog

## v.next

* Compounds: Batches are a grid you type new batches into; the Grok Batch ID opens the batch
* Compounds: A saved compound shows its stored identity (parent, formula, InChIKey, matching rule)
* Compounds: ⋯ > Open as table opens the listed compounds, structure filters included
* Compounds: Larger structures on cards and in the bulk registration preview
* Batches: A Grok Batch ID shows its batch and compound in the context panel, replacing "Databases | MolTrack"
* Batches and assay results: Lines name the compound and line up in columns; numbers are formatted
* Compound, batch and assay result rules also apply in child tabs and the `/domains` view
* Added `getEmsCompound`, `getEmsBatch`, `searchEms` and `registerEmsCompounds` for other packages and AI assistants
* Added MolTrack EMS: compounds, batches, additions, assays and results on a domain schema, each compound identified from its structure before it is saved
* Compounds: Search by structure (substructure, superstructure, similarity, exact) and by any column in a filter builder with nested conditions
* Compounds: Bulk registration previews every row (new, already registered, repeated, unreadable) and writes the Grok IDs back to the table
* Grok IDs open their compound from any table, from the context panel and from `?corporate_compound_id=` links

## 0.0.1 (2025-03-08)