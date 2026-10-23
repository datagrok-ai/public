# Mol Track changelog

## v.next

* Assay results: Showed concentrations and CLint values as formatted numbers instead of raw floats
* Compounds, batches and assay results show and check the same in child tabs and the `/domains` view as in their own lists
* Added MolTrack EMS: compounds, batches, additions, assays and results on a domain schema, each compound identified from its structure before it is saved
* Compounds: Search by structure (substructure, superstructure, similarity, exact) and by any column in a filter builder with nested conditions
* Compounds: Bulk registration previews every row (new, already registered, repeated, unreadable) and writes the Grok IDs back to the table
* Grok IDs open their compound from any table, from the context panel and from `?corporate_compound_id=` links

## 0.0.1 (2025-03-08)