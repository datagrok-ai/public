# Mol Track changelog

## v.next

* Added MolTrack EMS: compounds, batches, additions, assays and results on a domain schema, each compound identified from its structure before it is saved
* Compounds: Search by structure (substructure, superstructure, similarity, exact) and by any column in a filter builder with nested conditions
* Compounds: Bulk registration previews every row (new, already registered, repeated, unreadable) and writes the Grok IDs back to the table
* Home page shows how many compounds, batches, assays, runs, results and additions are registered, each count opening its list
* Grok IDs open their compound from any table, from the context panel and from `?corporate_compound_id=` links

## 0.0.1 (2025-03-08)