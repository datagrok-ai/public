# Mol Track changelog

## v.next

* Import: Added finding linked records by their fields: map a file column to Batch › EPA Batch ID, Assay Run › Cell Lot and the like
* Compounds: Added a table view with Output columns and Aggregations over batches and assay results under the filters
* Compounds: Fixed the Identity pane's form staying alive after going back from a compound or a new registration to the list
* App: Fixed reloading a view, opening a view link and Back from a record opening Compounds: the views live under `/apps/Moltrack`, where the platform serves the app
* Added Compound properties and Batch properties: add your own fields to compounds and batches, which then show on forms, grids and filters
* App: The header is u2's `AppHeader` (moved out of `@datagrok-libraries/utils`): the same icon, title, bullets, Learn more and rule
* Compounds: Added a batches grid to type new batches into
* Compounds: A saved compound shows its stored identity
* Compounds: Added ⋯ > Open as table
* Compounds: Larger structures on cards and in the registration preview
* Batches and assay results: Lines name the compound and line up in columns
* Compound, batch and assay result rules also apply in child tabs and `/domains`
* Rebuilt MolTrack without the Docker service: compounds, batches, additions, assays and results, each compound identified from its structure before it is saved
* Compounds: Search by structure (substructure, superstructure, similarity, exact) and by any column in a filter builder with nested conditions
* Compounds: Bulk registration previews every row (new, already registered, repeated, unreadable) and writes the Grok IDs back to the table
* Home page shows how many compounds, batches, assays, runs, results and additions are registered, each count opening its list
* Grok IDs open their compound from any table and from `?corporate_compound_id=` links

## 0.0.1 (2025-03-08)