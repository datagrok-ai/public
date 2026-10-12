# API tests to write, and the UI scenarios they let go

The audit of 2026-10-07 applied `CLAUDE.md` "Anything an API test can check" to every bdd project. Scenarios an
existing test already covered were removed then; these are the ones that wait for a test. Write the test, run it,
then remove the scenario or lines named, with one line in the feature's description pointing at the test.

| Behaviour | Where the test goes | What it lets go |
|---|---|---|
| A project saved with Data sync reopens with `.data-sync = success`; a snapshot reopens with no mark; a calculated-column chain survives a data-sync round trip and a source rename | ApiTests `src/dapi/projects.ts` (new) | UsageAnalysis `projects/projects-data-sync` S5 and the reopen halves of S3/S4; PowerPack `add-new-column/persistence-sources`, `persistence-northwind` |
| `getRecentEntities()` holds a project after `project.open()` | ApiTests `src/dapi/entities.ts` | the Recent claim's wait in `browse/browse-my-stuff` |
| A plain user does not list a connection another user owns and has not shared | datlas `test/users_groups/read_routes_authorization_test.dart` | `browse/browse-platform-and-databases` "Another user does not see a connection…" (two sign-ins) |
| A second root space of the same name, and a rename onto one, are refused; a copy outlives its deleted original | datlas `test/spaces/spaces_test.dart` | one of the two refusal scenarios (spaces-create / spaces-rename); the copy half of spaces-entity-ops |
| Sticky meta keyed by molecule survives a clone, a CSV round trip and a non-canonical SMILES; Clear values; no values after the schema is deleted (GROK-18980) | Chem `src/tests` (sticky meta) and datlas `test/sticky_meta` | `sticky-meta/persistence-and-delete` 3.2–3.4 |
| Database meta: one field deleted leaves the others; string-list Values round trip; a table's or column's annotation is absent on its namesakes | DBTests `src/db-annotations/db-annotations.ts` | `sticky-meta/database-meta` isolation lines |
| Every DBTests query that declares `meta.testExpectedRows` returns that count; a parameter choice filters | DBTests (new `expected-rows-test.ts`) | the `^France$` line of `queries/parameterized-queries` |
| SequenceTranslator: Convert over all of `cyclized.csv`, Combine Sequences, Combine Sense+Antisense, the Translator's outputs and bulk conversion, GROK-20958 | SequenceTranslator `src/tests` | value lines of `oligo-polytool`, `oligo-nucleotide-grid`, `oligo-toolkit-translator-structure` |
| Enrichment: apply, edit, delete, the FK lookup, a clashing name | PowerPack `src/tests/enrichment.ts` (new) | row-by-row checks and the replay half of `enrichment/data-enrichment` |
| A three-sheet workbook imports with its sheets, columns and values | PowerPack `src/tests/excel.ts` | the content claims of `io/xlsx-open`, `io/xlsx-shared-with-me` |
| Power search templates and syntax ("new", "use", "1+1", "Project[0-9]+") | PowerPack `src/tests/search.ts` | those rows of `search/power-search-enter` |
| Summary columns and the Forms viewer through a layout and a project; the Tags renderer (GROK-20888) | PowerGrid `src/tests` | the project halves of `grid/summary-columns`, `viewers/forms/forms-persistence`; the GROK-20888 scenario |
| Chem: activity cliffs at cutoffs 20/95, elemental analysis on SMARTS, chemical space on an all-empty column, R-groups "No R-Groups were found", BitBIRCH / similarity matrix / Cluster MCS, reactions, notation round trips (GROK-20956), the other fingerprints, a mixture renderer | Chem `src/tests` | the permutation scenarios of those Chem features |
| Bio: scan liabilities, similarity with a short/long reference (GROK-20963), sequence space and cliffs on HELM/MSA, a numbering round trip, multi-letter separator → HELM | Bio `src/tests` | `annotate/annotate`, `calculate/scoring` 3–5, the `analyze/other-notations` Outlines, `projects/round-trips` S2 |
| Peptides: exports, mutation-cliff counts on the 200-row subset, MCL threshold → cluster count, LST statistics | Peptides `src/tests` | `sar/export`, `sar/mutation-cliffs`, the `sar/similarity-threshold` Outline |
| Dendrogram centroid and median linkage with inversions (GROK-19595) | Dendrogram `src/tests/hierarchical-clustering-tests.ts` | `clustering/chem-dialog` centroid |
| EDA control comparisons on a filtered table (GROK-20795) | EDA `src/tests` | `analyze/filtered-group-comparison` S2 |
| The toolbox Search box with "and" / "or" (GROK-20229) | with the fix, beside xamgle `features/search.dart` | — (the cases left the feature already) |
