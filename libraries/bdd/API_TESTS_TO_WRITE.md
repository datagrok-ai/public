# API tests to write, and the UI scenarios they let go

The audit of 2026-10-07 applied `CLAUDE.md` "Anything an API test can check" to every bdd project. Scenarios an
existing test already covered were removed then; these are the ones that wait for a test. Write the test, run it,
then remove the scenario or lines named, with one line in the feature's description pointing at the test.

## Still to write

| Behaviour | Where the test goes | What it lets go |
|---|---|---|
| The toolbox Search box with "and" / "or" (GROK-20229) | with the fix, beside xamgle `features/search.dart` | — (the cases left the feature already); the fix is deferred |

## Written on 2026-10-18

| Behaviour | Test | What went |
|---|---|---|
| SequenceTranslator: Convert over all of `cyclized.csv`, Combine Sequences, Combine Sense+Antisense, the Translator's outputs and bulk conversion, GROK-20958 | SequenceTranslator `PolyTool: Convert` (`cyclized.csv`), `PolyTool: Combine Sequences`, `OligoRenderer: converters`, `Formats support`, `Translator: bulk`, `Translator: app` (`src/tests`) | the value lines of `oligo-polytool` (Convert's HELM, missing values and V3000; the combined row), `oligo-nucleotide-grid` (the combined duplexes of rows 1 and 34) and `oligo-toolkit-translator-structure` (the sample's translations, the HELM rows, the bulk values) |

## Written on 2026-10-12

| Behaviour | Test | What went |
|---|---|---|
| Projects with Data sync, snapshots, a calculated-column chain through a round trip and a source rename | ApiTests `Dapi: projects: data sync` (`src/dapi/projects.ts`) | `projects-data-sync` S5 and the reopen halves of S3/S4; the post-reopen value checks of PowerPack `persistence-sources`, `persistence-northwind` |
| `getRecentEntities()` after a project is saved and opened | ApiTests `Dapi: entities: recent` | nothing: `browse-my-stuff` keeps its wait, the tree reads Recent before the server lists the project (1 run in 3 without it); Recent is built from audit events, so the save puts the project there, not the open |
| A plain user does not list another user's unshared connection | datlas `test/users_groups/read_routes_authorization_test.dart` | `browse-platform-and-databases` "Another user does not see a connection…" |
| A second root space of the same name and a rename onto one are refused; a copy outlives its deleted original | datlas `test/spaces/spaces_test.dart` | the rename refusal scenario of `spaces-rename`; the survival half of `spaces-entity-ops` |
| Sticky meta keyed by molecule through a clone, a CSV round trip and another SMILES; Clear values; no values after the schema is deleted (GROK-18980) | ApiTests `Dapi: sticky meta: molecules`, datlas `test/sticky_meta` | `persistence-and-delete`'s CSV scenario and the reload of its Clear values scenario |
| Database meta: one field deleted leaves the others; string-list Values round trip; an annotation absent on its namesakes | DBTests `Database Meta: DbSchemaInfo` | `database-meta`'s namesake lines |
| Every DBTests query declaring `meta.testExpectedRows` returns that count; a choice parameter filters | DBTests `Queries: expected rows` (`src/expected-rows-test.ts`) | the `^France$` line of `queries/parameterized-queries` |
| Enrichment: apply, edit, delete, the FK lookup, a clashing name | PowerPack `src/tests/enrichment.ts` | the row-by-row checks and the replay half of `enrichment/data-enrichment` |
| A three-sheet workbook imports with its sheets, columns and values | PowerPack `src/tests/excel.ts` | the content claims of `io/xlsx-open`, `io/xlsx-shared-with-me` |
| Power search templates and syntax ("new", "use", "1+1", "Project[0-9]+") | PowerPack `src/tests/search.ts` | those rows of `search/power-search-enter` (the "Project[0-9]+" Enter stays: it is the one case with no suggestion menu) |
| Summary columns and the Forms viewer through a project; the Tags renderer (GROK-20888) | PowerGrid `src/tests` | the project halves of `grid/summary-columns`, `viewers/forms/forms-persistence`; the GROK-20888 scenario |
| Chem: cliffs at cutoff 20, elemental analysis on SMARTS, chemical space on an all-empty column, R-groups with none found, BitBIRCH / similarity matrix / Cluster MCS, reactions, notation round trips, the other fingerprints, a mixture renderer | Chem `src/tests` | the permutation scenarios of those Chem features (the cutoff-95 cliffs scenario stays: the only proof the dialog's cutoff reaches the analysis) |
| Bio: scan liabilities, similarity with a short/long reference (GROK-20963), sequence space and cliffs on HELM/MSA, a numbering round trip, multi-letter separator → HELM | Bio `src/tests` | `annotate/annotate`, `calculate/scoring` 3–5, the `analyze/other-notations` Outlines, `projects/round-trips` |
| Peptides: exports, mutation-cliff counts, MCL threshold → cluster count, LST statistics | Peptides `src/tests` | `sar/export`, `sar/mutation-cliffs`, the `sar/similarity-threshold` Outline, the LST statistics lines of `sar/tooltips` |
| Dendrogram centroid and median linkage with inversions (GROK-19595) | Dendrogram `hierarchicalClustering` tests | `clustering/chem-dialog` centroid |
| EDA control comparisons on a filtered table (GROK-20795) | EDA `Control comparisons: dialog` | `analyze/filtered-group-comparison` S2 |
| Predictive models: data checks, Ignore missing, One-hot on apply, a second model's own name | xamgle "Core: xamgle: Predictive modeling" | — (the EDA model features claim the gestures) |
| Predictive models: engine discovery, `/ml/zip`, blobs and images, save / suggested / delete | ApiTests `Dapi: models` (`src/dapi/models.ts`) | — (the EDA count checks prove the UI's Save and Delete) |
