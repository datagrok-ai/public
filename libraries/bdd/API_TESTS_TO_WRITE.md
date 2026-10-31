# API tests to write, and the UI scenarios they let go

The audit of 2026-10-07 applied `CLAUDE.md` "Anything an API test can check" to every bdd project. Scenarios an
existing test already covered were removed then; these are the ones that wait for a test. Write the test, run it,
then remove the scenario or lines named, with one line in the feature's description pointing at the test.

## Still to write

| Behaviour | Where the test goes | What it lets go |
|---|---|---|
| The toolbox Search box with "and" / "or" (GROK-20229) | with the fix, beside xamgle `features/search.dart` | — (the cases left the feature already); the fix is deferred |
| MPO: Score as on rising and falling curves and inside a flat part; a log10 curve with <-1 (empty) and >-1 (whole curve), >0 with Limit value gives 0; MPO Score by Profile on a qnum column uses the profile's qualifiers; the coverage count on a log10 curve at 0 or below, after a compute without coverage (no cached desirability reused), and in Chem:mpoCalculate's result table | Chem src/tests/mpo-tests.ts ('MPO: qnum scoring', 'MPO: coverage') | — (Chem bdd/features/mpo/MISSING.md items 6–7) |

## Written on 2026-10-30

| Behaviour | Test | What went |
|---|---|---|
| A layout with a logarithmic scatter plot saves and loads after its X column is removed (public#2867) | ApiTests `Shell: ViewLayout` | PR 4414 scatter-plot regression scenario (dropped before merge) |
| A histogram on a table with no rows renders (core#998) | ApiTests `AI: Viewers: Histogram extras` | PR 4414 histogram regression scenario |
| A datetime Split with Split Map month comes back from a layout (public#2562) | ApiTests `AI: Viewers: BarChart extras` | PR 4414 bar-chart regression scenario |
| A pie chart moved to a table with the same columns keeps its Category and angle aggregation (core#982) | ApiTests `AI: Viewers: PieChart JS API` | PR 4414 pie-chart regression scenario |
| A violin over a double column holding a raw NaN draws and takes a hover (core#960) | ApiTests `AI: Viewers: BoxPlot extras` | box-plot MISSING.md entry |
| The trellis plot's inner box plot look survives the view's table being replaced (public#3692) | ApiTests `AI: Viewers: Trellis Plot` | trellis-plot MISSING.md entry |
| A grid viewer bound to a second table keeps its table and sort through a project (public#2718) | ApiTests `Dapi: projects: viewers` | PR 4414 grid regression scenario |
| Compare Tables: the rows the diff holds and its counters | ApiTests `DataFrame: Methods` `compareTables`; ddt `table_comparison_test.dart` | PR 4464 compare-tables "With no value columns…" |
| Select > Duplicates with Superfluous rows only | ddt `data_utils_test.dart` `selectDuplicates (superfluousOnly)` | PR 4464 select-menu "Superfluous rows only…" |
| A formula naming its column as a string recalculates after an edit (core#1128) | ApiTests `DataFrame: Calculated columns` (unskipped) | — |
| A child of a scaffold drawn with an explicit H opens with no warning (public#2626) | Chem `scaffold tree` | PR 4414 scaffold-tree regression scenario |
| Use as filter after a second view removed its Structure card (public#2628) | Chem `substructure filters` | PR 4414 substructure-filter regression scenario |

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
