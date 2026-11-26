# Review: 01-getting-started (Datagrok basics)

Fresh-eyes review (2026-11-22, ~20:50-21:30), then the lead's 21:00 round ("no jumps", Browse sections opened to
show examples, global search clips 12+), which the reviewer took over for this section. Videos live in
`~/Desktop/grok/Core/howto-videos/01-getting-started/`.

| NN | Clip | Verdict |
|---|---|---|
| 01 | Find your way around Datagrok (3:34) | refilmed (reworked): every Browse section is opened to show what is in it (My stuff's buckets, Spaces, the apps, the Dashboards gallery searched for "cliffs", the Biologics connection with its queries, schemas and tables, Platform and its Plugins view, Files > Demo), and every move after the start is filmed: the molecules and sequences tables come from Browse (Files > Demo > chem > smiles.csv, Files > Demo > bio > peptides.csv) instead of appearing off camera. Remaining on screen: the unfiltered Dashboards gallery shows "cicd_tests" and "Auto Test- ..." cards for about 3 s; hovering My stuff pops the account card; smiles.csv shows a Bio menu next to Chem (product). |
| 02 | Open your data (0:40) | OK. The "drag a file onto the window" stop shows no drag (a browser drop can't be filmed); the table just appears under the line. |
| 03 | Work with the grid (0:46) | OK. |
| 04 | Format and color-code columns (0:36) | refilmed: "Two digits after the point" showed while "max two digits after comma" turned 73.20 into 73.2 and 93.00 into 93; now "At most two digits after the point". |
| 05 | The context panel (0:53) | refilmed (jump fixed, see the coordinator's 21:20 note): the table of molecules was opened off camera mid-clip; now it is opened from Browse > Files > Demo > chem > smiles.csv on camera. |
| 06 | Show and hide panels (0:28) | OK. |
| 07 | Get help (0:24) | OK. |
| 08 | Search from Home (0:26) | refilmed: the "cars" search listed some 18 test files (cars.d42, cars.kmz, cars_2_cols.d42 ...) and opened a 5-row mtcars file; it now searches "titanic" and opens the Demo titanic.csv (1,310 rows). Title made specific ("Search functions, help and files ..."), since 12 covers dashboards and queries. |
| 09 | Work with several tables side by side (0:21) | OK. After docking iris to the right, the left half shows Home until cars' tab is clicked (product); the clip names that click. |
| 10 | Find a column in a wide table (0:24) | OK. |
| 11 | Learn with interactive tutorials (0:30) | OK. |
| 12 | Search everything from one box (0:20) | new: "Threshold" finds dashboards, files and queries at once; pointing at a dashboard shows its preview, a double-click opens it. |
| 13 | Ask a question in the search box (0:28) | new: typing "biologics" suggests the sentences the queries answer; "biologics assays for E. coli" runs the query and shows its table; "adcs linked to GROKMOL-000001" lists the ADCs made with one drug. |
| 14 | Look up an ID (0:31) | new: GROKADC-000001 shows its record as a card; three GROKMOL- IDs show three cards; Open table puts them in a table whose IDs stay links (a click shows the record in the context panel). |

Verdicts: 7 OK, 4 refilmed (01, 04, 05, 08), 3 new (12, 13, 14). Nothing left that needs work in this section.

## Global search: what works on this stand

- Biologics (installed, its database filled) gives all three: dashboards named "ADCs by IC50 Threshold" and others,
  query patterns ("biologics assays for ${organism}", "adcs linked to ${drugID}", "biologics ADCs with IC50 ...",
  "biologics antibody profile for ... targeting ..."), and identifiers GROKADC-, GROKMOL-, GROKSEQ-, GROKORG-.
- Typing any text that is part of a query pattern opens the suggestion list over the results; clip 12 uses
  "Threshold" (in no pattern) so the list stays closed.
- Chembl's patterns need the ChEMBL database: not used. The identifier cards are titled with the semantic type's
  internal name ("1 DG_BIOLOGICS_DRUG_ID object"), and the opened table's column is named DG_BIOLOGICS_DRUG_ID:
  cosmetic, product side.
- The Projects list shows "ADCs by IC50 Threshold" twice (two copies on the stand).

## Data claims (checked 2026-11-22 ~22:40, one frame of the delivered video per line that states a fact)

All true, nothing changed. 01 (tour): My stuff opens on Recent, Favorites, Shared with me and the My ... buckets;
Spaces, Apps (Demo, Tutorials, Chem, Bio, ...) and the Dashboards gallery show real entries; Databases lists
ClickHouse, MS SQL, MariaDB, MySQL, Neo4j, Oracle, Postgres, Snowflake ...; the Biologics connection shows its
queries and Schemas; Files > Demo opens on its sample folders; demog opens with 11 columns and 5,850 rows; each
menu holds what its line names (Edit: Find and Replace, Add Rows, Add New Column; View: Layout, Edit Tooltip,
Presentation Mode; Select: All, Random Rows, Duplicates, Missing Values; Data: Aggregate Rows, Join, Link, Split
Column; ML: Cluster, Reduce Dimensionality, Impute); smiles.csv brings Chem (and Bio, the known product bug),
peptides.csv brings Bio only. Data's Join and Link Tables are greyed with one table open; they are there, so the
line stands. 12: "Threshold" finds two dashboards (the heading reads "Projects": the product's own word), two
files and two queries, all with the word in the name; the dashboard opens (559 rows). 13: "biologics" suggests the
four sentences; "biologics assays for E. coli" returns E. coli rows; "adcs linked to GROKMOL-000001" lists ADCs
with that drug's card under them. 14: GROKADC-000001 shows ADC_1's card; three GROKMOL- IDs give three cards; Open
table gives 3 rows; a click on GROKMOL-000002 shows Drug_2 in the context panel. 02: cars opens with 30 rows; DONE
gives the pasted 4 rows and 2 columns. 03: Ctrl+End lands on row 5,850; the header sorts AGE 89 down, then 18 up.
04: "At most two digits": 73.20 shows as 73.2, 93.00 as 93; Linear and Categorical color as named. 05: Stats lists
count, unique, min, max ...; Properties shows MW 349.43 and the rest for the clicked molecule. 08: results under
Scripts, Functions, Demos, Help, Samples; titanic.csv opens with 1,310 rows. 10: the status bar reads Columns: 88.
11: the tracks run from Exploratory Data Analysis to Cheminformatics and Bioinformatics. 06, 07, 09 state nothing
about data.

## Titles and descriptions

All plain, no banned words. 08's title and description were updated for "titanic"; 01's description and the plan
table describe the new tour; blocks for 12-14 added.

## Systematic

- Off-camera `Given` steps mid-scenario are the jumps the lead saw; section 01 now has none (checked with a scan of
  every feature for a `Given` after the first `When`/`Then`).
- A tour this long passes Playwright's 2-minute test timeout: film it with `--timeout=420000` (film.sh passes extra
  arguments through to the test runner).
- Home search and Browse show whatever test files and projects the stand has collected (18 cars.* files, Auto Test
  projects, duplicate dashboards); every clip that lists them depends on a clean stand.

## Merges after the lead's 23:00 round ("17 seconds is too short", 2026-11-22 ~23:20)

| NN | Clip | Absorbed | Checked |
|---|---|---|---|
| 03 | Work with the grid (1:28) | 10 Find a column in a wide table | spgi-100 is open from the start (its tab beside demog), switched to by a filmed tab click; pin, Columns: 88, search "herg", the grid goes to hERG Classification with Structure pinned |
| 06 | Arrange your workspace (0:53) | 09 Work with several tables side by side | cars and iris docked side by side, then the toggles; Presentation mode shows only the current view (cars), so the line says "just the current view" (re-rendered from the take) |
| 07 | Get help and learn (1:16) | 11 Learn with interactive tutorials | Context Help closed with its toggle before the Help menu; the tutorial opens its own "demog (2)" next to the clip's demog; closing it brings Browse and the context panel back (as in the old 11) |
| 08 | Search everything from Home (1:09) | 08 Search from Home, 12, 13, 14 | "scatter plot" (help link tooltip, no click), "Threshold" (dashboard preview, double-click opens it, the house tab back to Home), "biologics" with the suggestion list held by a hover, "biologics assays for E. coli" (E. coli rows), GROKADC-000001 card, three GROKMOL- cards, Open table (3 rows), a click on GROKMOL-000002 shows Drug_2. The "titanic" file search and the second sentence ("adcs linked to ...") were dropped as repetition |

All four: film.sh passed, audit 0 presses off the lit element, no off-camera Given after the first visible step,
frames checked at every line and result.
