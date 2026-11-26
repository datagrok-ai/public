# Review: 12-databases (Databases and queries)

Fresh-eyes review (2026-11-22, ~22:38-22:50) of the 7 delivered clips: feature, descriptions.md block, and contact
sheets at the key lines and results of each delivered take (every delivered video matched its latest take's length
before the review). Judged on: is the one idea clear to someone new to Datagrok, is each line true when it shows,
plain wording, jumps. The m:ss in descriptions.md match the delivered videos.

| NN | Clip | Verdict |
|---|---|---|
| 01 | Browse a database | OK. Type groups, the connection's hover card (server, database), queries + Schemas, the biologics schema's tables, the drugs preview ("50 out of 200"), its Columns in the context panel, the columns in the tree, Get All → 200 rows. After Get All the left pane turns into the Toolbox by itself (the app's doing, on camera). |
| 02 | Run a saved query | OK. The Organism and Target lists do come from the database (`choices: query(...)` in Biologics/queries/antibody_profile.sql). E. coli + HER2 gives 589 rows with the saved layout (filters, grid, form, pie chart, pivot table tab); PD-1 + REFRESH gives 712 rows in the same view, pie redrawn, so "Same view, new rows from the database" is true. |
| 03 | Write a SQL query | OK. The result (6 categories, 4+3+1+9+2+4 = 23 assay types) shows under the editor; SAVE, then Browse shows the query under Biologics with its hover card. The SQL appears in one go at the end of the "Type the SQL" step (the step replaces the editor's text). |
| 04 | Build a query without SQL | OK, description fixed: it listed Pivot among the builder's rows from the start, but Pivot (and Having) only show once something is aggregated; "the chain is pivoted into columns" is now "Pivot adds a column for each chain (HC and LC)". On screen: Group by → 3 rows (CHO, E.coli, HEK293), the picker's aggregation is avg, Pivot → HC and LC columns, as the lines say. |
| 05 | Add a parameter to a query | refilmed (0:52). In the delivered take both `the video says` lines ran ahead of the text they explain: "An --input line declares the parameter and its default" sat over an empty editor and faded as the line appeared, and "The SQL uses it as @category" faded before the select line with @category was typed. The dialog's default (Developability) was on screen for under a second before Legacy was typed over it. Now each line is a `# caption:` on a hover over the typed text (`user hovers over the text "--input"` / `"@category" in code editor`), and a new hover on the Category input carries "It starts with the default, Developability" before "Type another category". Description says the dialog starts with the default. 9 Legacy assay types, as checked. |
| 06 | Save viewers with a query | OK. Layout tab, run, box plot + histogram from the toolbox, SAVE (the tab renames to "Expression yields" only here), the editor closed on camera, Browse comes back by itself, Run → 1,200 rows with both viewers. The mid-scenario `Given the tooltip is hidden` (the closed tab's tooltip stuck on the page) changes nothing the viewer needs to follow. |
| 07 | See a whole schema at once | OK. 15 table boxes (adc ... targets), sequences and heavy_chain each show in the context panel (sections collapsed: the panel shows the name and its section headers), the drugs box's menu matches the tree's (Get All, Get Top 100, New SQL/Visual Query...), Get Top 100 → 100 rows. |

Verdicts: 5 OK (01, 02, 03, 06, 07), 1 description fixed (04), 1 refilmed (05).

## Data claims (every on-screen line and description sentence that states a result, checked against the frame)

All true: 01 preview of the first rows (50 of 200) and Get All 200; 02 choices from the database (choices queries),
589 → 712 rows after PD-1; 03 six categories; 04 3 rows, avg, HC/LC; 05 Developability default, 9 Legacy rows; 06 1,200
rows with box plot and histogram; 07 15 tables, 100 rows.

## Titles and descriptions

Plain and specific, no hype words; each names where things are (Browse > Databases > PostgresDart > Biologics) and the
demo database. Fixed: 04 (Pivot only after an aggregation), 05 (the default comes first). Checked against the lead's
23:00 rules: no card, line or description mentions how a video was made, and the shortest clip is 0:43.

## Systematic

- `user replaces the code of code editor with ...` and `user appends ... to code editor` put the text in at the press,
  at the end of the step: a `the video says` line placed before them shows over the editor as it was. A line about
  typed text belongs on a following hover (`user hovers over the text "<part>" in code editor` works in CodeMirror).
- A dialog that opens with defaults and is typed into right away shows the default for under a second; a hover on the
  input with a caption gives it time.
- Clips that run a query switch the left pane from Browse to the Toolbox by themselves (01 Get All, 02, 03, 05, 06,
  07): the app's doing, on camera, not a jump.

## For the lead

- Security check (from the maker, not verified here): the Biologics connection's hover card shows server localhost,
  database `datagrok`, port 5432 (01, 02, 03, 06), and its Schemas list `public`, `information_schema` and
  `pg_catalog` next to `biologics` (01, 04, 07). If the connection's login can read the platform's own tables in
  `public`, a viewer of these clips learns where to look.
- Product (the maker's list, seen in the takes): a new query's view tab keeps "New query" after SAVE (03; 06 renames
  it because it ran from the Layout tab); the "Save query on server" tooltip and a closed tab's tooltip stay on the
  page after the pointer leaves (03, 06); the connection's hover card rises late over its own children while they
  expand (01, 02, 03).
