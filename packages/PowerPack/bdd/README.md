# PowerPack behavioral tests

Gherkin features under `features/`, compiled by `@datagrok-libraries/bdd` (`public/libraries/bdd`)
into the Playwright specs under `generated/` — committed, never edited by hand. They translate the
TestTrack PowerPack section (`packages/UsageAnalysis/files/TestTrack/PowerPack/`), one feature per
dialog session or data source, each saying in its description what it claims and what it leaves out.

| Folder            | Features | What it claims |
|-------------------|----------|----------------|
| `add-new-column/` | add new column, formula editor, functions panel, formula refreshing, persistence over a file in My files, persistence over Northwind query results, result type and validation, edit column formula, script in formula, append tables | the Add New Column dialog from the toolbar and from Edit: tooltips, resize, a formula built from autocomplete and columns dragged from the grid's header and from the dialog's column list, the input history; autocomplete on a letter, Ctrl+Space and "$", the signature on hover, column references highlighted in --blue-2 (pasted, autocompleted, dropped, and the GROK-17004 formula); functions inserted by the plus icon and by a drag, the picked column passed when its type fits, the list reordered by the picked column's type and kept by "By name"; a chain of calculated columns recalculated from the Formula pane of the context panel; calculated columns following a rename and an edit of their source and coming back from a project saved with Data sync (GROK-17109); the result type (auto, picked by hand, a mismatch that disables OK, plain text as a template, string and plain text over a number (public#4295), an empty formula adding a plain empty column), a formula naming a missing column, a vector function inside arithmetic, a name already taken; Edit Column Formula from the Formula pane renaming and recalculating a column and the one over it, a grid edit recalculating a running total, the preview of a dialog-made column under its own name (public#4293); a user's JavaScript script in a formula; Append Tables keeping a formula, its OK disabled until a table is checked (core#784). `MISSING.md` holds the open questions |
| `annotation-regions/` | annotation regions, drawing, interaction, interaction on every viewer, titles, persistence, the Tools menu, the other 2D viewers | regions drawn and edited on the scatter plot, line chart and the other two-dimensional viewers, hovered and clicked, their titles and title strips, kept with the view — through the regions' hit areas and readings |
| `formula-lines/` | the Formula Lines dialog, the axis Annotations group, regressions, GitHub regressions | lines and bands added, edited, hidden and removed in the dialog and the look it writes; the Annotations group on a numeric axis; the fixed regressions, and #2487 (a datetime axis), #671 (the preview follows the viewer's axes and the selected line), #2747 (one dataframe line on a scatter plot and a line chart) |
| `home/` | home widgets, home as sharing user | the Home widgets, hiding them, Customize, a reload, the Spotlight tabs, the Community links, Usage and Reports; a second account without Usage and Reports, an unread notification on the server |
| `search/` | power search enter | Enter on the md's queries once the results have answered (`aria-busy=false`), with and without the suggestion menu, and the suggestions walked with the arrows |
| `navigation/` | direct link loading | a project opened by its direct link on a fresh page load and the same project from the Browse tree |
| `io/` | xlsx open, xlsx shared with me | a three-sheet workbook through My files, a drop, the file chooser, File > Open and Shared with me — exactly three table views |
| `enrichment/` | data enrichment | creating, applying, editing and deleting enrichments, several at once, a layout and a reopened project, the enrichments of a referenced column offered through a foreign key |
| `identifiers/` | identifiers, identifiers of a database with catalogs — both `@skip`, development only | Configure Identifiers... on a connection: the dialog that asks for the schema, the editor with an identifier previewed from its table, the saved configuration under the connection, opened again and removed; on a database with catalogs, the catalog asked for first, its schemas alone, and after a reload the identifier's row read from that catalog and not the connection's own |

`annotation-regions/` and `formula-lines/` came from `packages/UsageAnalysis/bdd/features/viewers/`:
the feature belongs to PowerPack (the Formula Lines dialog is `PowerPack:formulaLinesDialog`).

`add-new-column/persistence-northwind` is `@full-stand`: it needs the Postgres Northwind database
(the NorthwindTest connection of dev); a stand without it runs with `--grep-invert @full-stand`.

**For developers only.** Everything about identifiers and catalogs below — the two `identifiers/`
features, the catalog database `catalogs-db.mjs` makes, and the development-only package tests that use
it (PowerPack's catalog enrichment tests, ApiTests' `getCatalogs`) — is for a developer whose change
touches Configure Identifiers, the db-explorer library, catalogs or enrichments of a database with
catalogs. QA, and anyone else running the suites, can skip it: the features and tests are skipped and
need nothing set up. Do not make the fixture unless you are asked to or are working on that area.

Both features of `identifiers/` are `@skip`: they are for development and do not run in CI, as the
package tests of the same area (`src/tests/db-explorer.ts`, category DB Explorer, each with a
`skipReason`) do not. They write to the one file every identifiers configuration of a stand is kept
in. To work on Configure Identifiers or the db-explorer library, take `@skip` off the feature and the
`DEV_ONLY` reason off the tests, `grok-bdd compile`, run them on a local stand, and put both back
before committing.

`identifiers/identifiers-catalogs` also needs a database with catalogs behind a connection named
`CatalogsDemo`; unskipped, it still skips itself on a stand where the connection is missing or its
database does not answer. `features/identifiers/catalogs-db.mjs` makes the fixture on any machine with
Docker and Node (it needs nothing installed), from the two files beside it: `catalogs-fixture.sql`, the
one place the catalogs are written, and `catalogs-connection.json`, the connection with a placeholder for
the password. Every step is safe to repeat; every docker and sqlcmd call has a time limit. From the package's
`bdd` directory:

```bash
node features/identifiers/catalogs-db.mjs up                    # the container, the login, the catalogs
node features/identifiers/catalogs-db.mjs connection            # saves CatalogsDemo on the grok default server, tests it
node features/identifiers/catalogs-db.mjs connection --host <alias>   # another stand of ~/.grok/config.yaml (never dev)
node features/identifiers/catalogs-db.mjs status                # what is there and what is missing
node features/identifiers/catalogs-db.mjs up --reload           # makes the three working catalogs again
node features/identifiers/catalogs-db.mjs down                  # stops the container; nothing is deleted
node features/identifiers/catalogs-db.mjs down --purge --yes    # deletes the container and its volume
```

The database is one SQL Server container, which other projects may share. What the script keeps to:

- The container `dg-catalogs-db`, labelled `ai.datagrok.fixture=catalogs-db`, from
  `mcr.microsoft.com/mssql/server:2022-CU27-GDR1-ubuntu-22.04` (pinned), Developer edition,
  `--restart unless-stopped`, `--memory 8g` with `MSSQL_MEMORY_LIMIT_MB=6144`. Its databases live in the
  named volume `dg-catalogs-db-data` (`/var/opt/mssql`). SQL Server answers on `127.0.0.1:14330` only.
- The sa password is the container's `MSSQL_SA_PASSWORD`, generated by whoever creates the container. A
  script reads it with `docker inspect` and uses it inside the container only, never printed or on a
  command line:
  `docker exec -i dg-catalogs-db sh -c '/opt/mssql-tools18/bin/sqlcmd -S localhost -U sa -P "$MSSQL_SA_PASSWORD" -C -b -i /dev/stdin' < file.sql`.
  SQL Server keeps it in the volume, so a container made again on the old volume needs the same one:
  the script caches it, and without it in the cache it stops and says so.
- A project that keeps databases there has a login of its own that sees only them:
  `CREATE LOGIN ... CHECK_POLICY = OFF`, `DENY VIEW ANY DATABASE`, and `ALTER AUTHORIZATION` of each of its
  databases to the login; SQL Server then lists for the login its own databases, `master` and `tempdb`.
  This fixture's login is `catalogs_demo`, the owner of `research`, `archive`, `sales` and
  `tenant_001`..`tenant_200`; the script and `catalogs-fixture.sql` drop or change no other database.
  Another project may add its own databases to the same server under its own login.
- The script generates the login's password and keeps it in `features/identifiers/.catalogs-db/`
  (git-ignored): `secrets.json` holds it and the cached sa password, `catalogs-connection.json` the
  connection with it, which `connection` saves. Without the file, `up` sets a new password, and the
  connection has to be saved again.
- `down` stops the container: every database on it is down, the other projects' too, until `up`.
  `down --purge` deletes the container and its volume, and with them every database on the server, not
  only this fixture's: it says so and does nothing without `--yes`.

The feature counts 205 catalogs: the three it works on and the two hundred, which the login owns, and
`master` and `tempdb`. What else the server holds is not counted.

`host.docker.internal` in the connection is the Docker host as Grok Connect sees it from its container;
on Docker Desktop it reaches the port bound to 127.0.0.1. On Linux a container reaches no port bound to
the host's loopback, and `host.docker.internal` resolves only where the stand maps it to `host-gateway`:
attach the database to the network of the stand's Grok Connect (`docker network connect <network>
dg-catalogs-db`) and save the connection by the container's name,
`node features/identifiers/catalogs-db.mjs connection --server dg-catalogs-db --port 1433` (the address is
kept for later runs). After the database restarts, Grok Connect may keep dead pooled connections to it
for up to 30 minutes: restart the Grok Connect container rather than wait.

What the stand needs: PowerPack published from the same checkout, Chem (the SPGI features pick
Molecule functions), the demo files under `System:DemoFiles/` and a "My files" share for the
current user; `home/` and `io/xlsx-shared-with-me` also the library's sharing account, `home/` an
administrator account of its own (created once per stand), `enrichment/` the `System:Datagrok`
connection. `home/` signs in as those other accounts, and they see released packages only, not the
debug version the running account published: PowerPack (the widget hosts' names) and UsageAnalysis
(the Usage widget) go out with `grok publish <stand> --release`.

Run from the package's `bdd` directory. PowerPack is a package of the `public/` pnpm workspace, so
it resolves the library with nothing to link:

```bash
cd public && grok setup                               # once per checkout: the pnpm workspace
cd libraries/bdd && npm run build                     # the library; dist/ is not committed
cd ../../packages/PowerPack/bdd
npx grok-bdd run generated --workers=2 --reporter=line --grep-invert @full-stand
npx grok-bdd run generated/add-new-column             # one folder
```

The bindings are the package's own screen parts and checks: `bindings/add-new-column.ts` — the
Add New Column dialog's editor, hint, column list, functions list and preview; the column list's rows
found by the grid's own readings; the functions list's order and what its top rows take; the
highlighted column references and their computed color. `home.ts` — the Home widgets, their stored
settings, Spotlight's tip of the day and the administrator account the widget feature signs in as;
`search.ts` — the Home search and its suggestions; `enrichment.ts` — the Enrich pane, its editor and
the enrichment files on the server; `identifiers.ts` — the identifiers editor, its preview, and the
configurations the server keeps (one file, `System:AppData/PowerPack/db-explorer/configs.json`). What other suites use too is library vocabulary: signing in as
another account, a reload, files put into the user's files or a space, a dropped file, a direct link,
a named table's cells, element sizes, and the Formula Lines dialog — PowerPack's, but the NX chain of
UsageAnalysis drives it as well.
