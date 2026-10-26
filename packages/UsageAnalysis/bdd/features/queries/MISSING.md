# Queries edge cases — what the features cannot say yet

What the edge-case features of this folder (`query-parameters-edge-cases.feature`,
`query-editor-edge-cases.feature`) leave out, written for them: no library, binding or code change was made.
Each item names the case, why it is not a scenario yet, and what would make it one. Reproductions are on
localhost (core a6c9fcf85b), 2026-10-26.

## Still open

### Steps requested

1. **A query whose internal name has a dash breaks the dashboard built on it.** Confirmed on localhost: a query
   saved with the name `probe-dash-x` (not only the friendly name) is run, its result saved as a project with
   Data sync on, and the project reopened: the "Data loading error" dialog shows
   `Sub(Sub(Admin:probe, …), x())` and `Unable to get project asset "probe"`, and the table is not shown. The
   same query named `probe_under_x` reopens with fresh data. No ticket yet (`@ticket-pending`). The editor and
   the Rename dialog turn a typed `a-b` into `AB`, and the library's `a query … on … is on the server` step
   saves through `DataConnection.query`, which does the same, so no existing step makes this name.
   Step requested (platform bindings, next to `queryTextOnServer`):
   `Given a query with the internal name {string} on {string} reads {string}` — saves `conn.query(n, sql)`
   with `q.name = n` set before `grok.dapi.queries.save`, cleaned like the other query fixtures. The scenario
   it enables (on the System:Datagrok connection; the known-failure line is the dialog-absent claim):
   ```gherkin
   @known-failure:<line> @ticket-pending
   Scenario: A dashboard of a query whose internal name has a dash reopens with its data
     Given no project named "BDDDashProj{time}" is on the server
     And a query with the internal name "bdd-dash-{time}" on "System:Datagrok" reads "select name from entity_types where name = 'Project'"
     And Databases tree node inside browse tree is expanded
     And Databases---Postgres tree node inside browse tree is expanded
     And Databases---Postgres---Datagrok tree node inside browse tree is expanded
     When user double-clicks Databases---Postgres---Datagrok---bdd-dash-{time} tree node inside browse tree
     Then the "text of cell 1 of name" reading of grid should be "Project"
     When user clicks on Save ribbon item
     And user enters "BDDDashProj{time}" into Name text input in "Save project" dialog
     And user clicks on OK button in "Save project" dialog
     And user clicks on CANCEL button in "Share BDDDashProj{time}" dialog
     And user closes all views
     And user clicks on Dashboards tree node inside browse tree
     And user enters "BDDDashProj{time}" into gallery search
     And user double-clicks on BDDDashProj{time} gallery card
     Then "Data loading error" dialog should be absent
     And the "text of cell 1 of name" reading of grid should be "Project"
   ```

### Questions for owners

2. **The dash in a query name: refuse, normalize or quote?** Should the server refuse or normalize a query
   `name` with a dash at save, or should the creation script quote it? The help caution
   (`help/access/databases/databases.md`, under Query Editor) says a dash typed in the editor breaks the
   dashboard; on this stand the editor already turns `assay-results` into `AssayResults`, so the caution
   describes only names set through the JS API or a package's `--name:` header.
3. **An empty integer and an empty string parameter (core#229, open).** The ticket asks for one rule for both;
   which one (both refused unless nullable, as the help's "Initial values and optional parameters" reads)?
4. **Closing the query editor with unsaved SQL.** `data_query_view.dart` has no close guard; only the
   page-level "Warn on unsaved changes" exists. Should closing the view ask?
5. **A Transformations step on a column the SQL no longer returns.** For a dashboard, `help/datagrok/concepts/project/dashboard.md`
   says a "Data loading error" names the failed step; what should the query itself show when run?
6. **A post-process script that throws.** Is the result shown, and how is the error surfaced?

### API or server tests to write

7. **Two queries with one name (core#546)** — already requested in `projects/MISSING.md` §2 (ApiTests
   `dapi/queries.ts`).
8. **A nullable parameter left empty is passed as null** (`help/datagrok/concepts/functions/func-params-annotation.md`,
   "Initial values and optional parameters"): the rows the query returns for a null are the server's answer —
   a DBTests or ApiTests case. The feature claims only that the dialog lets it run.

## Won't do

- **A new query logs "FormatException: Unexpected end of JSON input" (core#603).** `query-lifecycle.feature` and
  `new-sql-query.feature` already run a new query in the editor and claim no errors logged; both pass on
  localhost, so the ticket does not reproduce here.
- **Visual query edge cases.** The ones found (a removed column under a saved step, a throwing post-process)
  are owner questions above; a Where value that matches no row is claimed in `visual-query-parameter.feature`
  (preview 0 rows).
- **CANCEL of the run dialog runs nothing** — `parameterized-queries.feature` claims CANCEL closes the dialog,
  and the script's `scripts-run.feature` claims a CANCEL followed by one run logs one call.
