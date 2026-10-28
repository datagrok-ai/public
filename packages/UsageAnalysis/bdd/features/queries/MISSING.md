# Queries edge cases — what the features cannot say yet

What the edge-case features of this folder (`query-parameters-edge-cases.feature`,
`query-editor-edge-cases.feature`) leave out, written for them: no library, binding or code change was made.
Each item names the case, why it is not a scenario yet, and what would make it one. Reproductions are on
localhost (core a6c9fcf85b), 2026-10-26.

## Still open

### Waiting for core#456

1. **Closing the query editor with unsaved SQL**: the confirmation is requested in core#456.

### Decided and written

A post-process script that throws shows an error balloon carrying the error and opens no table (decided
2026-10-28: post-processing is a step of the data transformation, so the SQL's raw result is not shown):
`query-editor-edge-cases`.
An empty value is valid by default, for a string as for an integer (core#229): `query-parameters-edge-cases`
claims both, the string half a known failure.

### API or server tests to write

2. **Two queries with one name (core#546)** — already requested in `projects/MISSING.md` §2 (ApiTests
   `dapi/queries.ts`).
3. **A nullable parameter left empty is passed as null** (`help/datagrok/concepts/functions/func-params-annotation.md`,
   "Initial values and optional parameters"): the rows the query returns for a null are the server's answer —
   a DBTests or ApiTests case. The feature claims only that the dialog lets it run.

## Won't do

- **A new query logs "FormatException: Unexpected end of JSON input" (core#603).** `query-lifecycle.feature` and
  `new-sql-query.feature` already run a new query in the editor and claim no errors logged; both pass on
  localhost, so the ticket does not reproduce here.
- **Visual query edge cases.** A Where value that matches no row is claimed in `visual-query-parameter.feature`
  (preview 0 rows).
- **CANCEL of the run dialog runs nothing** — `parameterized-queries.feature` claims CANCEL closes the dialog,
  and the script's `scripts-run.feature` claims a CANCEL followed by one run logs one call.
