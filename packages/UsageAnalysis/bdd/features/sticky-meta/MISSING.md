# Sticky Meta features — what is missing

The four features in this folder, against TestTrack `StickyMeta/` and the old Playwright specs. Triaged with
the lead on 2026-10-11; what is not planned is under "Won't do", with the reason.

## Resolved (2026-10-12)

- **Server fixtures for an entity type and a schema.** `a Sticky Meta schema {string} for molecules, with the entity
  type {string}, is on the server with the properties:` makes both through the API (removed now and at feature end,
  read back), and `the Sticky Meta schema {string} is deleted on the server` deletes the schema alone. The
  Backgrounds of `add-and-edit` and `persistence-and-delete` no longer walk the New Entity Type and New Schema
  dialogs, and their UI delete scenarios are gone; `persistence-and-delete` deletes the schema through the API before
  its GROK-18980 scenario, keeping the type as the ticket describes. `schema-and-type` keeps the UI walk and claims
  each OK and DELETE on the server (`the entity type … should be on the server`, `the Sticky Meta schema … should be
  on the server with the properties "rating: int, …"`, `no Sticky Meta schema named …`, `no entity type named …`),
  the property types included, which the Edit schema dialog shows in no readable control.
- **The values written onto the shared SPGI molecules.** Deleting a schema and its entity type removes the values
  and the molecules' entity records with them: after many runs the stand holds no property value without its
  property and no record without its type.
- **Database meta put back.** `database-meta` clears the Database meta of PostgresTest (NorthwindTest) through the
  API when it starts and when it ends.
- **`the context panel should show {string}` and a database column.** The step compares the names case aside, so
  `database-meta` claims the column as the current object (`the context panel should show "categoryid"`).
- **A known failure pinned to its line.** `@known-failure:<line>` passes only on a failure of the step at that
  feature line (library `harness.ts`, `compile.ts`).
- **Other schemas that match molecules.** The cell's dialog and pane put each schema's fields in a section of its
  own, named after the schema (`div-section--<schema>`, `df_properties.dart` `buildPropertiesValues`), and the
  features name their fields in their own section (`Rating input in "bdd-sm-keep-{time}" section in context panel`).
- **The intermittent "pane empty after a reload"** (2026-10-12, traced with logging in the client): after a reload
  the click on "cell 2 of Structure" landed on row 1, whose values an earlier scenario had cleared. On a fresh page
  the grid asked Chem for its molecule renderer while Chem was still loading, reported nothing pending, and grew the
  rows once the renderer arrived — under the point already aimed at. The grid now counts a package renderer still
  loading as pending work (`grid_core.dart` `detectRenderers`), `user opens … dataset` waits for the grid to settle,
  and the feature claims the current row after every cell click. 10 runs in a row passed since.
- **A reload with a schema view in front logged NullError.** `DbSchemaView` and `DbInfoView` restored from their
  address no longer read the connection before it resolves (`db_views.dart`), and a view that never loaded stays
  out of Recent.

## Still open

Nothing planned.

## Edge cases (`schema-edge-cases.feature`) — still open

Written for the edge-case feature, with no library, binding or code change. It claims the confirmation before a
taken schema name drops the existing schema's properties (CANCEL keeps them), a matching expression not in
"tag=value" form refused, an address naming no schema, and core#677 (a schema deleted from the Browse tree stays in
the tree; known failure). What follows has no source for its expected behaviour.

### Questions for owners

1. **OK in the confirmation for a taken schema name.** Does it replace the existing schema (its properties and their
   values dropped, as the message warns) or add a second schema under the same name? The dialog builds a new schema
   with a new id (`property_schema_meta.dart` `createNewCommand`), so the message and the save may disagree.
2. **Text typed into an `int` or `double` property in the cell's Sticky meta dialog.** Is it refused in the input,
   or saved as empty?
3. **Edit for all properties with no row selected.** Is the command offered, and what does Save write?
4. **The molecule column removed or renamed while its sticky columns are shown.** What the sticky columns and the
   column's Sticky meta pane should show afterwards.
5. **An entity type deleted while a schema is still associated with it** (the reverse order of core#316): refused,
   or the schema detached and its values removed?

## Won't do

- **The marker of a cell with metadata.** A status provider on the grid that reports every drawn dot is the kind
  of per-render signal the library does not take; the features claim the values through the marker's tooltip,
  built from the same cache the dot is drawn from.
- **The tooltip right after a fresh load.** It needs an `aria-busy` on the grid while `fetchForCells` runs; the
  features read the cell's Sticky meta pane first, whose fields appear once the values are read.
- **A fresh session of the running account.** `user signs in as themselves again` replays the session the feature
  started with, so TestTrack 3.3's "relogin" is claimed as a reload; a step minting a new session for the running
  account is not worth it for one case.
- **A balloon after the last check.** A short "quiet" wait for the error floor before the last scenario ends would
  catch a late balloon; not worth a wait in every feature.

### Notes on defects met

- **Emptied fields in the cell's dialog keep their values — by design.** Save does not delete an emptied field
  (`buildPropertiesValues` deletes one only for Database meta); GROK-15602 was resolved by Clear values for the
  schema, and `persistence-and-delete` claims both.
- **The values of a deleted schema — GROK-18980.** On 2026-10-01 and since, the cell's pane shows no section and no
  value once the schema is deleted, its type still there; `persistence-and-delete` claims that.

### Not translated

- A server restart (the primary copy-clone-delete case): a feature cannot cause one.
- The connection-level Database meta of CHEMBL (the primary database-meta case): a Postgres connection lists
  catalogs, and the platform builds no Database meta pane for such a connection.
