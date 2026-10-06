# Sticky Meta features — what is missing

The four features in this folder were written with the vocabulary the library has today (no binding
or core change). They pass on a local stand built from master. This file lists what they could not
say, what each missing piece must check, and where its signal would come from. Every item is a
gap of the bdd library or of the core, not of the features.

## 1. Server fixtures for an entity type and a schema (cleanup rule)

Every feature that needs a schema makes the entity type and the schema through the New Entity Type
and New Schema dialogs, and deletes both through the Schemas and Types galleries in its last
scenarios. The sweep the hard rule asks for is in the library now (`the Sticky Meta schema {string}
and entity type {string} are removed now and at feature end`, which also takes a `{time}` family's
members over an hour old); what is still wanted is making them and reading them back:

```gherkin
Given no Sticky Meta schema named "bdd-sm-cells-{time}" and no entity type named "bdd-sm-cells-type-{time}" are on the server
Given a Sticky Meta schema "bdd-sm-keep-{time}" for molecules is on the server with the properties:
  | rating | int    |
  | notes  | string |
Then the entity type "bdd-sm-type-{time}" should be on the server
Then the Sticky Meta schema "bdd-sm-schema-{time}" should be on the server with the properties "rating: int, notes: string, verified: bool, review_date: datetime"
Then no Sticky Meta schema named "bdd-sm-schema-{time}" should be on the server
```

- `grok.dapi.stickyMeta.deleteSchema(id)` takes an id the JS `Schema` has no getter for: the
  library's sweep reads it off the Dart entity (`grok_Entity_Get_Id(s.dart)`). The Tutorials steps
  `the entity type {string} should exist` and `the Sticky Meta schema {string} should exist` are the
  natural ones to promote next, with the negative and the property list.
- With the Given above, the Backgrounds of `add-and-edit` and `persistence-and-delete` shrink to it
  and their two delete scenarios go; `schema-and-type` keeps its UI walk, adds the server claims
  after each OK and DELETE, and claims the property types (the Edit schema dialog shows names in
  inputs, but the types not in a readable control — the old spec read them through the API).

- **The values written onto the shared SPGI molecules.** Sticky meta belongs to the molecule, not to
  the schema's lifetime: `add-and-edit` leaves rows 1-3 of spgi-100 with values, and
  `persistence-and-delete` leaves row 2's (4, good) — deleting the schema hides them (the feature
  claims that), but whether the server drops them is not read (GROK-18980 says they can outlive the
  schema). Wanted: the fixture Given above also deletes the values of the molecules a feature wrote,
  and a Then that reads a molecule's values from the server (`getAllValues`) to prove it.

## 2. Database meta put back at the feature's end

`database-meta` writes the Comment and LLM Comment of the `public` schema of NorthwindTest, the
Comment, LLM Comment and Row Count of `public.categories` and every field of its `categoryid` column, and its last scenario empties them through the pane. A
killed run leaves the values on a table every NorthwindTest feature reads. Wanted:

```gherkin
Given the Database meta of the "categories" table of "NorthwindTest" is put back at feature end
Given the Database meta of the "categoryid" column of the "categories" table of "NorthwindTest" is put back at feature end
Given the Database meta of the "public" schema of "NorthwindTest" is put back at feature end
```

Reading the entity's properties before the feature and writing them back (or deleting what was not
there) at its end: the entity is the `EntityRecord` the pane builds from
`<connection id>:<catalog>:<schema>:<table>[:<column>]` (`df_properties.dart`,
`_initContextPanelListener`, the `DatabaseProperties.isDatabaseEntity` branch).

## 3. The marker of a cell with metadata

The grid paints the dark-blue dot of a cell that has sticky meta and the blue circle of a sticky
column's header on its canvas (`df_properties.dart` `drawEllipsis`, `StickyMetaColumnRenderer`), and
reports neither. The features claim the values through the tooltip shown over the cell's top right
corner, which is built from the same cache the dot is drawn for. A reading would let a feature claim
the dot itself, and its absence:

- a status provider on the grid (`grid.statusProviders`, the JS map `DG.Widget.addStatusProvider`
  writes to) registered by `StickyMeta._initStickyRenderer`, reporting a `sticky meta marker of
  cell <row> of <column>` hit area (the 10×13 px corner `checkMouseOnIndication` hit-tests) for every
  drawn dot, a `sticky meta cells of <column>` reading (the rows with a dot) and `sticky column
  <column>` for a header with the circle.

## 4. The tooltip right after a fresh load

The marker tooltip is built from the cache at the moment of the hover (`renderMetaTooltip`), and the
cache is filled 300 ms after a draw (`debounce(grid.onAfterDrawContent)` → `fetchForCells`). Right
after a reload a hover reads an empty cache and shows "No sticky meta for this cell" for a cell
that has metadata — a negative claim there would pass on a not-yet-loaded cell. The features read the
cell's Sticky meta pane first (its fields appear only once the values are read). A signal would let
the tooltip be claimed on its own: an `aria-busy` on the grid while `fetchForCells` runs, or the
reading of item 3 appearing only once the visible rows are fetched.

## 5. `the context panel should show {string}` and a database column

The step compares the current object's `friendlyName` (`ColumnInfo "Categoryid"`, capitalized) with
the panel's text (`categoryid`, as the database has it), so it cannot pass for a column of a database
table. `database-meta` claims `context panel should contain text "categoryid"` instead, which does
not prove the column is the current object. Either the step compares case-insensitively, or it
reads the panel's header for the name.

## 6. Smaller things the review found

- **A fresh session of the running account.** `user signs in as themselves again` puts the session the
  feature started with back on the page (`signInAsSelf` replays the remembered token), so the
  "relogin" of TestTrack 3.3 is a page load under the same session. A step that mints a new session
  for the running account (as `user signs in as {string}` does for another login) would make it one.
- **A known failure pinned to its line** (a library wish met while this round still had one): a
  `@known-failure` scenario passes on a failure at any of its steps; the harness could take the line
  expected to fail (`@known-failure:<line>`) and report any other.
- **A balloon after the last check.** `no error or warning balloon should have been shown` reads the
  log once; a balloon a DELETE answers with after it is cleared at the next scenario's start and never
  read in the last one. A short "quiet" wait for the floor (the shell's pending requests) would close it.

## 7. Defects met

### Emptied fields in the cell dialog keep their values — by design

Save in the cell's Sticky meta dialog does not delete a field that was emptied: `buildPropertiesValues`
(`df_properties.dart`) deletes an emptied property only for Database meta. GROK-15602 ("cannot remove
values") was resolved by adding Clear values for the schema, so `persistence-and-delete` claims both:
an emptied field keeps its value, and Clear values removes them.

### The values of a deleted schema — GROK-18980

GROK-18980 (Open) reports the values of a deleted schema still shown on the molecule until the entity
type is deleted. On 2026-10-01 the cell's Sticky meta pane showed no section and no value once the
schema was deleted, the type still there; `persistence-and-delete` claims that. If the ticket's
screenshot shows another place, that place needs a claim too.

### A reload with a schema view in front logs NullError

Clicking a schema node (Databases > Postgres > NorthwindTest > Schemas > public) opens the schema's
view; reloading the page with that view in front logs `NullError: method not found:
'PackageEntityMixin_id' on null` from `DbSchemaView.saveStateMap` (db_views.dart:192), called by the
docking tabs and the history as the browse tree is opened again — the view restored from its address
(`DbSchemaView.fromPath`) has no connection until `loadStateMap` finds it. Seen twice in a row in
`database-meta`; that feature closes the views before each reload so its claims stay about the
Database meta. Not reproduced yet by a plain reload and a click outside the suite. The same
`NullError: method not found: 'PackageEntityMixin_id' on null` was reported on a project save as
GROK-15609 (Won't fix, no stack).

## 8. Other schemas that match molecules on the stand

The cell's Sticky meta dialog and pane show one section per schema matching the molecule
(`stickyMetaEditorForCell`), and a section is not a container: its header and its inputs are siblings
in one form. `Rating input in "Sticky meta" dialog`, `Save button in "Sticky meta" dialog` and
`"Add rating as a column" button` therefore name the fields of every matching schema. On a stand where
only the feature's schema matches molecules (a local or CI stand) they are unique; on dev on
2026-10-01 four more schemas matched — "Molecule meta", "Highlight", "TestSchema1" (the TestTrack
fixture with its own rating/notes/verified/review_date) and "schema for tutorial" (left by a hand walk
of the tutorial) — and `add-and-edit` and `persistence-and-delete` failed there on five Save buttons
and on TestSchema1's Rating, while `schema-and-type` and `database-meta` passed. Wanted, either:

- a container per schema in the editor, named after the schema (as `div-section--<name>`), so a
  feature scopes its phrases: `Rating input in "bdd-sm-cells-{time}" section in "Sticky meta" dialog`;
- or a Given that names the schemas matching molecules other than the feature's own and skips the
  feature with the reason (a stand-state gate), when cleaning the stand is not the feature's to do.

## Not translated

- A server restart (the primary copy-clone-delete case): a feature cannot cause one.
- The connection-level Database meta of CHEMBL (the primary database-meta case): a Postgres
  connection lists catalogs, and the platform builds no Database meta pane for such a connection.
