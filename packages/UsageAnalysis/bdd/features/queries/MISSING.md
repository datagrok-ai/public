# Queries — what the features leave out

## Still open

- core#229 — an emptied int parameter keeps OK enabled while an emptied string disables it; awaits a ruling (the
  help says every parameter is required unless nullable). Claim when ruled: ApiTests `Widgets: InputForm validation
  and events`.
- core#456 — closing the query editor with unsaved SQL asks nothing; claim when done: the confirmation on close
  (`query-lifecycle`).

## Won't do

- The run dialog's checks (min and max, text in an int, a validator, choices, a nullable left empty): ApiTests
  `Widgets: InputForm validation and events`, `Widgets: InputForm skipLogic`, `Property: Header parsing`.
- Rename... refusing a name with no letter: `DbEntityMeta.getRenameDialog`'s rule, the same for every entity.
- A post-process that throws: the query call's error, as for any function.
- CANCEL of the run dialog runs nothing: `parameterized-queries` claims it closes the dialog, `scripts-run` that a
  CANCEL then one run logs one call.
