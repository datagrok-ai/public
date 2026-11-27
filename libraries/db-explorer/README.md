# db-explorer

db-explorer is a [library](https://datagrok.ai/help/develop/develop#libraries) for the [Datagrok](https://datagrok.ai) platform.
This library will contain code for ease of exploring databases schemas and its content

## Links

The Links section of a card lists the rows of other tables that point at the card's row, or at a row the card points
at. Each group is a box headed by what its rows share with the card (`This drug`, `Same target · GROKTGT-000010`), and
each table is named, without configuration, by:

1. the configuration's `entities[table].label` (and `.short`, the code shown after it);
2. the semantic type of the identifier registered for the table (`DG_BIOLOGICS_DRUG_ID`: `Drug`), with the fixed part
   of its pattern as the code (`GROKMOL-`); a type that names the identifier rather than the records (`CHEMBL_ID`) is
   skipped;
3. the table's name without the words most tables of its schema share (`stg_`, `_v1`), codes in capitals, and without
   the name of the row it is listed under (`sequence_regions` under a sequence is `Regions`).

A table that points at the card twice, or at its own table, adds what is left of the column's name (`Hierarchy · as
parent` for ChEMBL's `parent_molregno`, `Child compounds` for a `parent_id` of the table itself), or the reference's
`reverseLabel`. The table and column behind a line, and the table's comment, are in its tooltip.

```ts
entities: {
  compounds: {label: 'Compound', short: 'CPD'},
  results_all: {more: true},   // under More: it repeats the other result tables
},
explicitReferences: [{table: 'compounds', column: 'parent_id', refTable: 'compounds', refColumn: 'id',
  label: 'Parent compound', reverseLabel: 'Child compounds'}],
```

A table's rows are read when its pane opens, and counted on its header then. Each row is one line: its key
(`uniqueColumns`, or its identifier's column), a link to its own card, and its name (`headerNames`, else a column named
like a name, else `name`, `type`); the rest of the row opens under it. When the configuration names card columns
(`customSelectedColumns`) for five tables or more, a table it says nothing about (no identifier, key, card or header
column) is taken for a view and folds under More.
