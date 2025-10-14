// Master-detail table linking feature.
//
// "linkTables" links two tables based on the key columns. The last parameter specifies what
// should be synchronized (current record, filter, selection, and mouse-over record),
// see constants definitions below.
//
// Once you run an example, it would make sense to open both tables side-by-side.
// Here, demog-types is the master table, and demog is the detail table. Specifying
// [CURRENT_ROW_TO_FILTER, MOUSE_OVER_ROW_TO_SELECTION] means that current record in the master
// table controls the filter in the second table. So, for instance, you can use UP and DOWN
// keys in the demog-types grid, and you would only see corresponding records in the other table.
//
// https://datagrok.ai/help/explore/link-tables

let demog = await grok.data.loadTable("https://public.datagrok.ai/demo/demog.csv");
let demogTypes = await grok.data.loadTable("https://public.datagrok.ai/demo/demog-types.csv");

grok.shell.addTableView(demog);
grok.shell.addTableView(demogTypes);
grok.data.linkTables(
  demogTypes,
  demog,
  ["sex", "race"],
  ["sex", "race"],
  [DG.SYNC_TYPE.CURRENT_ROW_TO_FILTER, DG.SYNC_TYPE.MOUSE_OVER_ROW_TO_SELECTION]
);

// Link types:
//   From SYNC_TYPE enum:
//     CURRENT_ROW_TO_ROW
//     CURRENT_ROW_TO_SELECTION
//     CURRENT_ROW_TO_FILTER
//     MOUSE_OVER_ROW_TO_SELECTION
//     MOUSE_OVER_ROW_TO_FILTER
//     FILTER_TO_FILTER
//     FILTER_TO_SELECTION
//     SELECTION_TO_FILTER
//     SELECTION_TO_SELECTION
