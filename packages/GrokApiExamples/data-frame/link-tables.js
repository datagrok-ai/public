// The following example shows tables link feature.
//
// "linkTables" links two tables based on the key columns. The last parameter specifies what
// should be synchronized (current record, filter, selection, and mouse-over record),
// see constants definitions below.
//
// https://datagrok.ai/help/features/tables-link

gr.loadDataFrame("/demo/demog.csv").then((demog) =>
    gr.loadDataFrame("/demo/demog-types.csv").then(function (demogTypes) {
        gr.addTableView(demog);
        gr.addTableView(demogTypes);
        gr.linkTables(demogTypes, demog,
            ['sex', 'race'], ['sex', 'race'],
            [CURRENT_ROW_TO_FILTER, MOUSE_OVER_ROW_TO_SELECTION]);
    }));

// Link types:
//
//   CURRENT_ROW_TO_ROW
//   CURRENT_ROW_TO_SELECTION
//   CURRENT_ROW_TO_FILTER
//   MOUSE_OVER_ROW_TO_SELECTION
//   MOUSE_OVER_ROW_TO_FILTER
//   FILTER_TO_FILTER
//   FILTER_TO_SELECTION
//   SELECTION_TO_FILTER
//   SELECTION_TO_SELECTION
