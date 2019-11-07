// The following example shows tables compare feature.
//
// https://datagrok.ai/help/dialogs/compare-tables

let t1 = gr.testData('demog', 10000);
t1.name = 'demog1';
t1.rows.removeAt(4, 4, false);
t1.set('height', 1, 150.0);
gr.addTableView(t1);

let t2 = gr.testData('demog', 10000);
t2.name = 'demog2';
t2.rows.removeAt(10, 3, false);
t2.set('weight', 2, 80.5);
gr.addTableView(t2);

let values = t1.cols.names().filter(name => name !== 'subj');
gr.compareTables(t1, t2, ['subj'], ['subj'], values, values);
