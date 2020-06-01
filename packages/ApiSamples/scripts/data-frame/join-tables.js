// The following example joins two tables, using the specified key and value columns.
//
// https://datagrok.ai/help/dialogs/join-tables

let t1 = DG.DataFrame.fromCsv(
`id1,id2,v1,v2
a1,b4,0,0
a2,b3,1,1
a3,b2,2,2
a4,b1,3,3`);
t1.name = 't1';

let t2 = DG.DataFrame.fromCsv(
`id3,id4,v3,v4
a1,b4,5,5
a2,b3,6,6
a3,b2,7,7
a4,b1,8,8
a5,b0,9,9
a6,b5,10,10`);
t2.name = 't2';

let tj = grok.data.joinTables(t1, t2, ['id1', 'id2'], ['id3', 'id4'], ['v1', 'v2'], ['v3', 'v4'], DG.JOIN_TYPE.INNER, false);

grok.shell.addTableView(t1);
grok.shell.addTableView(t2);
grok.shell.addTableView(tj);

// Join types:
//   From SYNC_TYPE enum:
//      INNER
//      OUTER
//      LEFT
//      RIGHT
