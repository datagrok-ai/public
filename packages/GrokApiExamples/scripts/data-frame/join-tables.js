// The following example joins two tables, using the specified key and value columns.
//
// https://datagrok.ai/help/dialogs/join-tables

let t1 = grok.DataFrame.fromCsv(
`id1,id2,v1,v2
a1,b4,0,0
a2,b3,1,1
a3,b2,2,2
a4,b1,3,3`);
t1.name = 't1';

let t2 = grok.DataFrame.fromCsv(
`id3,id4,v3,v4
a1,b4,5,5
a2,b3,6,6
a3,b2,7,7
a4,b1,8,8
a5,b0,9,9
a6,b5,10,10`);
t2.name = 't2';

let tj = grok.joinTables(t1, t2, ['id1', 'id2'], ['id3', 'id4'], ['v1', 'v2'], ['v3', 'v4'], JOIN_TYPE_INNER, false);

grok.addTableView(t1);
grok.addTableView(t2);
grok.addTableView(tj);

// Join types:
//
//   JOIN_TYPE_INNER
//   JOIN_TYPE_OUTER
//   JOIN_TYPE_LEFT
//   JOIN_TYPE_RIGHT
