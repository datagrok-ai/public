// An example of using Dart's future as JavaScript's promises

let df = grok.data.demo.demog();
let view = grok.shell.addTableView(df);
 
const tableInfo = df.getTableInfo(); 
const table = await grok.dapi.tables.save(tableInfo); 


let tableSearchRes = await grok.data.openTable(table.id);

grok.shell.info(tableSearchRes.name); 
await grok.dapi.tables.delete(table);