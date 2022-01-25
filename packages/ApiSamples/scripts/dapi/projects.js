let project = DG.Project.create();
let table = grok.data.demo.demog();
let tableInfo = table.getTableInfo();
project.addChild(tableInfo);

await grok.dapi.tables.uploadDataFrame(table);
await grok.dapi.tables.save(tableInfo);
await grok.dapi.projects.save(project);