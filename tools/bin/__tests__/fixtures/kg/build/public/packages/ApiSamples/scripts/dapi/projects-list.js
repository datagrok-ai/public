//api: DG.HttpDataSource.filter, DG.HttpDataSource.list
//help-url: https://datagrok.ai/help/datagrok/project
// Lists demo projects; the canonical sample for ~domains/bio project browsing.

let demo = grok.dapi.projects.filter('#demo');
let projects = await demo.list();
grok.shell.info(projects.length);
