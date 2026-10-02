//api: DG.HttpDataSource.filter, DG.HttpDataSource.order, DG.HttpDataSource.list, DG.HttpDataSource.count
//help-url: https://datagrok.ai/help/datagrok/project
// Lists #demo projects. Data-source verbs are immutable: each returns a new source, so a source
// kept in a variable stays a clean starting point.

let view = grok.shell.newView('projects');
let demo = grok.dapi.projects.filter('#demo').order('name');
let projects = await demo.list();
let recent = await demo.filter('#demo && createdOn > -1y').count();   // `demo` itself is unchanged
view.append(ui.divText(`${projects.length} demo projects, ${recent} created this year`));
view.append(ui.div(projects.map((p) => ui.renderCard(p)), 'grok-gallery-grid'));
