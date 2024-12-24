//tags: Project
//help-url: https://datagrok.ai/help/datagrok/project
// gets #demo project list and shows them in the view

let view = grok.shell.newView('projects');
let projects = await grok.dapi.projects
  .filter('#demo')
  .list();
view.append(ui.div(projects.map((p) => ui.renderCard(p)), 'grok-gallery-grid'));
