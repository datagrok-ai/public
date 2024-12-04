let card = (project) => ui.cards.summary(
  ui.image(project.pictureUrl, 50, 45),
  [
    ui.link(project.name, project),
    project.description
  ]
);

grok.shell.windows.showProperties = true;
let projects = await grok.dapi.projects
  .filter('#demo')
  .list({pageSize: 5});
grok.shell.o = ui.divV(projects.map(card));