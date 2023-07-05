let card = (project) => ui.cards.summary(
  ui.image(project.pictureUrl, 50, 45),
  [
    ui.link(project.name, project),
    project.description
  ]
);

grok.dapi.projects
  .filter('#demo')
  .list({pageSize: 5})
  .then((projects) => grok.shell.o = ui.divV(projects.map(card)));