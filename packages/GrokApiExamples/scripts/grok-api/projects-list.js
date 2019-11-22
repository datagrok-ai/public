// gets #demo project list and shows them in the view

let view = gr.newView('projects');
gr.dapi.projects
    .filter('#demo')
    .list()
    .then((projects) => view.append(ui.div(projects.map((p) => ui.renderCard(p.d)), 'grok-gallery-grid')));
