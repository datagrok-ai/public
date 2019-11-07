// Creating HTML-driven tables
// For high-performance grid, check out samples under /grid

gr.newView('tables', [
    ui.tableFromMap({
        user: gr.user.toMarkup(),
        project: gr.project.toMarkup(),
        time: new Date(),
    })
]);
