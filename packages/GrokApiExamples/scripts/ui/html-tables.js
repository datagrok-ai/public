// Creating HTML-driven tables
// For high-performance grid, check out samples under /grid

grok.newView('tables', [
    ui.tableFromMap({
        user: grok.user.toMarkup(),
        project: grok.project.toMarkup(),
        time: new Date(),
    })
]);
