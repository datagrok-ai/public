// Creating HTML-driven tables
// For high-performance grid, check out samples under /grid

grok.shell.newView('tables', [
    ui.tableFromMap({
        user: grok.shell.user.toMarkup(),
        project: grok.shell.project.toMarkup(),
        time: new Date(),
    })
]);
