class UsageAnalysisPackage extends GrokPackage {

    init() { }

    //tags: app
    startApp() {
        let view = gr.newView('Usage');
        let acc = ui.accordion();

        function addPane(paneName, queryName, f) {
            acc.addPane((paneName), () => {
                let host = ui.div([], 'usage-analysis-card');
                gr.query(queryName, {'date': 'today'}).then((t) => host.appendChild(f(t)));
                return host;
            });
        }


        acc.addPane(('Users'), () => {
            let host = ui.div();
            gr.query('UniqueUsersByDate', {'date': 'today'})
                .then(t => {
                    let ids = Array.from(t.getCol('id').values());
                    gr.dapi.getEntities(ids).then((users) => host.appendChild(ui.list(users)));
                });
            return host;
        });

        addPane('Usage', 'EventsOnDate', (t) => Viewer.scatterPlot(t).root);
        addPane('Errors', 'ErrorsOnDate', (t) => Viewer.grid(t).root);
        addPane('Event Types', 'EventsSummaryOnDate', (t) => Viewer.barChart(t).root);
        addPane('Test Tracking', 'ManualActivityByDate', (t) => Viewer.grid(t).root);

        view.root.appendChild(acc.root);
    }
}