class UsageAnalysisPackage extends GrokPackage {

    init() { }

    //tags: app
    startApp() {
        let view = gr.newView('Usage');
        var results = ui.div();
        let date = ui.stringInput('Date', 'today');
        date.addPatternMenu('datetime');

        function showUsage() {
            let acc = ui.accordion();

            function addPane(paneName, queryName, f) {
                acc.addPane((paneName), () => {
                    let host = ui.div([], 'usage-analysis-card');
                    gr.query(queryName, {'date': date.value}).then((t) => host.appendChild(f(t)));
                    return host;
                });
            }

            while (results.firstChild)
                results.removeChild(results.firstChild);

            acc.addPane(('Users'), () => {
                let host = ui.div();
                gr.query('UniqueUsersByDate', {'date': date.value})
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

            results.appendChild(acc.root);
        }

        date.onChanged(this.debounce(showUsage, 750));

        showUsage();

        view.append(ui.inputs([date]));
        view.root.appendChild(results);
    }

    debounce(fn, time) {
        let timeout;

        return function() {
            const functionCall = () => fn.apply(this, arguments);

            clearTimeout(timeout);
            timeout = setTimeout(functionCall, time);
        }
    }
}
