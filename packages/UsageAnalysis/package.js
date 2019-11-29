class UsageAnalysisPackage extends GrokPackage {

    init() { }

    //tags: app
    startApp() {
        let view = gr.newView('Usage');
        let results = ui.div();

        let date = ui.stringInput('Date', 'today');
        date.root.classList.add('usage-analysis-date', 'pure-form');
        date.addPatternMenu('datetime');

        let users = TagEditor.create();
        users.acceptsDragDrop = (x) => x instanceof User;
        users.doDrop = (user) => users.addTag(user.login);

        let addUser = ui.div([ui.iconFA('plus', () => {
            gr.dapi.users.list().then((allUsers) => {
                Menu.popup()
                    .item('Show info', () => gr.balloon.info('Info'))
                    .separator()
                    .items(allUsers.map(u => u.login), (item) => users.addTag(item))
                    .show();
            });
        })], 'usage-analysis-users-plus');

        function showUsage() {
            let acc = ui.accordion();

            function addPane(paneName, queryName, f) {
                acc.addPane(paneName, () => {
                    let host = ui.div([], 'usage-analysis-card');
                    host.appendChild(ui.loader());
                    gr.query('UsageAnalysis:' + queryName, {'date': date.value}).then((t) => {
                        host.removeChild(host.firstChild);
                        host.appendChild(f(t));
                    });
                    return host;
                });
            }

            while (results.firstChild)
                results.removeChild(results.firstChild);

            acc.addPane('Users', () => {
                let host = ui.div();
                host.appendChild(ui.loader());
                gr.query('UsageAnalysis:UniqueUsersByDate', {'date': date.value})
                    .then(t => {
                        let ids = Array.from(t.getCol('id').values());
                        gr.dapi.getEntities(ids).then((users) => {
                            host.removeChild(host.firstChild);
                            host.appendChild(ui.list(users));
                        });
                    });
                return host;
            });

            addPane('Usage', 'EventsOnDate', (t) => Viewer.scatterPlot(t, {'color': 'user'}).root);
            addPane('Errors', 'ErrorsOnDate', (t) => Viewer.grid(t).root);
            addPane('Event Types', 'EventsSummaryOnDate', (t) => Viewer.barChart(t).root);
            addPane('Error Types', 'ErrorsSummaryOnDate', (t) => Viewer.barChart(t).root);
            addPane('Test Tracking', 'ManualActivityByDate', (t) => Viewer.grid(t).root);

            results.appendChild(acc.root);
        }

        date.onChanged(this.debounce(showUsage, 750));

        showUsage();

        view.root.appendChild(results);

        let usersLabel = ui.div(null, 'usage-analysis-users-title');
        usersLabel.innerText = 'Users';
        view.toolbox = ui.divV([
            date.root,
            ui.divV([
                ui.divH([usersLabel, addUser]),
                users.root
            ], 'usage-analysis-users')

        ]);
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
