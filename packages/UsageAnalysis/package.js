class UsageAnalysisPackage extends GrokPackage {

    init() { }

    //tags: app
    startApp() {
        let acc = null;
        let view = gr.newView('Usage');
        let results = ui.div();
        let panesExpanded = {};

        let date = ui.stringInput('Date', 'today');
        date.root.classList.add('usage-analysis-date', 'pure-form');
        date.addPatternMenu('datetime');

        let users = TagEditor.create();
        users.acceptsDragDrop = (x) => x instanceof User;
        users.doDrop = (user) => users.addTag(user.login);

        let addUser = ui.div([ui.iconFA('plus', () => {
            gr.dapi.users.order('login').list().then((allUsers) => {
                Menu.popup()
                    .items(allUsers.map(u => u.login), (item) => users.addTag(item))
                    .show();
            });
        })], 'usage-analysis-users-plus');

        function showUsage() {
            if (acc !== null)
                for (let p of acc.panes)
                    panesExpanded[p.name] = p.expanded;

            acc = ui.accordion();

            function addPane(paneName, queryName, f, supportUsers = true) {
                if (!(paneName in panesExpanded))
                    panesExpanded[paneName] = false;
                acc.addPane(paneName, () => {
                    let host = ui.div([], 'usage-analysis-card');
                    host.appendChild(ui.loader());
                    let params = {'date': date.value};
                    let selectedUsers = users.tags;
                    if (supportUsers && selectedUsers.length !== 0) {
                        queryName += 'AndUsers';
                        params['users'] = selectedUsers;
                    }
                    gr.query('UsageAnalysis:' + queryName, params).then((t) => {
                        host.removeChild(host.firstChild);
                        host.appendChild(f(t));
                    });
                    return host;
                }, panesExpanded[paneName]);
            }

            while (results.firstChild)
                results.removeChild(results.firstChild);

            acc.addPane('Unique users', () => {
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
            addPane('Test Tracking', 'ManualActivityByDate', (t) => Viewer.grid(t).root, false);

            results.appendChild(acc.root);
        }

        date.onChanged(this.debounce(showUsage, 750));
        users.onChanged(showUsage);

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
