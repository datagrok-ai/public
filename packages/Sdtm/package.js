class SDTMPackage extends GrokPackage {
    //tags: app
    startApp(context) {
        let emptyTable = DataFrame.create();

        let limit = ui.intInput('Limit', 100);

        let view = View.create();
        view.name = 'SDTM:LB:Preview';
        view.description = 'SDTM LB domain viewer';
        view.root.className = 'grok-view grok-table-view stdm-result';

        function removeChildren(node) {
            while (node.firstChild)
                node.removeChild(node.firstChild);
        }

        function updatePreview(t) {
            let preview = Grid.create(t);
            removeChildren(view.root);
            view.root.appendChild(preview.root);
            ui.setUpdateIndicator(view.root, false);
        }

        updatePreview(emptyTable);
        ui.setUpdateIndicator(view.root, true);

        function selector(query, itemsColumn, caption, updateHandler) {
            let host = ui.divV([], 'stdm-selector-host');
            host.appendChild(ui.loader());
            grok.query(query, {}).then(t => {
                let items = t.getCol(itemsColumn).toList();
                let input = TagEditor.create();
                input.addTag(items[0]);
                let add = ui.divH([ui.span([caption]),
                    ui.iconFA('plus', () => {
                    Menu.popup().items(items, (item) => input.addTag(item)).show();
                })], 'stdm-selector');
                removeChildren(host);
                host.appendChild(add);
                host.appendChild(input.root);

                function update() {
                    ui.setUpdateIndicator(view.root, true);
                    if (host.childNodes.length > 2)
                        host.removeChild(host.childNodes[2]);
                    updateHandler(input.tags, host);
                }

                input.onChanged(() => {
                    update();
                });

                update();
            });
            return host;
        }

        let studiesParam = [];
        let testsParam = [];

        function aggregate(t) {
            grok.setVar('table', t);
            return grok.scriptSync('Aggregate(table, pivots = ["lbtest"], aggregations = ["values(lborres)"], groupByFields = ["studyid", "usubjid", "lbdy"])');
        }

        let selectors = selector('SDTM:Studies', 'Study', 'Studies', (studies, host) => {
            if (studies.length > 0) {
                studiesParam = studies;
                let testsSelector = selector('SDTM:Tests', 'Test', 'Tests', (tests) => {
                    testsParam = tests;
                    grok.query('SDTM:ResultsPreview', {
                        'studies': (studiesParam.length > 0) ? 'in (' + studiesParam.join(',') + ')' : '',
                        'tests': (testsParam.length > 0) ? 'in (' + testsParam.join(',') + ')' : '',
                        'limit': limit.value
                    }).then(t => {
                        updatePreview(aggregate(t));
                    });
                });
                host.appendChild(testsSelector);
            } else
                updatePreview(emptyTable);
        });

        function playIcon(handler) {
            let i = document.createElement('i');
            i.classList.add('grok-icon');
            i.classList.add('fas');
            i.classList.add('fa-play');
            i.style.color = '#3B9314';
            i.addEventListener("click", handler);
            return i;
        }

        let play = playIcon(() => {
            grok.callFunc('SDTM:Results', {
                'studies': (studiesParam.length > 0) ? 'in (' + studiesParam.join(',') + ')' : '',
                'tests': (testsParam.length > 0) ? 'in (' + testsParam.join(',') + ')' : '',
            }, true).then(t => {
                t.name = 'SDTM:LB:Results';
                grok.addTableView(aggregate(t));
            });
        });

        view.setRibbonPanels([[play]]);

        view.toolbox = ui.div([
            limit.root,
            selectors
        ], 'stdm-controls,pure-form');

        grok.addView(view);
    }
}
