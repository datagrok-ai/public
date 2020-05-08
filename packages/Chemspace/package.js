class ChemspacePackage extends GrokPackage {

    static token = null;

    //tags: app
    //name: Chemspace
    async startApp(context) {
        const catalogToParam = {'BB': 'smarts/bb', 'SCR': 'smarts/sc', 'REAL': 'advanced/1'};
        const modeToParam = {'Similar': 'sim', 'Substructure': 'sub'};

        let token = await ChemspacePackage.getApiToken();

        let molecule = ui.moleculeInput('', 'c1ccccc1O');
        let mode = ui.choiceInput('Mode', 'Similar', ['Similar', 'Substructure']);
        let similarity = ui.choiceInput('Similarity', '0.6', ['0.2', '0.4', '0.6', '0.8']);
        let catalog = ui.choiceInput('Catalog', 'SCR', ['BB', 'SCR', 'REAL']);
        let filtersHost =  ui.div([molecule.root, mode.root, similarity.root, catalog.root],
            'chemspace-controls,pure-form');

        let emptyTable = DataFrame.create();
        let view = grok.addTableView(emptyTable);
        view.name = 'Chemspace';
        view.basePath = '';
        view.description = 'Chemspace search viewer';
        view.root.className = 'grok-view grok-table-view chemspace';

        function update() {
            ui.setUpdateIndicator(view.root, true);

            function setDataFrame(t) {
                view.dataFrame = t;
                let order = ['smiles', 'CS-id', 'id', 'similarity', 'iupac_name', 'molformula', 'mw',
                    'molweight', 'cas', 'hac', 'logp', 'rotb', 'hba', 'hbd',
                    'ring_count', 'fsp3', 'tpsa', 'mfcd', 'price_category', 'vendor_id', 'link'
                ];
                if (t.rowCount > 0) {
                    let names = t.columns.names();
                    order = order.filter(o => names.includes(o));
                    view.grid.columns.setOrder(order);
                }
                ui.setUpdateIndicator(view.root, false);
            }

            ChemspacePackage
                .queryMultipart(`search/${catalogToParam[catalog.value]}/${modeToParam[mode.value]}`,
                    molecule.value,
                    mode.value === 'Similar' ? {'simThreshold' : parseFloat(similarity.value) * 100} : null,
                    token)
                .then(t => setDataFrame(t))
                .catch(err => setDataFrame(emptyTable));
        }

        update();

        molecule.onChanged(() => update());
        mode.onChanged(() => {
            similarity.enabled = mode.value === 'Similar';
            update();
        });
        similarity.onChanged(() => update());
        catalog.onChanged(() => update());

        let acc = view.toolboxPage.accordion;
        acc.addPane('Chemspace', () => filtersHost, true, acc.panes[0]);
    }

    //name: Chemspace Samples
    //description: Chemspace Samples
    //tags: panel, widgets
    //input: string smiles {semType: Molecule}
    //output: widget result
    //condition: true
    chemspaceSamplesPanel(smiles) {
        let panels = ui.div();
        ChemspacePackage.getApiToken().then(() => {
            panels.appendChild(this.createSearchPanel('Similar', smiles));
            panels.appendChild(this.createSearchPanel('Substructure', smiles));
        });
        return new Widget(panels);
    }

    //description: Creates search panel
    createSearchPanel(panelName, smiles) {
        const searchModeToCommand = {'Similar': 'smarts/bb/sim', 'Substructure': 'smarts/bb/sub'};

        let headerHost = ui.divH([ui.h2(panelName)], 'chemspace-panel-header');
        let compsHost = ui.divH([ui.loader()]);
        let panel = ui.divV([headerHost, compsHost], 'chemspace-panel');

        ChemspacePackage
            .queryMultipart(`search/${searchModeToCommand[panelName]}`, smiles,
                panelName === 'Similar' ? {'simThreshold': 40} : null, ChemspacePackage.token)
            .then(t => {
                compsHost.removeChild(compsHost.firstChild);
                if (t.rowCount === 0) {
                    compsHost.appendChild(ui.divText('No matches'));
                    return;
                }

                function getTooltip(n) {
                    let props = {
                        'ID': t.get('CS-id', n),
                        'IUPAC': t.get('iupac_name', n),
                        'Formula': t.get('molformula', n),
                        'MW': t.get('molweight', n)
                    };
                    return ui.divV([ui.tableFromMap(props), ui.divText('Click to open in the store.')]);
                }

                for (let n = 0; n < Math.min(t.rowCount, 20); n++) {
                    let smiles = t.get('smiles', n);
                    let mol = ui.svgMol(smiles, 150, 75);
                    ui.tooltip(mol, () => getTooltip(n));
                    mol.addEventListener('click', function () {
                        window.open(t.get('link', n), '_blank');
                    });
                    compsHost.appendChild(mol);
                }
                headerHost.appendChild(ui.iconFA('arrow-square-down', () => {
                    t.name = `Chemspace ${panelName}`;
                    grok.addTableView(t);
                }, 'Open compounds as table'));
                compsHost.style.overflowY = 'auto';
            })
            .catch(err => {
                compsHost.removeChild(compsHost.firstChild);
                let div = ui.divText('No matches');
                ui.tooltip(div, `${err}`);
                compsHost.appendChild(div);
            });

        return panel;
    }

    //description: Gets access token
    static async getApiToken() {
        if (ChemspacePackage.token === null) {
            let t = await grok.query('Chemspace:AuthToken', null, true, 100);
            ChemspacePackage.token = t.get('access_token', 0);
        }
        return ChemspacePackage.token;
    }

    //description: Perform query with multipart form data
    static queryMultipart(path, smiles, params, token) {
        // TODO: Deprecate after WebQuery 'multipart/form-data' support
        return new Promise(function (resolve, reject) {
            let host = 'https://api.chem-space.com';
            let xhr = new XMLHttpRequest();
            let formData = new FormData();
            formData.append('SMILES', smiles);
            let queryParams = params !== null ? `?${Object.keys(params).map(key => key + '=' + params[key]).join('&')}` : '';
            xhr.open('POST', `${host}/v2/${path}${queryParams}`);
            xhr.setRequestHeader('Authorization', `Bearer ${token}`);
            xhr.setRequestHeader('Accept', 'application/json; version=2.6');
            xhr.onload = function () {
                if (this.status >= 200 && this.status < 300) {
                    let list = JSON.parse(xhr.responseText)['items'];
                    if (list.length > 0)
                        resolve(DataFrame.fromJson(JSON.stringify(list)));
                    else
                        reject();
                } else
                    reject();
            };
            xhr.onerror = () => reject();
            xhr.send(formData);
        });
    }
}
