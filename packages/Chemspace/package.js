class ChemspacePackage extends GrokPackage {

    //tags: app
    //name: Chemspace
    async startApp(context) {
        let token = await ChemspacePackage.getApiToken();
        ChemspacePackage.queryMultipart(token).then(t => {
            gr.addTableView(t);
        })

        // let molecule = ui.moleculeInput('', 'c1ccccc1O');
        // let searchMode = ui.choiceInput('Mode', 'Similar', ['Exact', 'Similar', 'Substructure']);
        // let currency = ui.choiceInput('Currency', 'USD', ['USD', 'EUR']);
        // let similarity = ui.choiceInput('Similarity', '0.8', ['0.2', '0.4', '0.6', '0.8']);
        // let catalog = ui.choiceInput('Catalog', '', ['', 'BB', 'SCR', 'REAL']);
        // let filtersHost =  ui.div([molecule.root, searchMode.root, currency.root, similarity.root, catalog.root],
        //     'chemspace-controls,pure-form');
        //
        // let emptyTable = DataFrame.create();
        // let view = grok.addTableView(emptyTable);
        // view.name = 'Chemspace';
        // view.basePath = '';
        // view.description = 'Chemspace search viewer';
        // view.root.className = 'grok-view grok-table-view chemspace';
        //
        // function update() {
        //     ui.setUpdateIndicator(view.root, true);
        //     grok.callQuery('Chemspace:Search', {
        //         'code': `search_${molecule.value}_${ChemspacePackage.searchModeToCommand(searchMode.value)}`,
        //         'currency': currency.value,
        //         'sim': parseFloat(similarity.value),
        //         'mode': catalog.value
        //     }, true, 100).then(fc => {
        //         let data = JSON.parse(fc.getParamValue('stringResult'))['data'];
        //         view.dataFrame = data !== null ? ChemspacePackage.dataToTable(data, 'chemspace') : emptyTable;
        //         ui.setUpdateIndicator(view.root, false);
        //     });
        // }
        //
        // update();
        //
        // molecule.onChanged(() => update());
        // searchMode.onChanged(() => {
        //     similarity.enabled = searchMode.value === 'Similar';
        //     update();
        // });
        // currency.onChanged(() => update());
        // similarity.onChanged(() => update());
        // catalog.onChanged(() => update());
        //
        // let acc = view.toolboxPage.accordion;
        // acc.addPane('Chemspace', () => filtersHost, true, acc.panes[0]);
    }

    //name: Chemspace
    //description: Chemspace Samples
    //tags: panel, widgets
    //input: string smiles {semType: Molecule}
    //output: widget result
    //condition: true
    chemspacePanel(smiles) {
        let panels = ui.div([
            this.createSearchPanel('Exact', smiles),
            this.createSearchPanel('Similar', smiles),
            this.createSearchPanel('Substructure', smiles)
        ]);
        return new Widget(panels);
    }

    //description: Creates search panel
    createSearchPanel(panelName, smiles) {
        const currency = 'USD';
        let headerHost = ui.divH([ui.h2(panelName)], 'chemspace-panel-header');
        let compsHost = ui.divH([ui.loader()]);
        let panel = ui.divV([headerHost, compsHost], 'chemspace-panel');
        grok.callQuery('EnamineStore:Search', {
            'code': `search_${smiles}_${ChemspacePackage.searchModeToCommand(panelName)}`,
            'currency': currency
        }, true, 100).then(fc => {
            compsHost.removeChild(compsHost.firstChild);
            let data = JSON.parse(fc.getParamValue('stringResult'))['data'];
            if (data === null) {
                compsHost.appendChild(ui.divText('No matches'));
                return;
            }
            for (let comp of data) {
                let smiles = comp['smile'];
                let mol = ui.svgMol(smiles, 150, 75);
                let id = comp['Id'];
                let props = {
                    'ID': id,
                    'Formula': comp['formula'],
                    'MW': comp['mw'],
                    'Availability': comp['availability'],
                    'Delivery': comp['deliveryDays']
                };
                for (let pack of comp['packs'])
                    props[`${pack['amount']} ${pack['measure']}`] = `${pack['price']} ${currency}`;
                ui.tooltip(mol, ui.divV([ui.tableFromMap(props), ui.divText('Click to open in the store.')]));
                mol.addEventListener('click', function() {
                    window.open(comp['productUrl'], '_blank');
                });
                compsHost.appendChild(mol);
            }
            headerHost.appendChild(ui.iconFA('arrow-square-down', () =>
                grok.addTableView(ChemspacePackage.dataToTable(data, `EnamineStore ${panelName}`)), 'Open compounds as table'));
            compsHost.style.overflowY = 'auto';
        }).catch(err => {
            compsHost.removeChild(compsHost.firstChild);
            let div = ui.divText('No matches');
            ui.tooltip(div, `${err}`);
            compsHost.appendChild(div);
        });
        return panel;
    }

    // description: Converts JSON data into DataFrame
    static dataToTable(data, name) {
        let columns = [
            Column.fromStrings('smiles', data.map(comp => comp['smile'])),
            Column.fromStrings('ID', data.map(comp => comp['Id'])),
            Column.fromStrings('Formula', data.map(comp => comp['formula'])),
            Column.fromFloat32Array('MW', new Float32Array(data.map(comp => comp['mw']))),
            Column.fromInt32Array('Availability', new Int32Array(data.map(comp => comp['availability']))),
            Column.fromStrings('Delivery', data.map(comp => comp['deliveryDays']))
        ];
        let currency = null;
        let packsArrays = new Map();
        for (let n = 0; n < data.length; n++) {
            let packs = data[n]['packs'];
            for (let m = 0; m < packs.length; m++) {
                let pack = packs[m];
                let name = `${pack['amount']} ${pack['measure']}`;
                if (!packsArrays.has(name))
                    packsArrays.set(name, new Float32Array(data.length));
                packsArrays.get(name)[n] = pack['price'];
                if (currency === null && pack['currencyName'] !== null)
                    currency = pack['currencyName'];
            }
        }
        for (let name of packsArrays.keys()) {
            let column = Column.fromFloat32Array(name, packsArrays.get(name));
            column.semType = 'Money';
            column.setTag('format', `money(${currency === 'USD' ? '$' : '€'})`);
            columns.push(column);
        }
        let table = DataFrame.fromColumns(columns);
        table.name = name;
        return table;
    }

    //description: Converts search mode friendly name to command
    static searchModeToCommand(name) {
        const dict = {'Exact': 'exact', 'Similar': 'sim', 'Substructure': 'sub'};
        return dict[name];
    }

    //description: Gets access token
    static async getApiToken() {
        let t = await grok.query('Chemspace:AuthToken', null, true, 100);
        return t.get('access_token', 0);
    }

    //description: Perform query with multipart form data
    static queryMultipart(token) {
        // TODO: Deprecate after WebQuery 'multipart/form-data' support
        return new Promise(function (resolve, reject) {
            let host = 'https://api.chem-space.com';
            let formData = new FormData('SMILES=СOc1ccc(cc1)-c1nnc(SC)o1');
            let xhr = new XMLHttpRequest();
            xhr.setRequestHeader('Authorization', `Bearer ${token}`);
            xhr.setRequestHeader('Accept', 'application/json; version=2.6');
            xhr.open('POST', `${host}/v2/search/smiles`);
            xhr.onload = function () {
                if (this.status >= 200 && this.status < 300)
                    resolve(DataFrame.fromJson(xhr.responseText));
                else
                    reject({status: this.status, statusText: xhr.statusText});
            };
            xhr.onerror = function () {
                reject({status: this.status, statusText: xhr.statusText});
            };
            xhr.send(formData);
        });

    }
}
