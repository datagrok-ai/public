class ChemPackage extends DG.Package {
    
    //name: chemExportFunc
    //tags: autostart
    chemInit() { }
    
    /** Guaranteed to be executed exactly once before the execution of any function below */
    async init() {
        await initRDKit();
        console.log('RDKit initialized');
        Module.prefer_coordgen(true);
        this.STORAGE_NAME = 'rdkit_descriptors';
        this.KEY = 'selected';
    }

    //name: SubstructureFilter
    //description: RDKit-based substructure filter
    //tags: filter
    //output: filter result
    substructureFilter() {
        return new SubstructureFilter();
    }

    _svgDiv(mol) {
        let root = ui.div();
        root.innerHTML = mol.get_svg();
        return root;
    }

    //name: getCLogP
    //input: string smiles {semType: Molecule}
    //output: double cLogP
    getCLogP(smiles) {
        let mol = Module.get_mol(smiles);
        return JSON.parse(mol.get_descriptors()).CrippenClogP;
    }

    //name: RDKitInfo
    //tags: panel, widgets
    //input: string smiles {semType: Molecule}
    //output: widget result
    rdkitInfoPanel(smiles) {
        let mol = Module.get_mol(smiles);
        return new DG.Widget(ui.divV([
            this._svgDiv(mol),
            ui.divText(`${this.getCLogP(smiles)}`)
        ]));
    }

    //name: rdkitCellRenderer
    //tags: cellRenderer, cellRenderer-Molecule
    //meta-cell-renderer-sem-type: Molecule
    //output: grid_cell_renderer result
    rdkitCellRenderer() { return new RDKitCellRenderer(); }

    /*

    static _morganFP(smiles, fp_length = 128, fp_radius = 2) {
        let mol = Module.get_mol(smiles);
        let mfp = mol.get_morgan_fp(fp_radius, fp_length);
        mol.delete();
        return mfp;
    }

    //name: getFingerprints
    //input: column molColumn
    //output: column result [fingerprints]
    getFingerprints(molColumn) {
        let fps = molColumn.toList().map((smiles) => DG.BitSet.fromString(ChemPackage._morganFP(smiles)).d);
        return DG.Column.fromList('object', 'fingerprints', fps);
    }
    
    */
    
    //name: similarityScoring
    //input: column molStringsColumn
    //input: string molString
    //input: bool sorted
    //output: dataframe result
    similarityScoring(molStringsColumn, molString, sorted) {
        return chemSimilarityScoring(molStringsColumn, molString, {'sorted' : sorted});
    }
    
    //name: substructureSearch
    //input: column molStringsColumn
    //input: string molString
    //output: column result
    substructureSearch(molStringsColumn, molString) {
        return DG.Column.fromList('object', 'bitset',
          [chemSubstructureSearch(molStringsColumn, molString)]);
    }

    //name: molColumnPropertyPanel
    //input: column molColumn
    //tags: panel
    //output: widget result
    molColumnPropertyPanel(molColumn) { return getMolColumnPropertyPanel(molColumn); }

    //tags: app
    descriptorsApp(context) {
        let defaultSmiles = 'O=C1CN=C(c2ccccc2N1)C3CCCCC3';
        let sketcherValue = defaultSmiles;

        let windows = grok.shell.windows;
        windows.showToolbox = false;
        windows.showHelp = false;
        windows.showProperties = false;

        let table = DG.DataFrame.create();
        table.name = 'Descriptors';
        let view = grok.shell.addTableView(table);

        let dsDiv = ui.divV([], 'grok-prop-panel');
        dsDiv.appendChild(this.descriptorsWidget(defaultSmiles).root);

        let sketcher = grok.chem.sketcher((smiles, molfile) => {
            sketcherValue = smiles;
            RDKitDemoPackage.removeChildren(dsDiv);
            dsDiv.appendChild(this.descriptorsWidget(smiles).root);
        }, defaultSmiles);
        let addButton = ui.bigButton('ADD', async () => {
            this.getSelected().then(selected => {
                grok.chem.descriptors(DG.DataFrame.fromCsv(`smiles\n${sketcherValue}`), 'smiles', selected).then(t => {
                    let columnNames = table.columns.names();
                    if ((table.columns.length !== selected.length + 1) || selected.some(s => !columnNames.includes(s))) {
                        table = DG.DataFrame.create();
                        table.name = 'Descriptors';
                        view.dataFrame = table;
                        for (let col of t.columns.toList())
                            table.columns.addNew(col.name, col.type);
                    }
                    table.rows.addNew(t.columns.toList().map(c => c.get(0)));
                });
            });
        });
        addButton.style.marginTop = '12px';
        let skDiv = ui.divV([sketcher, addButton], 'grok-prop-panel,dlg-sketcher,pure-form');

        let skNode = view.dockManager.dock(skDiv, DG.DOCK_TYPE.RIGHT, null, 'Sketcher', 0.25);
        view.dockManager.dock(dsDiv, DG.DOCK_TYPE.DOWN, skNode, 'Descriptors', 0.5);

        grok.events.onViewRemoved.subscribe((v) => {
            if (v.name === view.name) {
                windows.showToolbox = true;
                windows.showHelp = true;
                windows.showProperties = true;
            }
        });
    }

    //name: Chem Descriptors
    //tags: panel, widgets
    //input: string smiles { semType: Molecule }
    //output: widget result
    descriptorsWidget(smiles) {
        let widget = new DG.Widget(ui.div());
        let result = ui.div();
        let selectButton = ui.bigButton('SELECT', async () => {
            RDKitDemoPackage.openDescriptorsDialog(await this.getSelected(), async (selected) => {
                await grok.dapi.userDataStorage.postValue(this.STORAGE_NAME, this.KEY, JSON.stringify(selected));
                update();
            });
        });
        selectButton.style.marginTop = '20px';

        let update = () => {
            RDKitDemoPackage.removeChildren(result);
            result.appendChild(ui.loader());
            this.getSelected().then(selected => {
                grok.chem.descriptors(DG.DataFrame.fromCsv(`smiles\n${smiles}`), 'smiles', selected).then(table => {
                    RDKitDemoPackage.removeChildren(result);
                    let map = {};
                    for (let descriptor of selected)
                        map[descriptor] = table.col(descriptor).get(0);
                    result.appendChild(ui.tableFromMap(map));
                });
            });
        }

        widget.root.appendChild(result);
        widget.root.appendChild(selectButton);

        update();

        return widget;
    }

    //description: Get selected descriptors
    async getSelected() {
        let str = await grok.dapi.userDataStorage.getValue(this.STORAGE_NAME, this.KEY);
        let selected = (str != null && str !== '') ? JSON.parse(str) : [];
        if (selected.length === 0) {
            selected = (await grok.chem.descriptorsTree())['Lipinski']['descriptors'].slice(0, 3).map(p => p['name']);
            await grok.dapi.userDataStorage.postValue(this.STORAGE_NAME, this.KEY, JSON.stringify(selected));
        }
        return selected;
    }

    //description: Open descriptors selection dialog
    static openDescriptorsDialog(selected, onOK) {
        grok.chem.descriptorsTree().then(descriptors => {
            let tree = ui.tree();
            tree.root.style.maxHeight = '400px';

            let groups = {};
            let items = [];

            for (let groupName in descriptors) {
                let group = tree.group(groupName, null, false);
                group.enableCheckBox();
                groups[groupName] = group;

                for (let descriptor of descriptors[groupName]['descriptors']) {
                    let item = group.item(descriptor['name'], descriptor);
                    item.enableCheckBox(selected.includes(descriptor['name']));
                    items.push(item);
                }
            }

            let clear = ui.button('NONE', () => {
                for (let g in groups) groups[g].checked = false;
                for (let i of items) i.checked = false;
            });

            ui.dialog('Chem Descriptors')
                .add(clear)
                .add(tree.root)
                .onOK(() => onOK(items.filter(i => i.checked).map(i => i.value['name'])))
                .show();
        });
    }

    //description: Removes all children from node
    static removeChildren(node) {
        while (node.firstChild)
            node.removeChild(node.firstChild);
    }
}