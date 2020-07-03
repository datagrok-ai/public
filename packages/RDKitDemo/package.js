class RDKitDemoPackage extends DG.Package {

    /** Guaranteed to be executed exactly once before the execution of any function below */
    async init() {
        await initRDKit();
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

    //tags: app
    descriptorsApp(context) {
        let windows = grok.shell.windows;
        windows.showToolbox = false;
        windows.showHelp = false;
        windows.showProperties = false;

        let sketcher = ui.div([grok.chem.sketcher(() => {})], 'dlg-sketcher,pure-form');
        let descriptors = ui.divV([], 'grok-prop-panel');
        let table = DG.DataFrame.create();
        table.name = 'Descriptors';

        descriptors.appendChild(this.descriptorsWidget('CC1=CC(=O)C=CC1=O').root);

        let view = grok.shell.addTableView(table);
        let skNode = view.dockManager.dock(sketcher, DG.DOCK_TYPE.RIGHT, null, 'Sketcher', 0.25);
        view.dockManager.dock(descriptors, DG.DOCK_TYPE.DOWN, skNode, 'Descriptors', 0.5);

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
        const STORAGE_NAME = 'rdkit_descriptors';
        const KEY = 'selected';

        async function getSelected() {
            let str = await grok.dapi.userDataStorage.getValue(STORAGE_NAME, KEY);
            let selected = (str != null && str !== '') ? JSON.parse(str) : [];
            if (selected.length === 0) {
                selected = (await grok.chem.descriptorsTree())['Lipinski']['descriptors'].slice(0, 3).map(p => p['name']);
                await grok.dapi.userDataStorage.postValue(STORAGE_NAME, KEY, JSON.stringify(selected));
            }
            return selected;
        }

        let widget = new DG.Widget(ui.div());
        let result = ui.div();
        let selectButton = ui.bigButton('SELECT', async () => {
            RDKitDemoPackage.openDescriptorsDialog(await getSelected(), async (selected) => {
                await grok.dapi.userDataStorage.postValue(STORAGE_NAME, KEY, JSON.stringify(selected));
                update();
            });
        });
        selectButton.style.marginTop = '20px';

        function update() {
            RDKitDemoPackage.removeChildren(result);
            result.appendChild(ui.loader());
            getSelected().then(selected => {
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
