class EnamineStorePackage extends GrokPackage {

    //name: Enamine Store
    //description: Enamine Store Samples
    //tags: panel, widgets
    //input: string smiles {semType: Molecule}
    //output: widget result
    //condition: true
    enamineStore(smiles) {
        // sub, sim, exact
        let host = ui.div([ui.loader()]);
        this.storeSearch(smiles,  'sim', t => {
            let grid = Grid.create(t);
            grid.root.style.width = '350px';
            grid.root.style.height = '400px';
            while (host.firstChild)
                host.removeChild(host.firstChild);
            host.appendChild(grid.root);
        })
        return new Widget(host);
    }

    storeSearch(smiles, mode, handler) {
        grok.query('EnamineStore:API', {
            'code': `search_${smiles}_${mode}`,
            'currency': 'USD'
        }, true, 100).then(handler);
    }
}
