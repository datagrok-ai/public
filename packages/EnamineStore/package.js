class EnamineStorePackage extends GrokPackage {

    //name: Enamine Store
    //description: Enamine Store Samples
    //tags: panel, widgets
    //input: string smiles {semType: Molecule}
    //output: widget result
    //condition: true
    enamineStore(smiles) {
        let host = ui.div([
            ui.h2('Exact'),
            this.storeSearch(smiles, 'exact'),
            ui.h2('Similar'),
            this.storeSearch(smiles, 'sim'),
            ui.h2('Substructure'),
            this.storeSearch(smiles, 'sub')
        ]);
        host.firstChild.style.marginTop = '6px';
        return new Widget(host);
    }

    storeSearch(smiles, mode) {
        const currency = 'USD';
        let host = ui.divH([ui.loader()]);
        grok.query('EnamineStore:Search', {
            'code': `search_${smiles}_${mode}`,
            'currency': currency
        }, true, 100).then(t => {
            host.removeChild(host.firstChild);
            if (t.rowCount === 0)
                host.appendChild(ui.divText('No matches'));
            for (let n = 0; n < t.rowCount; n++) {
                let smiles = t.get('data/smile', n);
                let mol = ui.svgMol(smiles, 150, 75);
                let id = t.get('data/Id', n);
                let tooltip = ui.tableFromMap({
                    'ID': id,
                    'Name': t.get('data/name', n),
                    'Formula': t.get('data/formula', n),
                    'MW': t.get('data/mw', n),
                    'Availability': t.get('data/availability', n),
                    'Delivery': t.get('data/deliveryDays', n),
                });
                // let packsHost = ui.divV([ui.h3('Packs'), ui.loader()]);
                // grok.query('EnamineStore:Packs', {'code': id, 'currency': currency}, true, 100).then(t => {
                //     packsHost.removeChild(packsHost.lastChild);
                //     let packs = {};
                //     for (let n = 0; n < t.rowCount; n++)
                //         packs[`${t.get('amount', n)} ${t.get('measure', n)}`] = `${t.get('price', n)} ${currency}`;
                //     packsHost.appendChild(ui.tableFromMap(packs));
                // });
                // tooltip.appendChild(packsHost);
                ui.tooltip(mol, tooltip);
                mol.addEventListener('click', function() {
                    window.open(t.get('data/productUrl', n), '_blank');
                });
                host.appendChild(mol);
            }
            host.style.overflowY = 'auto';
        });
        return host;
    }
}
