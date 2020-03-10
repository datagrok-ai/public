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
        grok.callQuery('EnamineStore:Search', {
            'code': `search_${smiles}_${mode}`,
            'currency': currency
        }, true, 100).then(fc => {
            host.removeChild(host.firstChild);
            console.log(fc.getParamValue('stringResult'));
            let data = JSON.parse(fc.getParamValue('stringResult'))['data'];
            if (data === null) {
                host.appendChild(ui.divText('No matches'));
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
                ui.tooltip(mol, ui.tableFromMap(props));
                mol.addEventListener('click', function() {
                    window.open(comp['productUrl'], '_blank');
                });
                host.appendChild(mol);
            }
            host.style.overflowY = 'auto';
        });
        return host;
    }
}
