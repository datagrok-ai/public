class RDKitDemoPackage extends GrokPackage {

    async init() {
        console.log('async init started');
        await initRDKit();
        console.log('async init finished');
    }

    //name: RDKitInfo
    //tags: panel, widgets
    //input: string smiles {semType: Molecule}
    //output: widget result
    rdkitInfoPanel(smiles) {
        var mol = Module.get_mol(smiles);
        let root = ui.div();
        var svg = mol.get_svg();
        root.innerHTML = svg;
        return new Widget(root);
    }


    //name: ShowMol
    //input: string smiles {semType: Molecule}
    showMol(smiles) {
        ui.dialog('Molecule')
            .add(this.webRoot)
            .add(ui.h1(smiles))
            .show();
        this._foo(smiles);
    }


    //name: AsyncTest
    //output: string userInfo
    async asyncTest() {
        var userInfo = await fetch('https://api.github.com/users/skalkin');
        return userInfo;
    }

    //name: Test
    async test() {
        //await initRDKit();

        var mol = Module.get_mol("c1ccccc1O");
        var descrs = JSON.parse(mol.get_descriptors());
        var clogP = descrs.CrippenClogP;

        ui.dialog('Molecule')
            .add(this.webRoot)
            .add(`LogP: ${clogP}`)
            .show();
    }
}
