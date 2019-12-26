class RDKitDemoPackage extends GrokPackage {

    /** Guaranteed to be executed exactly once before the execution of any function below */
    async init() {
        await initRDKit();
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
            .add(ui.h1(smiles))
            .show();
    }


    //name: Test
    async test() {
        var mol = Module.get_mol("c1ccccc1O");
        var descrs = JSON.parse(mol.get_descriptors());
        var clogP = descrs.CrippenClogP;

        ui.dialog('Molecule')
            .add(`LogP: ${clogP}`)
            .show();
    }
}
