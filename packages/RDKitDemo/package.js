class RDKitDemoPackage extends DG.Package {

    /** Guaranteed to be executed exactly once before the execution of any function below */
    async init() {
        await initRDKit();
    }

    _svgDiv(mol) {
        let root = ui.div();
        var svg = mol.get_svg();
        root.innerHTML = svg;
        return root;
    }

    //name: getCLogP
    //input: string smiles {semType: Molecule}
    //output: double cLogP
    getCLogP(smiles) {
        var mol = Module.get_mol(smiles);
        return JSON.parse(mol.get_descriptors()).CrippenClogP;
    }

    //name: RDKitInfo
    //tags: panel, widgets
    //input: string smiles {semType: Molecule}
    //output: widget result
    rdkitInfoPanel(smiles) {
        var mol = Module.get_mol(smiles);
        return new DG.Widget(ui.divV([
            this._svgDiv(mol),
            ui.divText(`${this.getCLogP(smiles)}`)
        ]));
    }
}
