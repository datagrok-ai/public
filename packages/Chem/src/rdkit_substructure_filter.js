/**
 * RDKit-based substructure filters that uses Datagrok's collaborative filtering.
 *
 * 1. On onRowsFiltering event, only FILTER OUT rows that do not satisfy this filter's criteria
 * 2. Call dataFrame.rows.requestFilter when filtering criteria changes.
 * */
class SubstructureFilter extends DG.Filter {

    constructor() {
        super();
        this.smiles = '';
        this.root = ui.div(null, 'grok-chem-substructure-filter');
        let sketcher = grok.chem.sketcher((smiles) => {
            this.smiles = smiles;
            this.dataFrame.temp.smarts = smiles;
            this.dataFrame.rows.requestFilter();
        })
        this.root.appendChild(sketcher);
    }

    attach(dFrame) {
        this.dataFrame = DG.toJs(dFrame);
        this.column = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
        this.dataFrame.onRowsFiltering.subscribe((_) => this.applyFilter());
    }

    applyFilter() {
        let subMol = Module.get_mol(this.smiles);

        for (let i of this.dataFrame.filter.getSelectedIndexes()) {
            let mol = Module.get_mol(this.column.get(i));
            let match = mol.get_substruct_match(subMol);
            if (match === "{}" )
                this.dataFrame.filter.set(i, false, false);
            mol.delete();
        }

        this.dataFrame.filter.fireChanged();
    }
}
