/**
 * RDKit-based molecule cell renderer.
 * */
class RDKitCellRenderer extends DG.GridCellRenderer {

    constructor() {
        super();
        this.molCache = new DG.LruCache();
        this.molCache.onItemEvicted = function (mol) { mol.delete(); };
    }

    get name() { return 'RDKit cell renderer'; }

    get cellType() { return DG.SEMTYPE.MOLECULE; }
    
    get defaultWidth() { return 200; }

    get defaultHeight() { return 150; }

    render(g, x, y, w, h, gridCell, cellStyle) {
        
        let value = gridCell.cell.value;
        if (value == null || value === '')
            return;

        let mol = this.molCache.getOrCreate(value, (s) => Module.get_mol(s));

        if (!mol.is_valid())
            return;
        
        let drawMolecule = function (rdkitMol) {
            //rdkitMol.draw_to_canvas_with_offset(g.canvas, x, -y, w, h);         
            const opts = {
              "clearBackground": false,
              "offsetx": Math.floor(x),
              "offsety": -Math.floor(y),
              "width": Math.floor(w),
              "height": Math.floor(h)
            }
            mol.draw_to_canvas_with_highlights(g.canvas, JSON.stringify(opts));
        }
    
        let molIsInMolBlock = function(molString, rdkitMol) {
            const smilesMolString = rdkitMol.get_smiles();
            if (smilesMolString === molString)
                return false;
            const cxsmilesMolString = rdkitMol.get_cxsmiles();
            if (cxsmilesMolString === molString)
                return false;
            const inchiMolString = rdkitMol.get_inchi();
            if (inchiMolString === molString)
                return false;
            return true;
        }
        
        let drawMoleculeWithScaffold = function(scaffoldMolString, rdkitMol) {
            let scaffoldMol = Module.get_mol(scaffoldMolString);
            if (!scaffoldMol.is_valid()) {
                drawMolecule(rdkitMol);
                return;
            }
            if (molIsInMolBlock(scaffoldMolString, scaffoldMol)) {
                const substructJson = rdkitMol.get_substruct_match(scaffoldMol);
                if (substructJson !== '{}') {
                    rdkitMol.generate_aligned_coords(scaffoldMol, true);
                }
                drawMolecule(rdkitMol);
            }
            scaffoldMol.delete();
        }
    
        const molCol = gridCell.tableColumn.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
        let singleScaffoldMolString = molCol ? molCol.tags['chem-scaffold'] : null;
        
        if (singleScaffoldMolString) {
            drawMoleculeWithScaffold(singleScaffoldMolString, mol);
        } else {
    
            let df = gridCell.tableColumn.dataFrame;
            let rowScaffoldCol = null;
            for (let j = 0; j < df.columns.length; ++j) {
                let col = df.columns.byIndex(j);
                let tags = col.tags;
                if (tags && tags['row-scaffold']) {
                    rowScaffoldCol = col;
                    break;
                }
            }
            
            if (rowScaffoldCol == null || rowScaffoldCol.name === gridCell.tableColumn.name) {
                // regular drawing
                drawMolecule(mol);
            } else {
                let idx = gridCell.tableRowIndex;
                let scaffoldMolString = df.get(rowScaffoldCol.name, idx);
                drawMoleculeWithScaffold(scaffoldMolString, mol);
            }
        }
    }
}