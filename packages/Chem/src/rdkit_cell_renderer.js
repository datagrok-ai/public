/**
 * RDKit-based molecule cell renderer.
 * */
class RDKitCellRenderer extends DG.GridCellRenderer {

    get name() { return 'RDKit cell renderer'; }

    get cellType() { return 'RDMolecule'; }

    render(g, x, y, w, h, gridCell, cellStyle) {
        
        let value = gridCell.cell.value;
        if (value == null || value === '')
            return;

        let mol = Module.get_mol(value);
        // let scaffold = gridCell.tableColumn.dataFrame.columns.byName('smiles').tags['chem-scaffold'];
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
        let smilesToMolBlock = function (smiles) {
            let qmol = Module.get_mol(smiles);
            let scaffoldInMolBlock = qmol.get_v3Kmolblock();
            qmol.delete();
            return scaffoldInMolBlock;
        }
        let drawAtomic = function (mol) {
            mol.draw_to_canvas_with_offset(g.canvas, x, -y, w, h);
        }
        let drawSimple = function (mol) {
            if (gridCell.tableColumn.name === 'scaffold') {
                // specially draw with MolBlock orientation
                let smol = Module.get_mol(smilesToMolBlock(value));
                mol.generate_aligned_coords(smol, true);
                drawAtomic(mol);
                smol.delete();
            } else {
                drawAtomic(mol);
            }
        }
        if (rowScaffoldCol == null) {
            // regular drawing
            drawSimple(mol);
        } else {
            let idx = gridCell.tableRowIndex;
            let scaffold = df.get(rowScaffoldCol.name, idx);
            // align to scaffold
            let smol = Module.get_mol(smilesToMolBlock(scaffold));
            try {
                if (mol.is_valid() && smol.is_valid()) {
                    const mdetails = mol.get_substruct_match(smol);
                    const match = JSON.parse(mdetails);
                    // draw_with_highlights(mol, match);
                    if (match.atoms && match.atoms.length) {
                        mol.generate_aligned_coords(smol, true);
                        drawAtomic(mol);
                    } else {
                        drawSimple(mol);
                    }
                }
            } catch {
                drawSimple(mol);
            }
            smol.delete();
        }
        mol.delete();
    }
}