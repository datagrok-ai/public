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
        let scaffold = gridCell.tableColumn.dataFrame.columns.byName('smiles').tags['chem-scaffold'];

        if (scaffold === null)
            // regular drawing
            mol.draw_to_canvas_with_offset(g.canvas, x, -y, w, h);
        else {
            // align to scaffold
            let qmol = Module.get_mol(scaffold);
            if (mol.is_valid() && qmol.is_valid()) {
                const mdetails = mol.get_substruct_match(qmol);
                const match = JSON.parse(mdetails);
                const useCoordgen = true;

                if (match.atoms && match.atoms.length)
                    mol.generate_aligned_coords(qmol, useCoordgen);
                // if (match.atoms && match.atoms.length)
                //     draw_with_highlights(mol, match);
                // else
                    mol.draw_to_canvas_with_offset(g.canvas, x, -y, w, h);
            }
        }

        mol.delete();
    }
}