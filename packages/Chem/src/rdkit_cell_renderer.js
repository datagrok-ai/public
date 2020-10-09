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
        mol.draw_to_canvas_with_offset(g.canvas, x, -y, w, h);
        mol.delete();
    }
}