export class FlagCellRenderer extends DG.GridCellRenderer {

  get name() {
    return 'Flag cell renderer';
  }

  get cellType() {
    return 'flag';
  }

  get defaultWidth() {
    return 50;
  }

  get defaultHeight() {
    return 50;
  }

  render(g, x, y, w, h, gridCell, cellStyle) {
    g.fillStyle = 'black';
    g.fillText('flag', x, y);
  }
}
