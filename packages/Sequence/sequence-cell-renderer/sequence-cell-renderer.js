class SequenceCellRenderer extends DG.GridCellRenderer {

  get cellType() { return 'nucleotides'; }

  render(g, x, y, w, h, cell, style) {
    g.font = '13px monospace';
    g.textBaseline = 'middle';
    let s = cell.cell.value;

    for (let i = 0; i < s.length; i++) {
      g.fillStyle = DG.Color.toHtml(DG.Color.getCategoricalColor(s.charCodeAt(i)));
      g.fillText(s[i], 2 + x + i * 9, y + h / 2);
    }
  }
}