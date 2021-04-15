export class PhyloTreeViewer extends DG.JsViewer {
  constructor() {
    super();
    this.radialLayout = this.bool('radialLayout', true);
    this.fontSize = this.string('fontSize', '6px');
    this.defaultSize = 400;
  }

  onTableAttached() {
    this.newickCol = this.dataFrame.columns.bySemType('newick');
    this.subs.push(this.dataFrame.onCurrentRowChanged.subscribe((_) => this.render()));
    this.render();
  }

  render() {
    $(this.root).empty();

    if (this.newickCol == null) {
      this.root.appendChild(ui.divText('Newick column not found.', 'd4-viewer-error'));
      return;
    }

    const newick = this.newickCol.get(this.dataFrame.currentRow.idx);
    const svg = d3.select(this.root).append("svg");

    const tree = d3.layout.phylotree()
      .svg(svg)
      .options({
        'left-right-spacing': 'fit-to-size',
        'top-bottom-spacing': 'fit-to-size',
      })
      .size([
        this.root.parentElement.clientHeight || this.defaultSize,
        this.root.parentElement.clientWidth || this.defaultSize
      ])
      .font_size(Number.parseInt(this.fontSize))
      .radial(this.radialLayout);
  
    tree(d3.layout.newick_parser(newick)).layout();
    this.root.firstChild.style = 'position: absolute; left: 50%; top: 50%; transform: translate(-50%, -50%);';
  }
}
