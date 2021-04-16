export class PhyloTreeViewer extends DG.JsViewer {
  constructor() {
    super();
    this.radialLayout = this.bool('radialLayout', true);
    this.fontSize = this.string('fontSize', '6px');
    this.defaultSize = 400;
  }

  onTableAttached() {
    this.newick = this.dataFrame.getTag('.newick');
    this.render();
  }

  onPropertyChanged() { this.render(); }

  render() {
    $(this.root).empty();

    if (this.newick == null) {
      this.root.appendChild(ui.divText('Newick tag not found.', 'd4-viewer-error'));
      return;
    }

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
  
    tree(d3.layout.newick_parser(this.newick)).layout();
    this.root.style = 'position: absolute; left: 0; right: 0; top: 0; bottom: 0;';
  }
}
