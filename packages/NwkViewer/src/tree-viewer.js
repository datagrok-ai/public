export class PhyloTreeViewer extends DG.JsViewer {
  constructor() {
    super();
    this.radialLayout = this.bool('radialLayout', true);
    this.fontSize = this.string('fontSize', '6px');
    this.defaultSize = 400;
    this.root.style = 'position: absolute; left: 0; right: 0; top: 0; bottom: 0;';
    this.nodeSelectionSourceColumnName = this.string('nodeSelectionSourceColumnName', 'node');
  }

  onTableAttached() {
    this.newick = this.dataFrame.getTag('.newick');
    this.parsedNewick = JSON.parse(this.dataFrame.getTag('.newickJson'));
    this.nodeSourceColumn = this.dataFrame.col(this.nodeSelectionSourceColumnName);
    this.subs.push(this.dataFrame.onCurrentRowChanged.subscribe(() => this.render()));
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
      .font_size(parseInt(this.fontSize))
      .radial(this.radialLayout);

    tree(this.parsedNewick).layout();

    if (this.nodeSourceColumn) {
      const node = this.nodeSourceColumn.get(this.dataFrame.currentRow.idx);
      if (node) {
        tree.modify_selection(node => false);
        tree.modify_selection(tree.path_to_root(tree.get_node_by_name(node)));
      }
    }
  }
}
