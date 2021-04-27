export class PhyloTreeViewer extends DG.JsViewer {
  constructor() {
    super();
    this.radialLayout = this.bool('radialLayout', false);
    this.fontSize = this.string('fontSize', '9px');
    this.selection = this.string('selection', 'path to root & descendants', { choices: [
      'none', 'path to root', 'descendants', 'path to root & descendants'
    ]});
    this.tooltipOffset = 10;
    this.defaultSize = 400;
    this.root.style = 'position: absolute; left: 0; right: 0; top: 0; bottom: 0;';
    this.tree = d3.layout.phylotree();
    this.nodes = new Map();
  }

  onTableAttached() {
    this.newick = this.dataFrame.getTag('.newick');
    this.parsedNewick = JSON.parse(this.dataFrame.getTag('.newickJson'));
    this.nodeIdColumn = this.dataFrame.col('id');
    this.nodeNameColumn = this.dataFrame.col('node');
    this.parentNameColumn = this.dataFrame.col('parent');

    this.subs.push(this.dataFrame.onCurrentRowChanged.subscribe(() => this.render(false)));
    this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render()));
    this.render();

    d3.select(this.root).selectAll('.node > text')
      .on('mouseover', d => {
        ui.tooltip.show(
          ui.span([`${d.name}, parent: ${d.parent.name}`]),
          d3.event.x + this.tooltipOffset,
          d3.event.y + this.tooltipOffset);
      })
      .on('mouseout', () => ui.tooltip.hide());

    d3.select(this.root).selectAll('.internal-node')
      .on('mouseover', d => {
        ui.tooltip.show(
          ui.span([d.name + (d.name ? `, ` : '') + `children: ${d.children.length}`]),
          d3.event.x + this.tooltipOffset,
          d3.event.y + this.tooltipOffset);
      })
      .on('mouseout', () => ui.tooltip.hide());
  }

  onPropertyChanged() { this.render(); }

  render(redraw = true) {
    if (redraw) {
      $(this.root).empty();
      // this.nodes.clear();

      if (this.newick == null) {
        this.root.appendChild(ui.divText('Newick tag not found.', 'd4-viewer-error'));
        return;
      }

      const svg = d3.select(this.root).append("svg");

      this.tree
        .svg(svg)
        .options({
          'left-right-spacing': 'fit-to-size',
          'top-bottom-spacing': 'fit-to-size',
          zoom: true,
        })
        .size([
          this.root.parentElement.clientHeight || this.defaultSize,
          this.root.parentElement.clientWidth || this.defaultSize
        ])
        .font_size(parseInt(this.fontSize))
        .radial(this.radialLayout);

      this.tree(this.parsedNewick).layout();

      if (!this.nodeIdColumn) {
        const nodes = this.tree.get_nodes();

        if (nodes.length === this.dataFrame.rowCount) {
          this.dataFrame.columns.addNewInt('id').init(i => {
            const node = nodes[i];

            if (node.name === this.nodeNameColumn.get(i) && (!node.parent ||
              node.parent.name === this.parentNameColumn.get(i))) {
              this.nodes.set(node.id, node);
              return node.id;
            }

            return null;
          });
        } else {
          console.log('Failed to add `id` column due to node count mismatch:',
          this.dataFrame.rowCount, 'rows and', nodes.length, 'nodes');
        }

        this.nodeIdColumn = this.dataFrame.col('id');
      }
    }

    if (!this.nodeIdColumn) return;

    const nodeId = this.nodeIdColumn.get(this.dataFrame.currentRow.idx);
    if (!nodeId) return;

    this.tree.modify_selection(() => false);
    if (this.selection === 'none') return;

    const node = this.nodes.get(nodeId);
    if (!node || node.name === 'root' || node.depth === 0) return;

    const selection = (this.selection === 'path to root') ?
    this.tree.path_to_root(node) : (this.selection === 'descendants') ?
    this.tree.select_all_descendants(node, true, true) : (this.selection === 'path to root & descendants') ?
    this.tree.path_to_root(node).concat(this.tree.select_all_descendants(node, true, true)) : [];

    this.tree.modify_selection(selection);
  }
}
