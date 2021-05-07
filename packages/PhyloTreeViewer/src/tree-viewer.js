import { waitForElm } from './utils.js';

export class PhyloTreeViewer extends DG.JsViewer {
  constructor() {
    super();
    this.radialLayout = this.bool('radialLayout', false);
    this.fontSize = this.string('fontSize', '9px');
    this.selection = this.string('selection', 'path to root & descendants', { choices: [
      'none', 'path to root', 'descendants', 'path to root & descendants',
      'incident branch', 'internal branches', 'terminal branches',
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

    this.subs.push(DG.debounce(this.dataFrame.onCurrentRowChanged, 50).subscribe(() => this.render(false)));
    this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render(false)));
    this.render();
  }

  onPropertyChanged() { this.render(false); }

  _createNodeMap(nodes) {
    this.nodes.clear();

    if (this.nodeIdColumn) {
      for (let i = 0; i < this.dataFrame.rowCount; i++) {
        const node = nodes[i];
        this.nodes.set(node.id, node);
      }
    } else {
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

  _updateSelection() {
    if (this.selection === 'none') return;

    const nodeId = this.nodeIdColumn.get(this.dataFrame.currentRow.idx);
    if (!nodeId) return;

    const node = this.nodes.get(nodeId);
    if (!node || node.name === 'root' || node.depth === 0) return;

    const selection = (this.selection === 'path to root') ?
    this.tree.path_to_root(node) : (this.selection === 'descendants') ?
    this.tree.select_all_descendants(node, true, true) : (this.selection === 'path to root & descendants') ?
    this.tree.path_to_root(node).concat(this.tree.select_all_descendants(node, true, true)) :
    (this.selection === 'incident branch') ? [node] : (this.selection === 'internal branches') ?
    this.tree.select_all_descendants(node, false, true) : (this.selection === 'terminal branches') ?
    this.tree.select_all_descendants(node, true, false) : [];

    this.tree.modify_selection(selection);
  }

  _centerLayout(width, height) {
    const container = d3.select('g.phylotree-container');
    const g = container.append('g').attr('class', 'phylotree-layout');
    const selection = container.selectAll('path.branch, g.internal-node, g.node');
    const data = selection.data();
    selection.each(function() { g.append(() => this); });
    g.selectAll('path.branch, g.internal-node, g.node').data(data);

    const bbox = document.querySelector('.phylotree-layout').getBBox();
    const x = (width - bbox.width) / 2 - Math.abs(bbox.x);
    const y = (height - bbox.height) / 2 - Math.abs(bbox.y);
    g.attr('transform', `translate(${x}, ${y})`);
  }

  render(computeData = true) {
    $(this.root).empty();

    if (this.newick == null) {
      this.root.appendChild(ui.divText('Newick tag not found.', 'd4-viewer-error'));
      return;
    }

    const svg = d3.select(this.root).append("svg");
    const width = this.root.parentElement.clientWidth || this.defaultSize;
    const height = this.root.parentElement.clientHeight || this.defaultSize;
    const margin = 20;

    this.tree
      .svg(svg)
      .options({
        'label-nodes-with-name': true,
        'left-right-spacing': this.radialLayout ? 'fixed-step' : 'fit-to-size',
        'top-bottom-spacing': this.radialLayout ? 'fixed-step' : 'fit-to-size',
        'is-radial': this.radialLayout,
        'max-radius': Math.min(width, height) / 2 - margin,
        zoom: true,
      })
      .size([height, width])
      .font_size(parseInt(this.fontSize));

    this.tree.modify_selection(() => false);
    this.tree(this.parsedNewick).layout();
    // if (computeData) this.tree(this.parsedNewick);
    // this.tree.layout();
    
    if (computeData) this._createNodeMap(this.tree.get_nodes());
    if (this.nodeIdColumn) this._updateSelection();

    waitForElm(`node-${this.nodeNameColumn.get(this.dataFrame.rowCount - 1)}`)
    .then(() => {
      d3.select(this.root).selectAll('g.internal-node')
        .on('mouseover', d => {
          ui.tooltip.show(
            ui.span([d.name + (d.name ? `, ` : '') + `children: ${d.children.length}`]),
            d3.event.x + this.tooltipOffset,
            d3.event.y + this.tooltipOffset);
        })
        .on('mouseout', () => ui.tooltip.hide());

      d3.select(this.root).selectAll('g.node > text')
        .on('mouseover', d => {
          ui.tooltip.show(
            ui.span([`${d.name}, parent: ${d.parent.name}`]),
            d3.event.x + this.tooltipOffset,
            d3.event.y + this.tooltipOffset);
        })
        .on('mouseout', () => ui.tooltip.hide());

      d3.select(this.root).selectAll('path.branch')
      .on('click', d => {
        if (!d.selected) {
          for (let i = 0; i < this.dataFrame.rowCount; i++) {
            if (this.nodeIdColumn.get(i) === d.target.id) {
              this.dataFrame.currentRow = i;
              return;
            }
          }
        }
      });

      // Layout fix
      if (this.radialLayout) this._centerLayout(width, height);
    })
    .catch(e => console.log(e));
  }
}
