import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


import {_package} from '../package';
import $ from 'cash-dom';
import {phylotree as Phylotree, PhylotreeNode, TreeRender} from 'phylotree';
import * as d3 from 'd3';
import {default as newickParser} from 'phylotree/src/formats/newick';
import {NodeType, TreeDefaultPalette} from '@datagrok-libraries/bio/src/trees';
import {intToHtml} from '@datagrok-libraries/utils/src/color';

type ParsedNewickType = { error: string, json: any, };

enum PROPS {
  // -- Layout --
  leftRightSpacing = 'leftRightSpacing',
  topBottomSpacing = 'topBottomSpacing',
  radialLayout = 'radialLayout',
  margin = 'margin',
  showScale = 'showScale',
  scaleBarFontSize = 'scaleBarFontSize',
  showLabels = 'showLabels',
  alignTips = 'alignTips',
  drawSizeBubbles = 'drawSizeBubbles',

  // -- Style --
  strokeColor = 'strokeColor',
  fillColor = 'fillColor',
  strokeWidth = 'strokeWidth',
  nodeSize = 'nodeSize',
  brush = 'brush',
  fontSize = 'fontSize',
  labelNodesWithName = 'labelNodesWithName',

  // -- Behaviour --
  selection = 'selection',
  zoom = 'zoom',
  bootstrap = 'bootstrap',
  collapsible = 'collapsible',
  selectable = 'selectable',
}

enum PROPS_CATS {
  LAYOUT = 'Layout',
  STYLE = 'Style',
  BEHAVIOUR = 'Behaviour',
}

enum SelectionType {
  None = 'none',
  PathToRoot = 'path to root',
  Descendants = 'descendants',
  PathToRootAndDescendants = 'path to root & descendants',
  IncidentBranch = 'incident branch',
  InternalBranches = 'internal branches',
  TerminalBranches = 'terminal branches'
}

enum SpacingType {
  FixedStep = 'fixed-step',
  FitToSize = 'fit-to-size',
}

/**
 * https://www.npmjs.com/package/phylotree
 * https://github.com/veg/phylotree.js
 * http://phylotree.hyphy.org/
 */
export class PhyloTreeViewer extends DG.JsViewer {
  // -- Layout properties --
  [PROPS.leftRightSpacing]: string;
  [PROPS.topBottomSpacing]: string;
  [PROPS.radialLayout]: boolean;
  [PROPS.showScale]: boolean;
  [PROPS.scaleBarFontSize]: number;
  [PROPS.showLabels]: boolean;
  [PROPS.alignTips]: boolean;
  [PROPS.margin]: number;
  [PROPS.fontSize]: string;

  // -- Style --
  [PROPS.strokeColor]: number;
  [PROPS.fillColor]: number;
  [PROPS.strokeWidth]: number;
  [PROPS.nodeSize]: number;
  [PROPS.brush]: boolean;
  [PROPS.drawSizeBubbles]: boolean;

  // -- Behaviour properties --
  [PROPS.selection]: string;
  // binarySelectable: boolean;
  [PROPS.labelNodesWithName]: boolean;
  [PROPS.zoom]: boolean;
  [PROPS.bootstrap]: boolean;
  [PROPS.collapsible]: boolean;
  [PROPS.selectable]: boolean;

  tooltipOffset: number;
  defaultSize: number;
  tree: Phylotree;

  /** Map nodes by node: {@link PhylotreeNode } id */
  nodes: Map<string, PhylotreeNode>;

  newick: string | null = null;
  parsedNewick: ParsedNewickType;

  nodeIdColumn: DG.Column | null;
  nodeNameColumn: DG.Column;
  parentNameColumn: DG.Column;

  newickCol: DG.Column | null;


  constructor() {
    super();

    // -- Layout --
    this.leftRightSpacing = this.string(PROPS.leftRightSpacing, SpacingType.FixedStep,
      {category: PROPS_CATS.LAYOUT, choices: Object.values(SpacingType)});
    this.topBottomSpacing = this.string(PROPS.topBottomSpacing, SpacingType.FixedStep,
      {category: PROPS_CATS.LAYOUT, choices: Object.values(SpacingType)});
    this.radialLayout = this.bool(PROPS.radialLayout, false,
      {category: PROPS_CATS.LAYOUT});
    this.margin = this.float(PROPS.margin, 20,
      {category: PROPS_CATS.LAYOUT});
    this.showScale = this.bool(PROPS.showScale, true,
      {category: PROPS_CATS.LAYOUT});
    this.scaleBarFontSize = this.float(PROPS.scaleBarFontSize, 12,
      {category: PROPS_CATS.LAYOUT});
    this.showLabels = this.bool(PROPS.showLabels, false,
      {category: PROPS_CATS.LAYOUT});
    this.alignTips = this.bool(PROPS.alignTips, false,
      {category: PROPS_CATS.LAYOUT});
    this.drawSizeBubbles = this.bool(PROPS.drawSizeBubbles, true,
      {category: PROPS_CATS.LAYOUT});

    // -- Style --
    this.strokeColor = this.int(PROPS.strokeColor, TreeDefaultPalette.Main);
    this.fillColor = this.int(PROPS.fillColor, TreeDefaultPalette.Main);
    this.strokeWidth = this.float(PROPS.strokeWidth, window.devicePixelRatio,
      {category: PROPS_CATS.STYLE, editor: 'slider', min: 0.0, max: 16, step: 0.1});
    this.nodeSize = this.float(PROPS.nodeSize, window.devicePixelRatio * 5,
      {category: PROPS_CATS.STYLE, editor: 'slider', min: 0.0, max: 16, step: 0.1});
    this.brush = this.bool(PROPS.brush, true,
      {category: PROPS_CATS.STYLE});
    this.fontSize = this.string(PROPS.fontSize, '9px',
      {category: PROPS_CATS.STYLE});
    this.labelNodesWithName = this.bool(PROPS.labelNodesWithName, true,
      {category: PROPS_CATS.STYLE});

    // -- Behaviour --
    this.selection = this.string(PROPS.selection, SelectionType.PathToRootAndDescendants,
      {category: PROPS_CATS.BEHAVIOUR, choices: Object.values(SelectionType)});
    this.zoom = this.bool(PROPS.zoom, false,
      {category: PROPS_CATS.BEHAVIOUR});
    this.bootstrap = this.bool(PROPS.bootstrap, false,
      {category: PROPS_CATS.BEHAVIOUR});
    this.collapsible = this.bool(PROPS.collapsible, true,
      {category: PROPS_CATS.BEHAVIOUR});
    this.selectable = this.bool(PROPS.selectable, true,
      {category: PROPS_CATS.BEHAVIOUR});
    // this.binarySelectable = this.bool('binarySelectable', true,
    //   {category: PROPS_CATS.BEHAVIOUR});

    this.tooltipOffset = 10;
    this.defaultSize = 400;

    // this.root.style = 'position: absolute; left: 0; right: 0; top: 0; bottom: 0;';
    this.root.style.position = 'absolute';
    this.root.style.left = '0';
    this.root.style.right = '0';
    this.root.style.top = '0';
    this.root.style.bottom = '0';

    //eslint-disable-next-line new-cap
    this.tree = new Phylotree(':1;', {type: 'nwk'});
    this.nodes = new Map();
  }

  onTableAttached() {
    this.newick = this.dataFrame.getTag('.newick');
    if (this.newick) {
      const parsedNewick2: NodeType = newickParser(this.newick).json;

      this.parsedNewick = JSON.parse(this.dataFrame.getTag('.newickJson')!);
      phylotreeNormalize(this.parsedNewick.json);
      // With changes in newickToDf, the column 'name' is for the node identifier
      this.nodeIdColumn = this.dataFrame.col('id'); // onTableAttached
      this.nodeNameColumn = this.dataFrame.getCol('node');
      this.parentNameColumn = this.dataFrame.getCol('parent');
      this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => {
        this.render(false); // onTableAttached onSizeChanged
      }));
    } else {
      this.newickCol = this.dataFrame.columns.bySemType('newick');
    }
    this.subs.push(DG.debounce(this.dataFrame.onCurrentRowChanged, 50).subscribe((_) => {
      this.render(false); // onTableAttached onCurrentRowChanged
    }));
    this.render(); // onTableAttached
  }

  onPropertyChanged(property: DG.Property | null) {
    this.render(false);
  }

  _createNodeMap(nodeList: PhylotreeNode[]) {
    this.nodes.clear();

    if (!this.nodeIdColumn) {
      if (nodeList.length === this.dataFrame.rowCount) {
        this.nodeIdColumn = this.dataFrame.columns.addNewInt('id').init((i) => {
          const node = nodeList[i];

          if (node.name === this.nodeNameColumn.get(i) && (!node.parent ||
            node.parent.name === this.parentNameColumn.get(i))) {
            this.nodes.set(node.id, node);
            return node.id;
          }

          return null;
        });
      } else {
        _package.logger.warning('PhyloTreeViewer: Failed to add `id` column due to node count mismatch: ' +
          `${this.dataFrame.rowCount} rows and ${nodeList.length} nodes`);
      }
    }

    if (this.nodeIdColumn) {
      for (let i = 0; i < this.dataFrame.rowCount; i++) {
        const node = nodeList[i];
        this.nodes.set(node.id, node);
      }
    }
  }

  _updateSelection() {
    if (this.selection === 'none') return;

    const nodeId = this.nodeIdColumn!.get(this.dataFrame.currentRow.idx);
    if (!nodeId) return;

    const node = this.nodes.get(nodeId);
    if (!node || node.name === 'root' || node.depth === 0) return;

    /* eslint-disable max-len*/
    const selection = (this.selection === 'path to root') ? this.tree.path_to_root(node) :
      (this.selection === 'descendants') ? this.tree.select_all_descendants(node, true, true) :
        (this.selection === 'path to root & descendants') ? this.tree.path_to_root(node).concat(this.tree.select_all_descendants(node, true, true)) :
          (this.selection === 'incident branch') ? [node] :
            (this.selection === 'internal branches') ? this.tree.select_all_descendants(node, false, true) :
              (this.selection === 'terminal branches') ? this.tree.select_all_descendants(node, true, false) : [];
    /* eslint-enable max-len */

    this.tree.modify_selection(selection);
  }

  _centerLayout(width: number, height: number) {
    const container = d3.select('g.phylotree-container');
    const g = container.append('g').attr('class', 'phylotree-layout');
    const selection = container.selectAll('path.branch, g.internal-node, g.node');
    const data = selection.data();
    selection.each(function() {
      // eslint-disable-next-line
      g.append(() => this); // this is some specific object, not the viewer
    });
    g.selectAll('path.branch, g.internal-node, g.node').data(data);

    //@ts-ignore
    const bbox = document.querySelector('.phylotree-layout').getBBox();
    const x = (width - bbox.width) / 2 - Math.abs(bbox.x);
    const y = (height - bbox.height) / 2 - Math.abs(bbox.y);
    g.attr('transform', `translate(${x}, ${y})`);
  }

  render(computeData = true) {
    $(this.root).empty();

    if (this.newick == null && this.newickCol == null) {
      this.root.appendChild(ui.divText('Newick not found.', 'd4-viewer-error'));
      return;
    }

    // const svg = d3.select(this.root).append('svg');
    // const width = this.root.parentElement.clientWidth || this.defaultSize;
    // const height = this.root.parentElement.clientHeight || this.defaultSize;
    // const margin = 20;

    // this.tree
    //   .svg(svg)
    //   .options({
    //     'label-nodes-with-name': true,
    //     'left-right-spacing': this.radialLayout ? 'fixed-step' : 'fit-to-size',
    //     'top-bottom-spacing': this.radialLayout ? 'fixed-step' : 'fit-to-size',
    //     'is-radial': this.radialLayout,
    //     'max-radius': Math.min(width, height) / 2 - margin,
    //     zoom: true,
    //   })
    //   .size([height, width])
    //   .font_size(parseInt(this.fontSize));
    if (computeData) {
      //this.tree.update(d3.hierarchy(this.parsedNewick.json));
      this.tree = new Phylotree('', {type: () => this.parsedNewick});
      const treeRoot = this.tree.getRootNode();
      const nodeList = getNodeList(treeRoot);
      this._createNodeMap(nodeList);
    }
    const cw: number = this.root.clientWidth;
    const ch: number = this.root.clientHeight;
    const display: TreeRender = this.tree.render({
      'is-radial': this.radialLayout,
      'left-right-spacing': this.radialLayout ? 'fixed-step' : 'fit-to-size',
      'top-bottom-spacing': this.radialLayout ? 'fixed-step' : 'fit-to-size',
      'label-nodes-with-name': this.labelNodesWithName,
      'max-radius': Math.min(cw, ch) / 2 - this.margin,
      'height': Math.max(ch, 0),
      'width': Math.max(cw - this.margin, 0),
      'layout': 'left-to-right',
      'branches': 'step',
      'scaling': true,
      'bootstrap': this.bootstrap,
      'color-fill': true,
      'internal-names': false,
      'selectable': this.selectable,
      'restricted-selectable': false,
      'collapsible': this.collapsible,
      'left-offset': this.margin,
      'show-scale': this.showScale,
      'scale-bar-font-size': this.scaleBarFontSize,
      'draw-size-bubbles': this.drawSizeBubbles,
      //TypeError: Cannot read properties of undefined (reading 'options')
      'binary-selectable': false,
      'annular-limit': 0.38196601125010515,
      'compression': 0.2,
      'align-tips': this.alignTips,
      'maximum-per-node-spacing': 100,
      'minimum-per-node-spacing': 2,
      'maximum-per-level-spacing': 100,
      'minimum-per-level-spacing': 10,
      'transitions': null,
      'brush': this.brush,
      'reroot': true,
      'hide': true,
      'zoom': this.zoom,
      'show-menu': true,
      'show-labels': this.showLabels,
      'node-span': null
    });


    // let treeEl: SVGSVGElement;
    // if (this.radialLayout){
    //
    // } else {
    //   treeEl = display
    // }
    display.modifySelection(() => false);
    const treeEl: SVGSVGElement = display.show();
    const strokeColorValue = intToHtml(this.strokeColor);
    const fillColorValue = intToHtml(this.fillColor);
    $(treeEl).find('g.phylotree-container path.branch').get()
      .forEach((pathEl) => {
        pathEl.style.strokeWidth = `${this.strokeWidth}px`;
        pathEl.style.stroke = strokeColorValue;
      });
    $(treeEl).find('g.phylotree-container g.node circle').get()
      .forEach((circleEl) => {
        // @ts-ignore
        circleEl.style.r = `${this.nodeSize}px`;
        circleEl.style.strokeWidth = `${this.strokeWidth}px`;
        circleEl.style.stroke = strokeColorValue;
        circleEl.style.fill = fillColorValue;
      });
    $(treeEl).find('g.phylotree-container g.internal-node circle').get()
      .forEach((circleEl) => {
        // @ts-ignore
        circleEl.style.r = this.strokeWidth / 2;
        circleEl.style.strokeWidth = '0';
        circleEl.style.fill = fillColorValue;
      });
    this.root.append(treeEl);

    const idx = this.dataFrame.currentRowIdx;
    //this.tree.update(this.newick ? this.parsedNewick : newickParser(this.newickCol.get((idx == -1) ? 1 : idx)));
    // if (computeData) this.tree(this.parsedNewick);
    // this.tree.layout();

    if (this.nodeIdColumn) this._updateSelection();

    // waitForElm(`node-${this.nodeNameColumn.get(this.dataFrame.rowCount - 1)}`)
    //   .then(() => {
    //     d3.select(this.root).selectAll('g.internal-node')
    //       .on('mouseover', d => {
    //         ui.tooltip.show(
    //           ui.span([d.name + (d.name ? `, ` : '') + `children: ${d.children.length}`]),
    //           d3.event.x + this.tooltipOffset,
    //           d3.event.y + this.tooltipOffset);
    //       })
    //       .on('mouseout', () => ui.tooltip.hide());
    //
    //     d3.select(this.root).selectAll('g.node > text')
    //       .on('mouseover', d => {
    //         ui.tooltip.show(
    //           ui.span([`${d.name}, parent: ${d.parent.name}`]),
    //           d3.event.x + this.tooltipOffset,
    //           d3.event.y + this.tooltipOffset);
    //       })
    //       .on('mouseout', () => ui.tooltip.hide());
    //
    //     d3.select(this.root).selectAll('path.branch')
    //       .on('click', d => {
    //         if (!d.selected) {
    //           for (let i = 0; i < this.dataFrame.rowCount; i++) {
    //             if (this.nodeIdColumn.get(i) === d.target.id) {
    //               this.dataFrame.currentRowIdx = i;
    //               return;
    //             }
    //           }
    //         }
    //       });
    //
    //     // Layout fix
    //     if (this.radialLayout) this._centerLayout(width, height);
    //   })
    //   .catch(e => console.log(e));
  }
}

function isLeaf<TNode extends NodeType>(node: TNode): boolean {
  return !node.children || node.children.length == 0;
}

function getNodeList<TNode extends NodeType>(node: TNode | null): TNode[] {
  if (!node) return [];

  if (isLeaf(node)) {
    return [node]; // node is a leaf
  } else {
    const childNodeListList: TNode[][] = node.children!
      .map((child) => { return getNodeList(child as TNode); });
    return ([] as TNode[]).concat(
      [node],
      ...childNodeListList);
  }
}

/** phylotree requires node attribute children */
function phylotreeNormalize<TNode extends NodeType>(node: TNode | null): void {
  if (node) {
    if (node.children == undefined) node.children = [];
    for (const childNode of node.children) phylotreeNormalize(childNode);
  }
}
