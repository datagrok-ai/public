import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PhylocanvasGL, TreeTypes, Shapes} from '@phylocanvas/phylocanvas.gl';
import {TREE_TAGS} from '../consts';
import {Unsubscribable} from 'rxjs';
// // import {TreeAnalyzer, PhylocanvasTreeNode, getVId} from './utils/tree-stats';
// import {Subscription, Unsubscribable} from 'rxjs';

// TODO: add test for these properties existing.
/**
 * Represents a single tree node.
 */
export interface PhylocanvasTreeNode {
  branchLength: 0;
  children: PhylocanvasTreeNode[];
  id: string;
  isCollapsed: boolean;
  isHidden: boolean;
  isLeaf: boolean;
  name: string;
  postIndex: number;
  preIndex: number;
  totalLeaves: number;
  totalNodes: number;
  totalSubtreeLength: number;
  visibleLeaves: number;
}

enum PROPS_CATS {
  APPEARANCE = 'Appearance',
  BEHAVIOR = 'Behavior',
  LAYOUT = 'Layout',
}

export enum TreeTypesNames {
  Radial = 'Radial',
  /** Rectangular edges, leaves listed __vertically__ */
  Rectangular = 'Rectangular',
  Polar = 'Polar',
  Diagonal = 'Diagonal',
  /** Rectangular edges, leaves listed __horizontally__ */
  Orthogonal = 'Orthogonal',
}

const TreeTypesTypes: { [streeType: string]: string } = {
  [TreeTypesNames.Radial]: TreeTypes.Radial,
  [TreeTypesNames.Rectangular]: TreeTypes.Rectangular, // rectangular edges, leaves listed vertically
  [TreeTypesNames.Polar]: TreeTypes.Circular,
  [TreeTypesNames.Diagonal]: TreeTypes.Diagonal,
  [TreeTypesNames.Orthogonal]: TreeTypes.Hierarchical, // rectangular edges, leaves listed horizontally
};

type StylesType = { [nodeName: string]: { [prop: string]: any } };

export class PhylocanvasGlViewer extends DG.JsViewer {
  alignLabels: boolean;

  treeType: string;
  lineWidth: number;
  nodeSize: number;
  strokeColor: number;
  showShapes: boolean;
  nodeShape: string;
  haloRadius: number;
  haloWidth: number;

  branchZoom: number;
  fillColor: number;
  showShapeBorders: boolean;
  shapeBorderWidth: number;
  shapeBorderAlpha: number;
  highlightColor: number;
  showBranchLengths: boolean;
  showEdges: boolean;
  showLabels: boolean;
  showLeafLabels: boolean;
  showInternalLabels: boolean;
  fontFamily: string;
  fontSize: number;

  padding: number;
  treeToCanvasRatio: number;

  interactive: boolean;
  stepZoom: number;
  zoom: number;

  newick: string | null = null;
  styles: StylesType | null = null;

  parsedNewick: Object;

  nodeIdColumn: DG.Column;
  nodeNameColumn: DG.Column;
  parentNameColumn: DG.Column;

//   title: string;
//
  treeDiv: HTMLDivElement | null;
  phylocanvasViewer: PhylocanvasGL | null;

  private viewSubs: Unsubscribable[] = [];

  constructor() {
    super();

    this.alignLabels = this.bool('alignLabels', false);

    this.treeType = this.string('treeType', TreeTypesNames.Rectangular,
      {category: PROPS_CATS.APPEARANCE, choices: Object.values(TreeTypesNames)});
    this.lineWidth = this.float('lineWidth', 1, {category: PROPS_CATS.APPEARANCE,});
    this.nodeSize = this.float('nodeSize', 14,
      {category: PROPS_CATS.APPEARANCE, editor: 'slider', min: 1, max: 32,});
    this.strokeColor = this.int('strokeColor', 0x222222, {category: PROPS_CATS.APPEARANCE});
    this.showShapes = this.bool('showShapes', true, {category: PROPS_CATS.APPEARANCE});
    this.nodeShape = this.string('nodeShape', 'Circle',
      {category: PROPS_CATS.APPEARANCE, choices: Object.keys(Shapes)});
    this.haloRadius = this.float('haloRadius', 12,
      {category: PROPS_CATS.APPEARANCE, editor: 'slider', min: 1, max: 32});
    this.haloWidth = this.float('haloWidth', 4,
      {category: PROPS_CATS.APPEARANCE, editor: 'slider', min: 1, max: 32});
    this.fillColor = this.int('fillColor', 0x333333, {category: PROPS_CATS.APPEARANCE,});
    this.showShapeBorders = this.bool('showShapeBorders', false, {category: PROPS_CATS.APPEARANCE,});
    this.shapeBorderWidth = this.float('shapeBorderWidth', 1,
      {category: PROPS_CATS.APPEARANCE, editor: 'slider', min: 0, max: 50.});
    this.shapeBorderAlpha = this.float('shapeBorderAlpha', 0.14,
      {category: PROPS_CATS.APPEARANCE, editor: 'slider', min: 0.0, max: 1.0,});
    this.highlightColor = this.int('highlightColor', 0x3C7383, {category: PROPS_CATS.APPEARANCE});
    this.showBranchLengths = this.bool('showBranchLengths', false, {category: PROPS_CATS.APPEARANCE});
    this.showEdges = this.bool('showEdges', true, {category: PROPS_CATS.APPEARANCE});
    this.showLabels = this.bool('showLabels', false, {category: PROPS_CATS.APPEARANCE});
    this.showLeafLabels = this.bool('showLeafLabels', false, {category: PROPS_CATS.APPEARANCE});
    this.showInternalLabels = this.bool('showInternalLabels', false, {category: PROPS_CATS.APPEARANCE});
    this.fontFamily = this.string('fontFamily', 'monospace', {category: PROPS_CATS.APPEARANCE});
    this.fontSize = this.float('fontSize', 16,
      {category: PROPS_CATS.APPEARANCE, editor: 'slider', min: 2, max: 48,});

    this.padding = this.float('padding', 16,
      {category: PROPS_CATS.LAYOUT, editor: 'slider', min: 0, max: 50,});
    this.treeToCanvasRatio = this.float('treeToCanvasRatio', 0.25,
      {category: PROPS_CATS.BEHAVIOR, editor: 'slider', min: -4, max: 4, step: 0.25});

    this.interactive = this.bool('interactive', false, {category: PROPS_CATS.BEHAVIOR});
    this.branchZoom = this.float('branchZoom', 0,
      {category: PROPS_CATS.BEHAVIOR, editor: 'slider', min: -4, max: 4, step: 0.1});
    this.stepZoom = this.float('stepZoom', 0,
      {category: PROPS_CATS.BEHAVIOR, editor: 'slider', min: -4, max: 4, step: 0.1});
    this.zoom = this.float('zoom', 0,
      {category: PROPS_CATS.BEHAVIOR, editor: 'slider', min: -4, max: +4, step: 0.1});

    this.subs.push(
      ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));
  }

  public async init() {
    // await this.buildView();
    // .catch((ex) => {
    //   console.error(`TreeBrowser.init() > buildView() error:\n${ex ? ex.toString() : 'none'}`);
    // });
  }

//   /**
//    * Returns filtered rows in MLB grid.
//    * @param {DG.DataFrame} mlbDf
//    * @param {{}} vIdIndex
//    * @return {string[]} Rows is left after filtering.
//    */
//   private static _calcFilteredItems(mlbDf: DG.DataFrame, vIdIndex: {}): string[] {
//     let res: string [] = [];
//     if (mlbDf) {
//       const filteredIndices = new Set(mlbDf.filter.getSelectedIndexes());
//       const isntFiltered = (x: number) => filteredIndices.has(x);
//       res = Object.keys(vIdIndex).filter((_, i) => isntFiltered(i));
//     }
//     return res;
//   }

  private rootOnSizeChanged(args: any): void {
    console.debug('PhyloTreeViewer: PhylocanvasGlViewer.rootOnSizeChanged()');
    this.calcSize();
  }


  calcSize(): void {
    if (!this.phylocanvasViewer)
      return;

    // this.treeDiv.innerText = `${this.root.clientWidth} x ${this.root.clientHeight}`;
    const cw: number = this.root.clientWidth;
    const ch: number = this.root.clientHeight;
    console.debug(`PhylocanvasGlViewer: TreeBrowser.calcSize( ${cw.toString()} x ${ch.toString()} )`);

    const width = this.root.clientWidth;
    const height = this.root.clientHeight;
    if (this.treeDiv) {
      this.treeDiv.style.width = `${width}px`;
      this.treeDiv.style.height = `${height}px`;

      this.phylocanvasViewer.setProps({size: this.treeDiv.getBoundingClientRect()});
    }
  }

  private async destroyView() {
    this.viewSubs.forEach((s: Unsubscribable) => { s.unsubscribe(); });

    if (this.phylocanvasViewer) {
      this.phylocanvasViewer!.destroy();
      this.phylocanvasViewer = null;

      this.styles = null;
    }

    if (this.treeDiv) {
      this.treeDiv.remove();
      this.treeDiv = null;
    }
  }

  private async buildView() {
    //const color: string = `#bbff${Math.ceil(128 + Math.random() * 127).toString(16)}`;
    if (!this.treeDiv) {
      this.treeDiv = ui.div([], {
        style: {
          //backgroundColor: color,
          width: '100px',
          height: '100px',
        }
      });
      this.root.appendChild(this.treeDiv);
    }
    // const treeNode = mlbView.dockManager.dock(treeDiv, DG.DOCK_TYPE.DOWN);

    if (!this.phylocanvasViewer) {
      const props: { [p: string]: any } = {
        interactive: this.interactive,
        showLabels: this.showLabels,
        showLeafLabels: this.showLeafLabels,
        nodeShape: this.nodeShape,
        fontFamily: this.fontFamily,
        fontSize: this.fontSize,
        nodeSize: this.nodeSize,
        size: this.treeDiv.getBoundingClientRect(),
        source: this.newick,
        type: TreeTypesTypes[this.type],
      };
      if (this.styles)
        props['styles'] = this.styles;
      this.phylocanvasViewer = new PhylocanvasGL(this.treeDiv, props);
      this.phylocanvasViewer.deck.setProps({useDevicePixels: true});
      this.calcSize();
    }
    // this.phyloTreeViewer.selectNode = this.tvSelectNode.bind(this);
    // this.phyloTreeViewer.handleHover = this.tvHandleHover.bind(this);

    // this._viewSubs.push(DG.debounce(this.dataFrame.onCurrentRowChanged, 50).subscribe((_) => this.render(false)));
  }

  override async onTableAttached() {
    console.debug(`PhyloTreeViewer: PhylocanvasGlViewer.onTableAttached( dataFrame = ${!this.dataFrame ? 'null' : 'value'} )`);
    window.setTimeout(async () => {
      await this.destroyView();

      this.newick = this.dataFrame.getTag(TREE_TAGS.NEWICK);
      if (this.newick) {
        this.parsedNewick = JSON.parse(this.dataFrame.getTag('.newickJson')!);
        //this.nodeIdColumn = this.dataFrame.getCol('id');
        this.nodeNameColumn = this.dataFrame.getCol('node');
        this.parentNameColumn = this.dataFrame.getCol('parent');

        // this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render(false)));
      }
      // else {
      //   this.newickCol = this.dataFrame.columns.bySemType('newick');
      // }

      await this.buildView();
    }, 0 /* next event cycle */);
  }

  override onPropertyChanged(property: DG.Property | null) {
    super.onPropertyChanged(property);

    if (!property) {
      console.warn('PhyloTreeViewer: PhylocanvasGlViewer.onPropertyChanged() No property value');
      return;
    }

    if (this.phylocanvasViewer) {
      switch (property.name) {
      case 'treeType':
        const treeTypeValue = TreeTypesTypes[this.treeType];
        this.phylocanvasViewer.setProps({type: treeTypeValue});
        break;

      case 'nodeShape':
        const nodeShapeValue: string = Shapes[this.nodeShape];
        this.phylocanvasViewer.setProps({nodeShape: nodeShapeValue});
        break;

      case 'strokeColor':
        const strokeColorValue: string = `#${(this.strokeColor & 0xFFFFFF).toString(16).padStart(6, '0')}`;
        this.phylocanvasViewer.setProps({strokeColour: strokeColorValue});
        console.debug('PhyloTreeView: PhylocanvasGlViewer.onPropertyChanged() ' +
          `${property.name} = ${strokeColorValue} .`);
        break;

      case 'fillColor':
        const fillColorValue: string = `#${(this.fillColor & 0x00FFFFFF).toString(16).padStart(6, '0')}`;
        this.phylocanvasViewer.setProps({fillColour: fillColorValue});
        console.debug('PhyloTreeView: PhylocanvasGlViewer.onPropertyChanged() ' +
          `${property.name} = ${fillColorValue} .`);
        break;

      case 'highlightColor':
        const highlightColorValue = `#${(this.highlightColor & 0xFFFFFF).toString(16).padStart(6, '0')}`;
        this.phylocanvasViewer.setProps({highlightColour: highlightColorValue});
        console.debug('PhyloTreeView: PhylocanvasGlViewer.onPropertyChanged() ' +
          `${property.name} = ${highlightColorValue} .`);
        break;

      case 'lineWidth':
      case 'nodeSize':
      case 'showShapes':
      case 'haloRadius':
      case 'haloWidth':
      case 'branchZoom':
      case 'showShapeBorders':
      case 'shapeBorderWidth':
      case 'shapeBorderAlpha':
      case 'showBranchLengths':
      case 'showEdges':
      case 'showLabels':
      case 'showLeafLabels':
      case 'showInternalLabels':
      case 'fontFamily':
      case 'fontSize':
      case 'padding':
      case 'treeToCanvasRatio':
      case 'interactive':
      case 'stepZoom':
      case 'zoom':
        const propValue: any = this[property.name];
        this.phylocanvasViewer.setProps({[property.name]: propValue});
        console.debug('PhyloTreeView: PhylocanvasGlViewer.onPropertyChanged() ' +
          `${property.name} = ${propValue.toString()} .`);
        break;

      default:
        console.warn('PhyloTreeView: PhylocanvasGlViewer.onPropertyChanged() ' +
          `Unexpected property '${property.name}'.`);
      }
    }
  }

  setStyles(styles: { [nodeName: string]: { [prop: string]: any } }): void {
    this.styles = styles;
    if (this.phylocanvasViewer)
      this.phylocanvasViewer?.setProps({styles: this.styles});
  }
}

interface TreeMap {
  index: number[];
  items: string[];
}

type TreeLeaves = { [key: string]: TreeMap };
type TreeLeavesMap = { [key: string]: number[] };
