import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

import * as rxjs from 'rxjs';
import {TREE_TAGS} from '../consts';
import {Observable, Subject, Unsubscribable} from 'rxjs';
import $ from 'cash-dom';
import {PickingInfo} from '@deck.gl/core/typed';
import {MjolnirPointerEvent} from 'mjolnir.js';

// TODO: add test for these properties existing.
/** Represents a single tree node */
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

const TreeTypesTypes: { [streeType: string]: string } = {
  [bio.TreeTypesNames.Radial]: bio.TreeTypes.Radial,
  [bio.TreeTypesNames.Rectangular]: bio.TreeTypes.Rectangular, // rectangular edges, leaves listed vertically
  [bio.TreeTypesNames.Polar]: bio.TreeTypes.Circular,
  [bio.TreeTypesNames.Diagonal]: bio.TreeTypes.Diagonal,
  [bio.TreeTypesNames.Orthogonal]: bio.TreeTypes.Hierarchical, // rectangular edges, leaves listed horizontally
};


export class PhylocanvasGlViewer extends DG.JsViewer implements bio.IPhylocanvasGlViewer {
  private viewed: boolean = false;

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

  parsedNewick: Object;

  nodeIdColumn: DG.Column;
  nodeNameColumn: DG.Column;
  parentNameColumn: DG.Column;

//   title: string;
//
  treeDiv: HTMLDivElement | null;
  viewer: bio.PhylocanvasGL | null;
  /** Container to store prop values while phylocanvasViewer is not created yet */
  phylocanvasProps: { [p: string]: any } = {};

  private _onAfterRender = new Subject<{ gl: WebGLRenderingContext }>();

  private _onHover = new Subject<{ info: PickingInfo, event: MjolnirPointerEvent }>();

  private viewSubs: Unsubscribable[] = [];

  get onAfterRender(): Observable<{ gl: WebGLRenderingContext }> { return this._onAfterRender; }

  get onHover(): Observable<{ info: PickingInfo, event: MjolnirPointerEvent }> { return this._onHover; }

  constructor() {
    super();

    this.alignLabels = this.bool('alignLabels', false);

    this.treeType = this.string('treeType', bio.TreeTypesNames.Rectangular,
      {category: PROPS_CATS.APPEARANCE, choices: Object.values(bio.TreeTypesNames)});
    this.lineWidth = this.float('lineWidth', 1, {category: PROPS_CATS.APPEARANCE,});
    this.nodeSize = this.float('nodeSize', 14,
      {category: PROPS_CATS.APPEARANCE, editor: 'slider', min: 1, max: 32,});
    this.strokeColor = this.int('strokeColor', 0x222222, {category: PROPS_CATS.APPEARANCE});
    this.showShapes = this.bool('showShapes', true, {category: PROPS_CATS.APPEARANCE});
    this.nodeShape = this.string('nodeShape', 'Circle',
      {category: PROPS_CATS.APPEARANCE, choices: Object.keys(bio.Shapes)});
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

  // It breaks creating viewer with DataFrame.plot.fromType(...)
  // public override get dataFrame(): DG.DataFrame {
  //   throw new Error('Not supported while onTableAttached() problem');
  // }
  //
  // public override set dataFrame(value: DG.DataFrame) {
  //   throw new Error('Not supported while onTableAttached() problem');
  // }

  private _nwkDf: DG.DataFrame;

  public get nwkDf(): DG.DataFrame {
    return this._nwkDf;
  }

  public set nwkDf(value: DG.DataFrame) {
    console.debug('PhyloTreeViewer: PhylocanvasGlViewer.onTableAttached() ' +
      `this.dataFrame = ${!this.nwkDf ? 'null' : 'value'} )`);

    if (this.viewed) {
      this.destroyView();
      this.viewed = false;
    }

    this._nwkDf = value;

    this.newick = this.nwkDf.getTag(TREE_TAGS.NEWICK);
    if (this.newick) {
      this.parsedNewick = JSON.parse(this.nwkDf.getTag('.newickJson')!);
      //this.nodeIdColumn = this.dataFrame.getCol('id');
      this.nodeNameColumn = this.nwkDf.getCol('node');
      this.parentNameColumn = this.nwkDf.getCol('parent');
    }
    // else {
    //   this.newickCol = this.dataFrame.columns.bySemType('newick');
    // }

    if (!this.viewed) {
      this.buildView();
      this.viewed = true;
    }
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
    if (!this.viewer)
      return;

    // this.treeDiv.innerText = `${this.root.clientWidth} x ${this.root.clientHeight}`;
    const cw: number = this.root.clientWidth;
    const ch: number = this.root.clientHeight;
    console.debug(`PhyloTreeViewer: PhylocanvasGlViewer.calcSize( ${cw.toString()} x ${ch.toString()} )`);

    const width = this.root.clientWidth;
    const height = this.root.clientHeight;
    if (this.treeDiv) {
      this.treeDiv.style.width = `${width}px`;
      this.treeDiv.style.height = `${height}px`;

      this.viewer.setProps({size: this.treeDiv.getBoundingClientRect()});
    }
  }

  private destroyView() {
    console.debug('PhyloTreeViewer: PhylocanvasGlViewer.destroyView() ');

    this.viewSubs.forEach((s: Unsubscribable) => { s.unsubscribe(); });
  }

  private buildView() {
    console.debug('PhyloTreeViewer: PhylocanvasGlViewer.buildView() ');

    //const color: string = `#bbff${Math.ceil(128 + Math.random() * 127).toString(16)}`;
    if (!this.treeDiv && !this.viewer) {
      this.treeDiv = ui.div([], {
        style: {
          //backgroundColor: color,
          width: '100px',
          height: '100px',
        }
      });
      this.root.appendChild(this.treeDiv);

      // default props required to prevent throwing exception
      const defaultNwkStr = '(NONE:1);'; // minimal tree required to not throw exception in PhylocanvasGL
      const defaultNwkRoot = bio.Newick.parse_newick(defaultNwkStr);
      const props: { [p: string]: any } = {
        source: {type: 'biojs', data: defaultNwkRoot},
      };

      Object.assign(props, this.phylocanvasProps);
      this.viewer = new bio.PhylocanvasGL(this.treeDiv, props);
      this.viewer.deck.setProps({
        useDevicePixels: true,
        onAfterRender: this.viewerDeckOnAfterRender.bind(this),
        onHover: this.viewerOnHover.bind(this),
      });
      this.viewer.view.style.backgroundImage = 'none';

      this.calcSize();
      let k = 11;

      // this.subs.push(
      //   rxjs.fromEvent<MouseEvent>(this.phylocanvasViewer.deck.canvas!, 'mousedown').subscribe((args: MouseEvent) => {
      //
      //   }));
    }
    const newickRoot = bio.Newick.parse_newick(this.newick!);
    this.viewer!.setProps({source: {type: 'biojs', data: newickRoot}});

    // this.phyloTreeViewer.selectNode = this.tvSelectNode.bind(this);
    //@ts-ignore
    //this.viewer.handleHover = this.viewerHandleHover.bind(this);

    // this._viewSubs.push(DG.debounce(this.dataFrame.onCurrentRowChanged, 50).subscribe((_) => this.render(false)));
  }

  override onTableAttached() {
    super.onTableAttached();
    this.nwkDf = this.dataFrame;
  }

  override detach() {
    if (this.viewed) {
      this.destroyView();
      this.viewed = false;
    }

    if (this.viewer) {
      this.viewer.destroy();
      this.viewer = null;
      // this.phylocanvasProps = {}; // the value is required to restore view in destroyView/buildView
    }

    if (this.treeDiv) {
      this.treeDiv.remove();
      this.treeDiv = null;
    }

    super.detach();
  }

  override onPropertyChanged(property: DG.Property | null) {
    super.onPropertyChanged(property);

    if (!property) {
      console.warn('PhyloTreeViewer: PhylocanvasGlViewer.onPropertyChanged() No property value');
      return;
    }

    const setProps = (props: { [p: string]: any }) => {
      Object.assign(this.phylocanvasProps, props);
      if (this.viewer)
        this.viewer.setProps(props);
    };

    switch (property.name) {
    case 'treeType':
      const treeTypeValue = TreeTypesTypes[this.treeType];
      setProps({type: treeTypeValue});
      break;

    case 'nodeShape':
      const nodeShapeValue: string = bio.Shapes[this.nodeShape];
      setProps({nodeShape: nodeShapeValue});
      break;

    case 'strokeColor':
      const strokeColorValue: string = `#${(this.strokeColor & 0xFFFFFF).toString(16).padStart(6, '0')}`;
      setProps({strokeColour: strokeColorValue});
      console.debug('PhyloTreeViewer: PhylocanvasGlViewer.onPropertyChanged() ' +
        `${property.name} = ${strokeColorValue} .`);
      break;

    case 'fillColor':
      const fillColorValue: string = `#${(this.fillColor & 0x00FFFFFF).toString(16).padStart(6, '0')}`;
      setProps({fillColour: fillColorValue});
      console.debug('PhyloTreeViewer: PhylocanvasGlViewer.onPropertyChanged() ' +
        `${property.name} = ${fillColorValue} .`);
      break;

    case 'highlightColor':
      const highlightColorValue = `#${(this.highlightColor & 0xFFFFFF).toString(16).padStart(6, '0')}`;
      setProps({highlightColour: highlightColorValue});
      console.debug('PhyloTreeViewer: PhylocanvasGlViewer.onPropertyChanged() ' +
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
      setProps({[property.name]: propValue});
      console.debug('PhyloTreeViewer: PhylocanvasGlViewer.onPropertyChanged() ' +
        `${property.name} = ${propValue.toString()} .`);
      break;

    default:
      console.warn('PhyloTreeViewer: PhylocanvasGlViewer.onPropertyChanged() ' +
        `Unexpected property '${property.name}'.`);
    }
  }

  setProps(updater: { [propName: string]: any }): void {
    Object.assign(this.phylocanvasProps, updater);
    if (this.viewer !== null)
      this.viewer.setProps(updater);
  }

  protected viewerDeckOnAfterRender({gl}: { gl: WebGLRenderingContext }): void {
    this._onAfterRender.next({gl});
  }

  protected viewerOnHover(info: PickingInfo, event: MjolnirPointerEvent): void {
    this._onHover.next({info, event});
  }
}

interface TreeMap {
  index: number[];
  items: string[];
}

type TreeLeaves = { [key: string]: TreeMap };
type TreeLeavesMap = { [key: string]: number[] };
