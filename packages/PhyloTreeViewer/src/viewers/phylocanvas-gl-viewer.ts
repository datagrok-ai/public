import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import * as rxjs from 'rxjs';
import $ from 'cash-dom';
import {TREE_TAGS} from '../consts';
import {Observable, Subject, Unsubscribable} from 'rxjs';
import {PickingInfo} from '@deck.gl/core/typed';
import {MjolnirPointerEvent} from 'mjolnir.js';
import {PhylocanvasGL} from '@phylocanvas/phylocanvas.gl';
import {intToHtml} from '@datagrok-libraries/utils/src/color';
import {IPhylocanvasGlViewer, TreeTypesNames} from '@datagrok-libraries/bio/src/viewers/phylocanvas-gl-viewer';
import {TreeDefaultPalette} from '@datagrok-libraries/bio/src/trees';
import {parseNewick, Shapes, TreeTypes} from '@datagrok-libraries/bio/src/trees/phylocanvas';
import {testEvent} from '@datagrok-libraries/test/src/test';

// TODO: add test for these properties existing.

const enum PROPS {
  // -- Data --
  nodeColumnName = 'nodeColumnName',
  parentColumnName = 'parentColumnName',

  // -- Style --
  alignLabels = 'alignLabels',
  treeType = 'treeType',
  strokeWidth = 'strokeWidth',
  nodeSize = 'nodeSize',
  strokeColor = 'strokeColor',
  showShapes = 'showShapes',
  nodeShape = 'nodeShape',
  haloRadius = 'haloRadius',
  haloWidth = 'haloWidth',
  fillColor = 'fillColor',
  showShapeBorders = 'showShapeBorders',
  shapeBorderWidth = 'shapeBorderWidth',
  shapeBorderAlpha = 'shapeBorderAlpha',
  highlightColor = 'highlightColor',
  showBranchLengths = 'showBranchLengths',
  showEdges = 'showEdges',
  showLabels = 'showLabels',
  showLeafLabels = 'showLeafLabels',
  showInternalLabels = 'showInternalLabels',
  fontFamily = 'fontFamily',
  fontSize = 'fontSize',
  padding = 'padding',
  treeToCanvasRatio = 'treeToCanvasRatio',
  interactive = 'interactive',
  branchZoom = 'branchZoom',
  stepZoom = 'stepZoom',
  zoom = 'zoom',
}

const enum PROPS_CATS {
  APPEARANCE = 'Appearance',
  BEHAVIOR = 'Behavior',
  LAYOUT = 'Layout',
  DATA = 'Data',
}

const TreeTypesTypes: { [treeType: string]: string } = {
  [TreeTypesNames.Radial]: TreeTypes.Radial,
  [TreeTypesNames.Rectangular]: TreeTypes.Rectangular, // rectangular edges, leaves listed vertically
  [TreeTypesNames.Polar]: TreeTypes.Circular,
  [TreeTypesNames.Diagonal]: TreeTypes.Diagonal,
  [TreeTypesNames.Orthogonal]: TreeTypes.Hierarchical, // rectangular edges, leaves listed horizontally
};

/** minimal tree required to not throw exception in PhylocanvasGL */
const defaultNwkRoot = parseNewick('(NONE:1);');

export class PhylocanvasGlViewer extends DG.JsViewer implements IPhylocanvasGlViewer {
  private viewed: boolean = false;

  [PROPS.alignLabels]: boolean;

  [PROPS.treeType]: string;
  [PROPS.strokeWidth]: number;
  [PROPS.nodeSize]: number;
  [PROPS.strokeColor]: number;
  [PROPS.showShapes]: boolean;
  [PROPS.nodeShape]: string;
  [PROPS.haloRadius]: number;
  [PROPS.haloWidth]: number;

  [PROPS.branchZoom]: number;
  [PROPS.fillColor]: number;
  [PROPS.showShapeBorders]: boolean;
  [PROPS.shapeBorderWidth]: number;
  [PROPS.shapeBorderAlpha]: number;
  [PROPS.highlightColor]: number;
  [PROPS.showBranchLengths]: boolean;
  [PROPS.showEdges]: boolean;
  [PROPS.showLabels]: boolean;
  [PROPS.showLeafLabels]: boolean;
  [PROPS.showInternalLabels]: boolean;
  [PROPS.fontFamily]: string;
  [PROPS.fontSize]: number;

  [PROPS.padding]: number;
  [PROPS.treeToCanvasRatio]: number;

  [PROPS.interactive]: boolean;
  [PROPS.stepZoom]: number;
  [PROPS.zoom]: number;

  [PROPS.nodeColumnName]: string;
  [PROPS.parentColumnName]: string;

  newick: string | null = null;

  parsedNewick: Object;

  nodeCol: DG.Column;
  parentCol: DG.Column;

  // title: string;

  treeDiv: HTMLDivElement | null;
  viewer: PhylocanvasGL | null;
  /** Container to store prop values while phylocanvasViewer is not created yet */
  phylocanvasProps: { [p: string]: any } = {
    strokeColour /* [PROPS.strokeColor] */: intToHtml(TreeDefaultPalette.Main),
    fillColour /* [PROPS.fillColor] */: intToHtml(TreeDefaultPalette.Main),
    highlightColour /* [PROPS.highlightColor] */: intToHtml(TreeDefaultPalette.Selection),
  };

  private _onAfterRender = new Subject<{ gl: WebGLRenderingContext }>();

  private _onHover = new Subject<{ info: PickingInfo, event: MjolnirPointerEvent }>();

  private viewSubs: Unsubscribable[] = [];

  get onAfterRender(): Observable<{ gl: WebGLRenderingContext }> { return this._onAfterRender; }

  get onHover(): Observable<{ info: PickingInfo, event: MjolnirPointerEvent }> { return this._onHover; }

  constructor() {
    super();

    // -- Data --
    this.nodeColumnName = this.string(PROPS.nodeColumnName, 'node',
      {category: PROPS_CATS.DATA});
    this.parentColumnName = this.string(PROPS.parentColumnName, 'parent',
      {category: PROPS_CATS.DATA});

    // -- Style --
    this.alignLabels = this.bool(PROPS.alignLabels, false,
      {category: PROPS_CATS.APPEARANCE});

    this.treeType = this.string(PROPS.treeType, TreeTypesNames.Rectangular,
      {category: PROPS_CATS.APPEARANCE, choices: Object.values(TreeTypesNames)});
    this.strokeWidth = this.float(PROPS.strokeWidth, 1,
      {category: PROPS_CATS.APPEARANCE, editor: 'slider', min: 0, max: 16, step: 0.1});
    this.nodeSize = this.float(PROPS.nodeSize, 14,
      {category: PROPS_CATS.APPEARANCE, editor: 'slider', min: 1, max: 16, step: 0.1});
    this.strokeColor = this.int(PROPS.strokeColor, TreeDefaultPalette.Main,
      {category: PROPS_CATS.APPEARANCE});
    this.showShapes = this.bool(PROPS.showShapes, true, {category: PROPS_CATS.APPEARANCE});
    this.nodeShape = this.string(PROPS.nodeShape, 'Circle',
      {category: PROPS_CATS.APPEARANCE, choices: Object.keys(Shapes)});
    this.haloRadius = this.float(PROPS.haloRadius, 12,
      {category: PROPS_CATS.APPEARANCE, editor: 'slider', min: 1, max: 32});
    this.haloWidth = this.float(PROPS.haloWidth, 4,
      {category: PROPS_CATS.APPEARANCE, editor: 'slider', min: 1, max: 32});
    this.fillColor = this.int(PROPS.fillColor, TreeDefaultPalette.Main, {category: PROPS_CATS.APPEARANCE});
    this.showShapeBorders = this.bool(PROPS.showShapeBorders, false, {category: PROPS_CATS.APPEARANCE});
    this.shapeBorderWidth = this.float(PROPS.shapeBorderWidth, 1,
      {category: PROPS_CATS.APPEARANCE, editor: 'slider', min: 0, max: 50.});
    this.shapeBorderAlpha = this.float(PROPS.shapeBorderAlpha, 0.14,
      {category: PROPS_CATS.APPEARANCE, editor: 'slider', min: 0.0, max: 1.0});
    this.highlightColor = this.int(PROPS.highlightColor, TreeDefaultPalette.Selection,
      {category: PROPS_CATS.APPEARANCE});
    this.showBranchLengths = this.bool(PROPS.showBranchLengths, false,
      {category: PROPS_CATS.APPEARANCE});
    this.showEdges = this.bool(PROPS.showEdges, true,
      {category: PROPS_CATS.APPEARANCE});
    this.showLabels = this.bool(PROPS.showLabels, false,
      {category: PROPS_CATS.APPEARANCE});
    this.showLeafLabels = this.bool(PROPS.showLeafLabels, false,
      {category: PROPS_CATS.APPEARANCE});
    this.showInternalLabels = this.bool(PROPS.showInternalLabels, false,
      {category: PROPS_CATS.APPEARANCE});
    this.fontFamily = this.string(PROPS.fontFamily, 'monospace',
      {category: PROPS_CATS.APPEARANCE});
    this.fontSize = this.float(PROPS.fontSize, 16,
      {category: PROPS_CATS.APPEARANCE, editor: 'slider', min: 2, max: 48});

    this.padding = this.float(PROPS.padding, 16,
      {category: PROPS_CATS.LAYOUT, editor: 'slider', min: 0, max: 50});
    this.treeToCanvasRatio = this.float(PROPS.treeToCanvasRatio, 0.25,
      {category: PROPS_CATS.BEHAVIOR, editor: 'slider', min: 0.01, max: 2, step: 0.01});

    this.interactive = this.bool(PROPS.interactive, true, {category: PROPS_CATS.BEHAVIOR});
    this.branchZoom = this.float(PROPS.branchZoom, 0,
      {category: PROPS_CATS.BEHAVIOR, editor: 'slider', min: -4, max: 4, step: 0.1});
    this.stepZoom = this.float(PROPS.stepZoom, 0,
      {category: PROPS_CATS.BEHAVIOR, editor: 'slider', min: -4, max: 4, step: 0.1});
    this.zoom = this.float(PROPS.zoom, 0,
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
      this.nodeCol = this.nwkDf.getCol(this.nodeColumnName);
      this.parentCol = this.nwkDf.getCol(this.parentColumnName);
    }
    // else {
    //   this.newickCol = this.dataFrame.columns.bySemType('newick');
    // }

    if (!this.viewed) {
      this.buildView();
      this.viewed = true;
    }
  }

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
      const props: { [p: string]: any } = {
        source: {type: 'biojs', data: defaultNwkRoot},
        interactive: true,
      };

      Object.assign(props, this.phylocanvasProps);
      this.viewer = new PhylocanvasGL(this.treeDiv, props);
      this.viewer.deck.setProps({
        useDevicePixels: true,
        onAfterRender: this.viewerDeckOnAfterRender.bind(this),
        onHover: this.viewerOnHover.bind(this),
      });

      this.onContextMenu.subscribe(this.onContextMenuHandler.bind(this));
      this.viewer.view.style.backgroundImage = 'none';

      // Preserve original selectNode method to call it bound from replacing method
      const viewerSelectNodeSuper = this.viewer.selectNode.bind(this.viewer);
      const viewerSelectNodeHandler = this.viewerOnSelected.bind(this);
      this.viewer.selectNode = (nodeOrId: any, append = false): void => {
        viewerSelectNodeSuper(nodeOrId, append);
        viewerSelectNodeHandler();
      };

      this.calcSize();

      // this.subs.push(rxjs.fromEvent<MouseEvent>(this.phylocanvasViewer.deck.canvas!, 'mousedown')
      //   .subscribe((args: MouseEvent) => {
      //
      //   }));
    }
    const newickRoot = this.newick ? parseNewick(this.newick) : defaultNwkRoot;
    this.viewer!.setProps({source: {type: 'biojs', data: newickRoot}});

    this.viewSubs.push(this.nwkDf.onSelectionChanged.subscribe(this.dfOnSelectionChanged.bind(this)));
    // this.phyloTreeViewer.selectNode = this.tvSelectNode.bind(this);
    //@ts-ignore
    //this.viewer.handleHover = this.viewerHandleHover.bind(this);

    // this._viewSubs.push(DG.debounce(this.dataFrame.onCurrentRowChanged, 50).subscribe((_) => this.render(false)));
  }

  private onContextMenuHandler(menu: DG.Menu) {
    menu.item('Reset view', () => {
      this.viewer?.fitInCanvas();
    });
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

    const setProps = (updater: { [p: string]: any }) => {
      Object.assign(this.phylocanvasProps, updater);
      if (this.viewer)
        this.viewer.setProps(updater);
    };

    switch (property.name) {
    case PROPS.treeType:
      const treeTypeValue = TreeTypesTypes[this.treeType];
      setProps({type: treeTypeValue});
      break;

    case PROPS.nodeShape:
      const nodeShapeValue: string = Shapes[this.nodeShape as keyof typeof Shapes];
      setProps({nodeShape: nodeShapeValue});
      break;

    case PROPS.strokeColor:
      const strokeColorValue: string = intToHtml(this.strokeColor);
      setProps({strokeColour: strokeColorValue});
      console.debug('PhyloTreeViewer: PhylocanvasGlViewer.onPropertyChanged() ' +
        `${property.name} = ${strokeColorValue} .`);
      break;

    case PROPS.fillColor:
      const fillColorValue: string = intToHtml(this.fillColor);
      setProps({fillColour: fillColorValue});
      console.debug('PhyloTreeViewer: PhylocanvasGlViewer.onPropertyChanged() ' +
        `${property.name} = ${fillColorValue} .`);
      break;

    case PROPS.highlightColor:
      const highlightColorValue: string = intToHtml(this.highlightColor);
      setProps({highlightColour: highlightColorValue});
      console.debug('PhyloTreeViewer: PhylocanvasGlViewer.onPropertyChanged() ' +
        `${property.name} = ${highlightColorValue} .`);
      break;

    case PROPS.strokeWidth:
    case PROPS.nodeSize:
    case PROPS.showShapes:
    case PROPS.haloRadius:
    case PROPS.haloWidth:
    case PROPS.branchZoom:
    case PROPS.showShapeBorders:
    case PROPS.shapeBorderWidth:
    case PROPS.shapeBorderAlpha:
    case PROPS.showBranchLengths:
    case PROPS.showEdges:
    case PROPS.showLabels:
    case PROPS.showLeafLabels:
    case PROPS.showInternalLabels:
    case PROPS.fontFamily:
    case PROPS.fontSize:
    case PROPS.padding:
    case PROPS.treeToCanvasRatio:
    case PROPS.interactive:
    case PROPS.stepZoom:
    case PROPS.zoom:
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
    this._onRendered.next();
  }

  protected viewerOnHover(info: PickingInfo, event: MjolnirPointerEvent): void {
    this._onHover.next({info, event});
  }

  /** Sync selection from nwkDf to PhylocanvasGL viewer */
  protected dfOnSelectionChanged(value: any) {
    const nodeCol: DG.Column = this.nwkDf.getCol(this.nodeColumnName);
    const selectedIds: string[] = wu.count(0).take(nodeCol.length)
      .filter((rowI) => this.nwkDf.selection.get(rowI))
      .map((rowI) => nodeCol.get(rowI))
      .toArray();
    this.viewer!.setProps({selectedIds: selectedIds});
  }

  /** Sync selection from PhylocanvasGL viewer to nwkDf */
  protected viewerOnSelected() {
    const nodeCol: DG.Column = this.nwkDf.getCol(this.nodeColumnName);
    const selectedIdSet = new Set(this.viewer!.props.selectedIds);
    this.nwkDf.selection.init((rowI) => selectedIdSet.has(nodeCol.get(rowI)));
  }

  // -- IRenderer --

  private _onRendered: Subject<void> = new Subject<void>();

  get onRendered(): Observable<void> { return this._onRendered; }

  invalidate(caller?: string): void {
    this.viewer!.render();
  }

  async awaitRendered(timeout: number | undefined = 5000): Promise<void> {
    await testEvent(this.onRendered, () => {}, () => {
      this.invalidate();
    }, timeout);
  }
}

interface TreeMap {
  index: number[];
  items: string[];
}

type TreeLeaves = { [key: string]: TreeMap };
type TreeLeavesMap = { [key: string]: number[] };
