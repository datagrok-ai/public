import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import * as bio from '@datagrok-libraries/bio';
import wu from 'wu';

import * as rxjs from 'rxjs';
import {JsViewer} from 'datagrok-api/dg';
import {Unsubscribable} from 'rxjs';
import {TREE_TAGS} from '../consts';
import {markupNode, MarkupNodeType} from './tree-renderers/markup';
import {LeafRangeGridTreeRenderer} from './tree-renderers/grid-tree-renderer';
import {CanvasTreeRenderer} from './tree-renderers/canvas-tree-renderer';
import {TreeRendererBase} from './tree-renderers/tree-renderer-base';

enum PROPS_CATS {
  APPEARANCE = 'Appearance',
  BEHAVIOR = 'Behavior',
  LAYOUT = 'Layout',
  DATA = 'Data',
}

enum Props {
  lineWidth = 'lineWidth',
  nodeSize = 'nodeSize',
  strokeColor = 'strokeColor',
  fillColor = 'fillColor',
  font = 'font',

  showGrid = 'showGrid',
  showLabels = 'showLabels',

  firstLeaf = 'firstLeaf',
  step = 'step',
  stepZoom = 'stepZoom',

  newick = 'newick',
}

export class Dendrogram extends DG.JsViewer {
  private viewed: boolean = false;

  [Props.lineWidth]: number;
  [Props.nodeSize]: number;
  [Props.strokeColor]: number;
  [Props.fillColor]: number;
  [Props.font]: string;

  [Props.showGrid]: boolean;
  [Props.showLabels]: boolean;

  [Props.firstLeaf]: number;
  [Props.step]: number;
  [Props.stepZoom]: number;

  constructor() {
    super();

    this.lineWidth = this.float(Props.lineWidth, 1, {category: PROPS_CATS.APPEARANCE,});
    this.nodeSize = this.float(Props.nodeSize, 14,
      {category: PROPS_CATS.APPEARANCE, editor: 'slider', min: 1, max: 32,});
    this.strokeColor = this.int(Props.strokeColor, 0x222222, {category: PROPS_CATS.APPEARANCE});
    this.fillColor = this.int(Props.fillColor, 0x333333, {category: PROPS_CATS.APPEARANCE,});

    this.showGrid = this.bool(Props.showGrid, false, {category: PROPS_CATS.APPEARANCE,});

    this.showLabels = this.bool(Props.showLabels, false, {category: PROPS_CATS.APPEARANCE});

    this.font = this.string(Props.font, 'monospace 10pt', {category: PROPS_CATS.APPEARANCE});

    this.stepZoom = this.float(Props.stepZoom, 0,
      {category: PROPS_CATS.BEHAVIOR, editor: 'slider', min: -4, max: 4, step: 0.1});

    this.step = this.float('step', 28,
      {category: PROPS_CATS.APPEARANCE, editor: 'slider', min: 0, max: 64, step: 0.1});

    // data
    this.newick = this.string('newick', ';', {category: PROPS_CATS.DATA});
  }

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
    // TODO: Logic to get newick
    this.newick = this.nwkDf.getTag(TREE_TAGS.NEWICK) ?? '';

    if (!this.viewed) {
      this.buildView();
      this.viewed = true;
    }
  }

  private _newick: string;

  public get newick(): string { return this._newick; }

  public set newick(value: string) {
    this._newick = value;
    if (this.renderer) {
      const treeRoot = bio.Newick.parse_newick(this._newick);
      markupNode(treeRoot);
      this.renderer!.treeRoot = treeRoot;
    }
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

    super.detach();
  }

  override onPropertyChanged(property: DG.Property | null) {
    super.onPropertyChanged(property);

    if (!property) {
      console.warn('PhyloTreeViewer: PhylocanvasGlViewer.onPropertyChanged() No property value');
      return;
    }

    switch (property.name) {
    case 'newick':
      if (this.viewed) {
        this.destroyView();
        this.viewed = true;
      }

      if (this.viewed) {
        this.viewed = false;
        this.buildView();
      }
      break;

    case Props.lineWidth:
    case Props.nodeSize:
    case Props.strokeColor:
    case Props.fillColor:
    case Props.font:
      this.render();
      break;
    }
  }

  // -- View --

  private viewSubs: Unsubscribable[];

  private treeDiv?: HTMLDivElement;
  private renderer?: TreeRendererBase<MarkupNodeType>;

  private destroyView(): void {
    console.debug('PhyloTreeViewer: Dendrogram.destroyView()');

    delete this.renderer;
    delete this.treeDiv;
    for (const sub of this.viewSubs) sub.unsubscribe();
  }

  private buildView(): void {
    console.debug('PhyloTreeViewer: Dendrogram.buildView()');

    this.treeDiv = ui.div([], {
      style: {
        width: `${this.root.clientWidth}px`,
        height: `${this.root.clientHeight}px`,
        backgroundColor: '#A0A0FF',
      }
    });
    this.root.appendChild(this.treeDiv);

    const treeRoot: MarkupNodeType = bio.Newick.parse_newick(this.newick);
    markupNode(treeRoot);
    const totalLength: number = treeRoot.subtreeLength!;
    this.renderer = new CanvasTreeRenderer(treeRoot, totalLength, this.treeDiv);

    this.viewSubs.push(
      ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));
  }

  // -- render --

  private render(): void {
    if (!this.renderer) return;

    this.renderer.render();
  }

  // -- Handle controls events --

  private rootOnSizeChanged() {
    console.debug('PhyloTreeViewer: Dendrogram.rootOnSizeChanged()');

    this.treeDiv!.style.width = `${this.root.clientWidth}px`;
    this.treeDiv!.style.height = `${this.root.clientHeight}px`;
  }
}

export class MyViewer extends DG.JsViewer {
  prop1: number;
  propMax: number;

  constructor() {
    super();

    this.prop1 = this.float('prop1', 1,
      {category: 'cat', editor: 'slider', min: 0, max: 64, step: 0.1});
    this.propMax = this.float('propMax', 16);
  }

  onPropertyChanged(property: DG.Property | null) {
    super.onPropertyChanged(property);

    if (property) {
      switch (property.name) {
      case 'propMax':
        //this.props.get('prop1');
        this.props.getProperty('prop1').options.max = this.propMax;
        break;
      default:
        throw new Error('Unhandled property changed \'${}\'.');
      }
    }
  }

  /*
  //name: Template
  //description: Hello world script
  //language: javascript

  const tv = grok.shell.tableView('Table');
  const v1 = await tv.dataFrame.plot.fromType('MyViewer', {});

  alert(v1);

  tv.dockManager.dock(v1, 'right');
  */
}
