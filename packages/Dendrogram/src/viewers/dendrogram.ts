import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';
import * as rxjs from 'rxjs';

import {Unsubscribable} from 'rxjs';
import {TAGS as treeTAGS} from '@datagrok-libraries/bio/src/trees';
import {HoverType, ITreePlacer, ITreeStyler, markupNode, MarkupNodeType, TreeStylerBase} from './tree-renderers/markup';
import {CanvasTreeRenderer} from './tree-renderers/canvas-tree-renderer';
import {TreeRendererBase} from './tree-renderers/tree-renderer-base';
import {isLeaf, ITreeHelper, parseNewick} from '@datagrok-libraries/bio';
import {RectangleTreeHoverType, RectangleTreePlacer} from './tree-renderers/rectangle-tree-placer';
import {TreeHelper} from '../utils/tree-helper';
import {toRgba, setAlpha} from '@datagrok-libraries/utils/src/color';
import {DendrogramColorCodingTreeStyler, DendrogramTreeStyler} from './tree-renderers/dendrogram-tree-styler';

export const LINE_WIDTH = 2;
export const NODE_SIZE = 4;
export const TRANS_ALPHA = 0.5;

export enum TreeColorNames {
  Main = 'Main',
  Light = 'Light',
  Current = 'Current',
  MouseOver = 'MouseOver',
  Selection = 'Selection',
}

export const TreeDefaultPalette: { [name: string]: number } = {
  [TreeColorNames.Main]: DG.Color.categoricalPalette[12],
  [TreeColorNames.Light]: DG.Color.categoricalPalette[13],
  [TreeColorNames.Current]: DG.Color.currentRow,
  [TreeColorNames.MouseOver]: DG.Color.mouseOverRows,
  [TreeColorNames.Selection]: DG.Color.selectedRows,
};

// Obtained with DG.Color.toRgb(DG.Color.categoricalPalette[12])
const categoricalPaletteList: string[] = [
  'rgb(31,119,180)',
  'rgb(255,187,120)',
  'rgb(44,160,44)',
  'rgb(152,223,138)',
  'rgb(214,39,40)',
  'rgb(255,152,150)',
  'rgb(148,103,189)',
  'rgb(197,176,213)',
  'rgb(140,86,75)',
  'rgb(196,156,148)',
  'rgb(227,119,194)',
  'rgb(247,182,210)',
  'rgb(127,127,127)',
  'rgb(199,199,199)',
  'rgb(188,189,34)',
  'rgb(219,219,141)',
  'rgb(23,190,207)',
  'rgb(158,218,229)'];

export enum PROPS_CATS {
  APPEARANCE = 'Appearance',
  BEHAVIOR = 'Behavior',
  LAYOUT = 'Layout',
  DATA = 'Data',
}

export enum PROPS {
  // -- Data --
  newick = 'newick',
  newickTag = 'newickTag',
  nodeColumnName = 'nodeColumnName',
  colorColumnName = 'colorColumnName',
  colorAggrType = 'colorAggrType',

  // -- Appearance --
  lineWidth = 'lineWidth',
  nodeSize = 'nodeSize',
  mainColor = 'mainColor',
  lightColor = 'lightColor',
  currentColor = 'currentColor',
  mouseOverColor = 'mouseOverColor',
  selectionsColor = 'selectionsColor',
  font = 'font',
  showGrid = 'showGrid',
  showLabels = 'showLabels',
  firstLeaf = 'firstLeaf',
  step = 'step',

  // -- Behavior --
  stepZoom = 'stepZoom',
}


const newickDefault: string = ';';

export interface IDendrogram {
  get placer(): RectangleTreePlacer<MarkupNodeType>;

  get renderer(): CanvasTreeRenderer<MarkupNodeType>;
}

export class Dendrogram extends DG.JsViewer implements IDendrogram {
  private viewed: boolean = false;

  // -- Data --
  [PROPS.newick]: string;
  [PROPS.newickTag]: string;
  [PROPS.nodeColumnName]: string;
  [PROPS.colorColumnName]: string;
  [PROPS.colorAggrType]: string;

  // -- Appearance --
  [PROPS.lineWidth]: number;
  [PROPS.nodeSize]: number;
  [PROPS.showGrid]: boolean;
  [PROPS.mainColor]: number;
  [PROPS.lightColor]: number;
  [PROPS.currentColor]: number;
  [PROPS.mouseOverColor]: number;
  [PROPS.selectionsColor]: number;

  [PROPS.font]: string;

  [PROPS.showLabels]: boolean;

  [PROPS.step]: number;

  // -- Behavior --
  [PROPS.stepZoom]: number;

  mainStyler: DendrogramTreeStyler;
  mainStylerOnTooltipShowSub?: rxjs.Unsubscribable;

  lightStyler: DendrogramTreeStyler;
  currentStyler: DendrogramTreeStyler;
  mouseOverStyler: DendrogramTreeStyler;
  selectionsStyler: DendrogramTreeStyler;

  constructor() {
    super();

    // -- Data --, not userEditable option is not displayed in Property panel, but can be set through setOptions()
    this.newick = this.string(PROPS.newick, newickDefault,
      {category: PROPS_CATS.DATA/*, userEditable: false*/});
    this.newickTag = this.string(PROPS.newickTag, null,
      {category: PROPS_CATS.DATA, choices: []});
    this.nodeColumnName = this.string(PROPS.nodeColumnName, null,
      {category: PROPS_CATS.DATA});
    this.colorColumnName = this.string(PROPS.colorColumnName, null,
      {category: PROPS_CATS.DATA});
    this.colorAggrType = this.string(PROPS.colorAggrType, null,
      {category: PROPS_CATS.DATA, choices: [DG.AGG.AVG, DG.AGG.MIN, DG.AGG.MAX, DG.AGG.MED, DG.AGG.TOTAL_COUNT]});

    this.lineWidth = this.float(PROPS.lineWidth, LINE_WIDTH,
      {category: PROPS_CATS.APPEARANCE, editor: 'slider', min: 0, max: 16, step: 0.1});
    this.nodeSize = this.float(PROPS.nodeSize, NODE_SIZE,
      {category: PROPS_CATS.APPEARANCE, editor: 'slider', min: 0, max: 16, step: 0.1});

    this.showGrid = this.bool(PROPS.showGrid, false, {category: PROPS_CATS.APPEARANCE});

    this.mainColor = this.int(PROPS.mainColor, TreeDefaultPalette[TreeColorNames.Main],
      {category: PROPS_CATS.APPEARANCE});

    this.lightColor = this.int(PROPS.lightColor, TreeDefaultPalette[TreeColorNames.Light],
      {category: PROPS_CATS.APPEARANCE});

    this.currentColor = this.int(PROPS.currentColor, TreeDefaultPalette[TreeColorNames.Current],
      {category: PROPS_CATS.APPEARANCE});

    this.mouseOverColor = this.int(PROPS.mouseOverColor, TreeDefaultPalette[TreeColorNames.MouseOver],
      {category: PROPS_CATS.APPEARANCE});

    this.selectionsColor = this.int(PROPS.selectionsColor, TreeDefaultPalette[TreeColorNames.Selection],
      {category: PROPS_CATS.APPEARANCE});

    this.showLabels = this.bool(PROPS.showLabels, false, {category: PROPS_CATS.APPEARANCE});

    this.font = this.string(PROPS.font, 'monospace 10pt', {category: PROPS_CATS.APPEARANCE});

    this.stepZoom = this.float(PROPS.stepZoom, 0,
      {category: PROPS_CATS.BEHAVIOR, editor: 'slider', min: -4, max: 4, step: 0.1});

    this.step = this.float(PROPS.step, 28,
      {category: PROPS_CATS.APPEARANCE, editor: 'slider', min: 0, max: 64, step: 0.1});

    this.mainStyler = this.getMainStyler();
    this.lightStyler = new DendrogramTreeStyler('light',
      this.lineWidth, this.nodeSize, false,
      toRgba(setAlpha(this.lightColor, TRANS_ALPHA)),
      toRgba(setAlpha(this.lightColor, TRANS_ALPHA)));
    this.currentStyler = new DendrogramTreeStyler('current',
      this.lineWidth, this.nodeSize, false,
      toRgba(setAlpha(this.currentColor, TRANS_ALPHA)),
      toRgba(setAlpha(this.currentColor, TRANS_ALPHA)));
    this.mouseOverStyler = new DendrogramTreeStyler('mouseOver',
      this.lineWidth, this.nodeSize, false,
      toRgba(setAlpha(this.mouseOverColor, TRANS_ALPHA)),
      toRgba(setAlpha(this.mouseOverColor, TRANS_ALPHA)));
    this.selectionsStyler = new DendrogramTreeStyler('selections',
      this.lineWidth, this.nodeSize, false,
      toRgba(setAlpha(this.selectionsColor, TRANS_ALPHA)),
      toRgba(setAlpha(this.selectionsColor, TRANS_ALPHA)));
  }

  private _newick: string = '';

  // effective tree value (to plot)
  private treeNewick: string | null = null;

  override onTableAttached() {
    super.onTableAttached();
    console.debug('Dendrogram: PhylocanvasGlViewer.onTableAttached() ' +
      `this.dataFrame = ${!this.dataFrame ? 'null' : 'value'} )`);

    // -- Editors --
    // update editors for properties dependent from viewer's dataFrame
    const dfTagNameList = wu<string>(this.dataFrame.tags.keys())
      .filter((tagName: string) => tagName.startsWith('.')).toArray();
    this.props.getProperty(PROPS.newickTag).choices = ['', ...dfTagNameList];

    this.setData();
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
      console.warn('Dendrogram: PhylocanvasGlViewer.onPropertyChanged() No property value');
      return;
    }

    switch (property.name) {
    case PROPS.lineWidth:
      this.mainStyler.lineWidth = this.lineWidth;
      this.lightStyler.lineWidth = this.lineWidth;
      this.currentStyler.lineWidth = this.lineWidth;
      this.mouseOverStyler.lineWidth = this.lineWidth;
      this.selectionsStyler.lineWidth = this.lineWidth;
      break;

    case PROPS.nodeSize:
      this.mainStyler.nodeSize = this.nodeSize;
      this.lightStyler.nodeSize = this.nodeSize;
      this.currentStyler.nodeSize = this.nodeSize;
      this.mouseOverStyler.nodeSize = this.nodeSize;
      this.selectionsStyler.nodeSize = this.nodeSize;
      break;

    case PROPS.showGrid:
      this.mainStyler.showGrid = this.showGrid;
      break;

    case PROPS.mainColor:
      this.mainStyler.strokeColor = toRgba(setAlpha(this.mainColor, TRANS_ALPHA));
      this.mainStyler.fillColor = toRgba(setAlpha(this.mainColor, TRANS_ALPHA));
      break;

    case PROPS.lightColor:
      this.lightStyler.strokeColor = toRgba(setAlpha(this.lightColor, TRANS_ALPHA));
      this.lightStyler.fillColor = toRgba(setAlpha(this.lightColor, TRANS_ALPHA));
      break;

    case PROPS.currentColor:
      this.currentStyler.strokeColor = toRgba(setAlpha(this.currentColor, TRANS_ALPHA));
      this.currentStyler.fillColor = toRgba(setAlpha(this.currentColor, TRANS_ALPHA));
      break;

    case PROPS.selectionsColor:
      this.selectionsStyler.strokeColor = toRgba(setAlpha(this.selectionsColor, TRANS_ALPHA));
      this.selectionsStyler.fillColor = toRgba(setAlpha(this.selectionsColor, TRANS_ALPHA));
      break;

    case PROPS.font:
      break;

    case PROPS.colorColumnName:
    case PROPS.colorAggrType:
      if (this.viewed)
        this.mainStylerOnTooltipShowSub!.unsubscribe();

      this._renderer!.mainStyler = this.mainStyler = this.getMainStyler();

      if (this.viewed)
        this.mainStylerOnTooltipShowSub = this.mainStyler.onTooltipShow.subscribe(this.stylerOnTooltipShow.bind(this));

      break;
    }

    // Rebuild view
    switch (property.name) {
    case PROPS.newick:
    case PROPS.newickTag:
      this.setData();
      break;
    }
  }

  // -- Data --

  private setData(): void {
    if (this.viewed) {
      this.destroyView();
      this.viewed = false;
    }

    // -- Tree data --
    // Tree newick data source priorities
    // this.newick (property)                  // the highest priority
    // this.dataFrame.getTag(this.newickTag)
    // this.dataFrame.getTag(TREE_TAGS.NEWICK) // the lowest priority
    let newickTag: string = treeTAGS.NEWICK;
    if (this.newickTag) newickTag = this.newickTag;
    this.treeNewick = this.dataFrame.getTag(newickTag);
    if (this.newick && this.newick != newickDefault) this.treeNewick = this.newick;

    if (!this.viewed) {
      this.buildView();
      this.viewed = true;
    }
  }

  // -- View --

  private viewSubs: Unsubscribable[] = [];

  private treeDiv?: HTMLDivElement;

  private _placer?: RectangleTreePlacer<MarkupNodeType>;

  get placer(): RectangleTreePlacer<MarkupNodeType> { return this._placer!; }

  private _renderer?: CanvasTreeRenderer<MarkupNodeType>;

  get renderer(): CanvasTreeRenderer<MarkupNodeType> { return this._renderer!; }

  private destroyView(): void {
    console.debug('Dendrogram: Dendrogram.destroyView()');

    this.mainStylerOnTooltipShowSub!.unsubscribe();
    for (const sub of this.viewSubs) sub.unsubscribe();
    this.viewSubs = [];

    this._renderer!.detach();
    delete this._renderer;

    delete this._placer;

    this.treeDiv!.remove();
    delete this.treeDiv;
  }

  private buildView(): void {
    console.debug('Dendrogram: Dendrogram.buildView()');

    this.treeDiv = ui.div([], {
      style: {
        width: `${this.root.clientWidth}px`,
        height: `${this.root.clientHeight}px`,
        backgroundColor: '#A0A0FF',
      }
    });
    this.treeDiv.style.setProperty('overflow', 'hidden', 'important');
    this.root.appendChild(this.treeDiv);

    const treeRoot: MarkupNodeType = parseNewick(this.treeNewick!) as MarkupNodeType;
    markupNode(treeRoot);
    const totalLength: number = treeRoot.subtreeLength!;
    this.mainStylerOnTooltipShowSub = this.mainStyler.onTooltipShow.subscribe(this.stylerOnTooltipShow.bind(this));
    this._placer = new RectangleTreePlacer<MarkupNodeType>(
      treeRoot.minIndex - 0.5, treeRoot.maxIndex + 0.5, totalLength);
    this._renderer = new CanvasTreeRenderer(
      treeRoot, this._placer,
      this.mainStyler, this.lightStyler,
      this.currentStyler, this.mouseOverStyler, this.selectionsStyler);
    this.viewSubs.push(this._renderer.onCurrentChanged.subscribe(this.rendererOnCurrentChanged.bind(this)));
    this.viewSubs.push(this._renderer.onMouseOverChanged.subscribe(this.rendererOnMouseOverChanged.bind(this)));
    this.viewSubs.push(this._renderer.onSelectionChanged.subscribe(this.rendererOnSelectionChanged.bind(this)));
    this._renderer.attach(this.treeDiv);

    this.viewSubs.push(ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));
    // this.viewSubs.push(this.mainStyler.onStylingChanged.subscribe(this.stylerOnStylingChanged.bind(this)));
    this.viewSubs.push(this.mainStyler.onTooltipShow.subscribe(this.stylerOnTooltipShow.bind(this)));

    this.viewSubs.push(this.dataFrame.onSelectionChanged.subscribe(this.dataFrameOnSelectionChanged.bind(this)));
    this.viewSubs.push(this.dataFrame.onCurrentRowChanged.subscribe(this.dataFrameOnCurrentRowChanged.bind(this)));
    this.viewSubs.push(this.dataFrame.onMouseOverRowChanged.subscribe(this.dataFrameOnMouseOverRowChanged.bind(this)));
  }

  // -- Handle controls events --

  private rootOnSizeChanged(): void {
    console.debug('Dendrogram: Dendrogram.rootOnSizeChanged()');

    this.treeDiv!.style.width = `${this.root.clientWidth}px`;
    this.treeDiv!.style.height = `${this.root.clientHeight}px`;
  }

  private rendererOnCurrentChanged(): void {
    window.setTimeout(() => {
      if (!this._renderer) return;

      const oldCurrentRowIdx: number = this.dataFrame.currentRowIdx;
      if (!this._renderer.currentNode) {
        this.dataFrame.currentRowIdx = -1;
      } else {
        if (this.nodeColumnName) {
          const rowCount = this.dataFrame.rowCount;
          const nodeCol: DG.Column = this.dataFrame.getCol(this.nodeColumnName);
          const nodeName = this._renderer.currentNode.name;

          const newCurrentRowIdx = wu.count(0).take(rowCount)
            .find((rowI: number) => nodeCol.get(rowI) == nodeName) ?? -1;

          if (newCurrentRowIdx != oldCurrentRowIdx)
            this.dataFrame.currentRowIdx = newCurrentRowIdx;
        }
      }
    });
  }

  private rendererOnMouseOverChanged(): void {
    window.setTimeout(() => {
      if (!this._renderer) return;

      const oldMouseOverRowIdx: number = this.dataFrame.mouseOverRowIdx;
      if (!this._renderer.mouseOverNode) {
        this.dataFrame.mouseOverRowIdx = -1;
      } else {
        if (this.nodeColumnName) {
          const rowCount = this.dataFrame.rowCount;
          const nodeCol: DG.Column = this.dataFrame.getCol(this.nodeColumnName);
          const nodeName = this._renderer.mouseOverNode.name;

          const newMouseOverRowIdx = wu.count(0).take(rowCount)
            .find((rowI: number) => nodeCol.get(rowI) == nodeName) ?? -1;

          if (newMouseOverRowIdx != oldMouseOverRowIdx)
            this.dataFrame.mouseOverRowIdx = newMouseOverRowIdx;
        }
      }
    });
  }

  private rendererOnSelectionChanged(): void {
    window.setTimeout(() => {
      if (!this._renderer) return;

      const oldSelection: DG.BitSet = this.dataFrame.selection.clone();
      if (this._renderer.selectedNodes.length == 0) {
        this.dataFrame.selection.init((rowI) => { return false; }, false);
      } else {
        if (this.nodeColumnName) {
          const nodeCol: DG.Column = this.dataFrame.getCol(this.nodeColumnName);
          const th = new TreeHelper();
          const nodeNameSet = new Set(
            this._renderer.selectedNodes
              .map((sn) => th.getNodeList(sn).map((n) => n.name))
              .flat());
          console.debug('Dendrogram: Dendrogram.rendererOnSelectionChanged(), ' +
            `nodeNameSet = ${JSON.stringify([...nodeNameSet])}`);

          this.dataFrame.selection.init(
            (rowI) => {
              const nodeName = nodeCol.get(rowI);
              return nodeNameSet.has(nodeName);
            },
            false);
        }
      }

      const newSelection: DG.BitSet = this.dataFrame.selection;
      let selectionChanged: boolean = oldSelection.length !== newSelection.length || // != -> true (changed)
        oldSelection.trueCount !== newSelection.trueCount;
      if (!selectionChanged) {
        for (let rowI: number = 0; rowI < oldSelection.length; rowI++) {
          if (oldSelection.get(rowI) !== newSelection.get(rowI)) {
            selectionChanged = true;
            break;
          }
        }
      }
      if (selectionChanged)
        this.dataFrame.selection.fireChanged();
    }, 0 /* next event cycle*/);
  }

  private dataFrameOnSelectionChanged(value: any): void {
    if (!this._renderer || !this._placer) return;

    /* Here we get selected rows from dataFrame.
     * Some of selected nodes can be in subtree of others.
     * We need to merge nodes to form selections object.
     * Node rows can be selected or deselected.
     */
    const th: TreeHelper = new TreeHelper();
    if (this.nodeColumnName) {
      const nodeList: MarkupNodeType[] = th.getNodeList(this._renderer.treeRoot);
      const nodeDict: { [name: string]: MarkupNodeType } = {};
      for (const node of nodeList) {
        if (node.name in nodeDict)
          throw new Error('Non unique key tree node name');
        nodeDict[node.name] = node;
      }

      const rowDict: { [name: string]: number } = {};
      const nodeCol: DG.Column = this.dataFrame.getCol(this.nodeColumnName);
      for (let rowI = 0; rowI < nodeCol.length; rowI++) {
        const nodeName: string = nodeCol.get(rowI);
        rowDict[nodeName] = rowI;
      }

      const selections: RectangleTreeHoverType<MarkupNodeType>[] = [];
      const selectedIndexes = this.dataFrame.selection.getSelectedIndexes();
      for (const selRowI of selectedIndexes) {
        const nodeName: string = nodeCol.get(selRowI);
        const node: MarkupNodeType = nodeDict[nodeName];

        // add new node to selections structure
        // if a new node is subnode of one of the nodes already in selection so skip add
        // if a there are nodes in selections which are subnodes of adding so we should remove existing and add new

        const toDeleteList: number[] = [];
        let addSkip: boolean = false;
        for (let selI = 0; selI < selections.length; selI++) {
          const selNode = selections[selI].node;

          // if (th.includes(selNode, node))
          if (
            ((selNode.minIndex ?? selNode.index) <= (node.minIndex ?? node.index)) &&
            ((node.maxIndex ?? node.index) <= (selNode.maxIndex ?? selNode.index))) {
            // The adding node is in the subtree of already selected node
            addSkip = true;
          } else if ( /* if (th.includes(node, selNode)) */
            ((node.minIndex ?? node.index) <= (selNode.minIndex ?? selNode.index)) &&
            ((selNode.maxIndex ?? selNode.index) <= (node.maxIndex ?? node.index))) {
            // The already selected node is in the subtree of the adding node
            toDeleteList.push(selI);
          }
        }
        for (let toDeleteI = toDeleteList.length - 1; toDeleteI >= 0; toDeleteI--)
          selections.splice(toDeleteList[toDeleteI]);

        if (!addSkip)
          selections.push({node: node, nodeHeight: this._placer.getNodeHeight(this._renderer.treeRoot, node)!});
      }

      function selectSubs(df: DG.DataFrame, rowDict: { [name: string]: number }, node: MarkupNodeType) {
        const rowI = rowDict[node.name];
        df.selection.set(rowI, true, false);
        for (const childNode of (node.children ?? []))
          selectSubs(df, rowDict, childNode);
      }

      for (const selection of selections)
        selectSubs(this.dataFrame, rowDict, selection.node);

      this._renderer.selections = selections;
    }
  }

  private dataFrameOnCurrentRowChanged(value: any): void {
    if (!this._renderer) return;

    if (this.nodeColumnName) {
      const nodeCol: DG.Column = this.dataFrame.getCol(this.nodeColumnName);
      const currentNodeName = nodeCol.get(this.dataFrame.currentRowIdx);

      const th: ITreeHelper = new TreeHelper();
      const currentNode: MarkupNodeType | null = th.getNodeList(this._renderer.treeRoot)
        .find((node) => currentNodeName == node.name) ?? null;
      const current: RectangleTreeHoverType<MarkupNodeType> | null = currentNode ? {
        node: currentNode,
        nodeHeight: this._placer!.getNodeHeight(this._renderer.treeRoot, currentNode)!
      } : null;

      this._renderer.current = current;
    }
  }

  private dataFrameOnMouseOverRowChanged(value: any): void {
    if (!this._renderer) return;

    if (this.nodeColumnName) {
      const nodeCol: DG.Column = this.dataFrame.getCol(this.nodeColumnName);
      const mouseOverNodeName: string = nodeCol.get(this.dataFrame.mouseOverRowIdx);

      const th: ITreeHelper = new TreeHelper();
      const mouseOverNode: MarkupNodeType | null = th.getNodeList(this._renderer.treeRoot)
        .find((node) => mouseOverNodeName == node.name) ?? null;
      const mouseOver: RectangleTreeHoverType<MarkupNodeType> | null = mouseOverNode ? {
        node: mouseOverNode,
        nodeHeight: this._placer!.getNodeHeight(this._renderer.treeRoot, mouseOverNode)!
      } : null;

      this._renderer.mouseOver = mouseOver;
    }
  }

  private stylerOnTooltipShow({node, e}: { node: MarkupNodeType, e: MouseEvent }): void {
    if (node) {
      const minMaxIndexStr: string = !isLeaf(node) ? ` (min: ${node.minIndex}, max: ${node.maxIndex})` : '';
      const tooltip = ui.divV([
        ui.div(`${node.name}`),
        ui.div(`index: ${node.index}${minMaxIndexStr}`),
        ui.div(`desc: ${node.desc}`)]);
      ui.tooltip.show(tooltip, e.clientX + 16, e.clientY + 16);
    } else {
      ui.tooltip.hide();
    }
  }

  // private stylerOnStylingChanged() {
  //   this.lightStyler.lineWidth = this.getHighlightStylerLineWidth();
  //   this.lightStyler.nodeSize = this.getHighlightStylerNodeSize();
  //
  //   this.selectedStyler.lineWidth = this.getSelectionStylerLineWidth();
  //   this.selectedStyler.nodeSize = this.getSelectionStylerNodeSize();
  // }

  // getHighlightStylerLineWidth(): number {
  //   return Math.max(this.mainStyler.lineWidth + 4, this.mainStyler.lineWidth * 1.8);
  // }
  //
  // getHighlightStylerNodeSize(): number {
  //   return Math.max(this.mainStyler.nodeSize + 4, this.mainStyler.nodeSize * 1.8);
  // }
  //
  // getSelectionStylerLineWidth(): number {
  //   return Math.max(this.mainStyler.lineWidth + 2, this.mainStyler.lineWidth * 1.4);
  // }
  //
  // getSelectionStylerNodeSize(): number {
  //   return Math.max(this.mainStyler.nodeSize + 2, this.mainStyler.nodeSize * 1.4);
  // }

  getMainStyler(): DendrogramTreeStyler {
    let res: DendrogramTreeStyler;
    if (this.colorColumnName) {
      const colorCol: DG.Column = this.dataFrame.getCol(this.colorColumnName);
      const nodeCol: DG.Column = this.dataFrame.getCol(this.nodeColumnName);
      res = new DendrogramColorCodingTreeStyler('main-color-coding',
        this.lineWidth, this.nodeSize, this.showGrid,
        nodeCol, colorCol, this.colorAggrType,
        toRgba(setAlpha(this.mainColor, TRANS_ALPHA)),
        toRgba(setAlpha(this.mainColor, TRANS_ALPHA)));
    } else {
      res = new DendrogramTreeStyler('main',
        this.lineWidth, this.nodeSize, this.showGrid,
        toRgba(setAlpha(this.mainColor, TRANS_ALPHA)),
        toRgba(setAlpha(this.mainColor, TRANS_ALPHA)));
    }
    return res;
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
