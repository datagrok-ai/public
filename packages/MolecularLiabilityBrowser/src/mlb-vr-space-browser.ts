import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {toRgba} from '@datagrok-libraries/utils/src/color';
import {reduceDimensinalityWithNormalization} from '@datagrok-libraries/ml/src/sequence-space';
import {StringMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {Unsubscribable} from 'rxjs';
import {ScatterPlotViewer} from 'datagrok-api/dg';
import {getVId} from './utils/tree-stats';
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';
import {cleanMlbNewick} from './mlb-tree';
import {MlbDataFrame, TreeDataFrame} from './types/dataframe';
import {isLeaf, NodeType} from '@datagrok-libraries/bio/src/trees';
import {getTreeHelper, ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {VdRegion} from '@datagrok-libraries/bio/src/vd-regions';
import {parseNewick} from '@datagrok-libraries/bio/src/trees/phylocanvas';

enum EMBED_COL_NAMES {
  X = 'Embed_X',
  Y = 'Embed_Y'
}

/**  */
type LineType = { a: DG.Point, b: DG.Point };

type PointNodeType = NodeType & { vId: string, point?: DG.Point, screenPoint?: DG.Point };

export enum TreeColorNames {
  Main = 'Main',
}

export const TreeDefaultPalette: { [name: string]: number } = {
  [TreeColorNames.Main]: DG.Color.categoricalPalette[12],
};

export function getZoomCoordinates(W0: number, H0: number, x1: number, y1: number, x2: number, y2: number) {
  const W1 = Math.abs(x1 - x2);
  const H1 = Math.abs(y1 - y2);
  const scaleW = W0 / W1;
  const scaleH = H0 / H1;
  const scale = Math.min(scaleW, scaleH);
  const W2 = (W0 / scale) * 5;
  const H2 = (H0 / scale) * 5;
  const left = x1 < x2 ? x1 : x2;
  const top = y1 > y2 ? y1 : y2;
  const zoomLeft = (left + W1 / 2) - W2 / 2;
  const zoomRight = zoomLeft + W2;
  const zoomTop = (top - H1 / 2) + H2 / 2;
  const zoomBottom = zoomTop - H2;
  return {zoomLeft: zoomLeft, zoomRight: zoomRight, zoomTop: zoomTop, zoomBottom: zoomBottom};
}

// export function checkCursorOnLine(event: any, canvas: any, lines: ILine[]): ILine | null {
//   const rect = canvas.getBoundingClientRect();
//   const x = event.clientX - rect.left;
//   const y = event.clientY - rect.top;
//   let closestLine = null;
//   let minDist = 0;
//   for (const line of lines) {
//     const dist =
//       Math.abs(Math.hypot(line.a[0] - x, line.a[1] - y) +
//         Math.hypot(line.b[0] - x, line.b[1] - y) - Math.hypot(line.a[0] - line.b[0], line.a[1] - line.b[1]));
//     if ((!minDist && dist < 2) || dist < minDist) {
//       minDist = dist;
//       closestLine = line;
//     }
//   }
//   return closestLine;
// }

export enum VrSpaceMethodName {
  UMAP = 'UMAP',
  tSNE = 't-SNE'
}

/**  */
export class MlbVrSpaceBrowser {
  private viewed: boolean = false;

  private th: ITreeHelper;

  private _mlbDf: MlbDataFrame;
  get mlbDf(): MlbDataFrame { return this._mlbDf; }

  private _treeDf: TreeDataFrame;
  get treeDf(): TreeDataFrame { return this._treeDf; }

  private _regions: VdRegion[] = [];
  get regions(): VdRegion[] { return this._regions; }

  private _method: VrSpaceMethodName = Object.values(VrSpaceMethodName)[0];
  get method(): VrSpaceMethodName { return this._method; }

  private scatterPlot: DG.ScatterPlotViewer;
  private context: CanvasRenderingContext2D;

  private vrRows: { [vId: string]: number };
  private embedCols: { [embed: string]: DG.Column };

  /** lines - with canvas coordinates */
  private treeList: PointNodeType[] = [];
  private currentTreeIdx: number = -1;
  private mouseOverTreeIdx: number = -1;

  get root(): HTMLElement { return this.scatterPlot!.root; }

  constructor() {}

  async init(): Promise<void> {
    this.th = await getTreeHelper();

    this._mlbDf = MlbVrSpaceBrowser.emptyMlbDf;
    this.scatterPlot = await this._mlbDf.plot.fromType(DG.VIEWER.SCATTER_PLOT, {
      'xColumnName': EMBED_COL_NAMES.X,
      'yColumnName': EMBED_COL_NAMES.Y,
      'lassoTool': true,
    }) as DG.ScatterPlotViewer;
    const canvas: HTMLCanvasElement = this.scatterPlot.getInfo()['canvas'];
    this.context = canvas.getContext('2d')!;

    this.viewSubs.push(this.scatterPlot.onBeforeDrawScene.subscribe(this.scatterPlotOnBeforeDrawScene.bind(this)));
    this.viewSubs.push(this.scatterPlot.onAfterDrawScene.subscribe(this.scatterPlotOnAfterDrawScene.bind(this)));
    const k = 11;
  }

  /**
   * @param mlbDf null for default empty data frame with VRs/MLB
   * @param treeDf null for default empty data frame with Trees
   * @param regions
   * @param method
   */
  async setData(
    mlbDf: MlbDataFrame | null, treeDf: TreeDataFrame | null, regions: VdRegion[], method: VrSpaceMethodName
  ): Promise<void> {
    if (this.viewed) {
      await this.destroyView();
      this.viewed = false;
    }

    this._mlbDf = mlbDf != null ? mlbDf : MlbVrSpaceBrowser.emptyMlbDf;
    this._treeDf = treeDf != null ? treeDf : MlbVrSpaceBrowser.emptyTreeDf;
    this._regions = regions;
    this._method = method;

    this.viewSubs = [];

    if (!this.viewed) {
      await this.buildView();
      this.viewed = true;
    }
  }

  // -- View --

  private viewSubs: Unsubscribable[] = [];

  async destroyView(): Promise<void> {
    for (const sub of this.viewSubs) sub.unsubscribe();

    if (this.scatterPlot) {
      checkEmbedCols(Object.values(EMBED_COL_NAMES), MlbVrSpaceBrowser.emptyMlbDf);
      this.scatterPlot.dataFrame = MlbVrSpaceBrowser.emptyMlbDf;
    }

    this.vrRows = {};
    this.embedCols = {};
  }

  async buildView(): Promise<void> {
    try {
      if (this.mlbDf.rowCount == 0) {
        const k = 11;
        // No VRs
      } else if (this.mlbDf.rowCount == 1) {
        const k = 11;
        // Single VR to plot
      } else {
        const seqCols: { [chain: string]: DG.Column } = Object.assign({}, ...['Heavy', 'Light']
          .map((chain) => ({[chain]: this.mlbDf.getCol(`${chain} chain sequence`)})));
        const seqList = MlbVrSpaceBrowser.getSeqList(this.mlbDf, seqCols, this.regions);

        const redDimRes = await reduceDimensinalityWithNormalization(
          seqList, this.method, StringMetricsNames.Manhattan, {});
        const embedColNameList = Object.values(EMBED_COL_NAMES);
        for (let embedI: number = 0; embedI < embedColNameList.length; embedI++) {
          const embedColData: Float32Array = redDimRes.embedding[embedI];
          const embedColName: string = embedColNameList[embedI];
          const embedCol: DG.Column | null = this.mlbDf.col(embedColName);
          if (embedCol) {
            // TODO: User DG.Column.setRawData()
            // embedCol.setRawData(embedColData);
            embedCol.init((rowI) => { return embedColData[rowI]; });
          } else {
            // Notification is required to reflect added data frame Embed_<X> columns to grid columns
            // MolecularLiabilityBrowser.setView() corrects grid columns' names with .replace('_', ' ');
            const notify: boolean = embedI == embedColNameList.length - 1; // notify on adding last Embed_<X> column
            this.mlbDf.columns.add(DG.Column.fromFloat32Array(embedColName, embedColData), notify);
          }
        }

        this.vrRows = {};
        const vrRowCount: number = this.mlbDf.rowCount;
        const vIdCol: DG.Column = this.mlbDf.getCol('v id');
        for (let vrRowI = 0; vrRowI < vrRowCount; vrRowI++) {
          const vId = vIdCol.get(vrRowI);
          this.vrRows[vId] = vrRowI;
        }

        this.embedCols = {};
        for (const embedColName of Object.values(EMBED_COL_NAMES)) {
          const embedCol: DG.Column = this.mlbDf.getCol(embedColName);
          this.embedCols[embedColName] = embedCol;
        }

        const treeCol: DG.Column = this.treeDf.getCol('TREE');
        this.treeList = MlbVrSpaceBrowser.buildTreeList(this.th, treeCol, this.vrRows, this.embedCols);
      }
      checkEmbedCols(Object.values(EMBED_COL_NAMES), this.mlbDf);
      this.scatterPlot.dataFrame = this.mlbDf;
      this.scatterPlot.props.title = `VR space (method: ${this.method})`;
      this.viewSubs.push(this.treeDf.onMouseOverRowChanged.subscribe(this.treeDfOnMouseOverRowChanged.bind(this)));
      this.viewSubs.push(this.treeDf.onCurrentRowChanged.subscribe(this.treeDfOnCurrentRowChanged.bind(this)));
    } catch (err: any) {
      const errMsg: string = errorToConsole(err);
      console.debug('MLB: MlbVrSpaceBrowser.buildView() error:\n' + errMsg);
      throw err;
    }
  }

  // -- Handle events --

  private scatterPlotOnBeforeDrawScene(): void {
    console.debug('MLB: MlbVrSpaceBrowser.scatterPlotOnBeforeDrawScene() ');
    this.context.save();
    try {
      const mainStyler = new EdgeStyler(1, '#888888');
      for (const treeRoot of this.treeList)
        renderAndUpdateNode(this.context, this.scatterPlot, treeRoot, mainStyler);
    } catch (err: any) {
      const errMsg = errorToConsole(err);
      console.error('MLB: MlbVrSpaceBrowser.scatterPlotOnBeforeDrawScene() error:\n' + errMsg);
      //throw err; // Do not throw to prevent disabling event handler
    } finally {
      this.context.restore();
    }
  }

  private scatterPlotOnAfterDrawScene(): void {
    //console.debug('MLB: MlbVrSpaceBrowser.scatterPlotOnAfterDrawScene() ');
    this.context.save();
    try {
      if (this.mouseOverTreeIdx != -1) {
        const mouseOverStyler = new EdgeStyler(1.5, 'black');
        const mouseOverTreeRoot = this.treeList[this.mouseOverTreeIdx];
        renderAndUpdateNode(this.context, this.scatterPlot, mouseOverTreeRoot, mouseOverStyler);
      }
      if (this.currentTreeIdx != -1) {
        const currentStyler = new EdgeStyler(1.2, 'darkgreen');
        const currentTreeRoot = this.treeList[this.currentTreeIdx];
        renderAndUpdateNode(this.context, this.scatterPlot, currentTreeRoot, currentStyler);
      }
    } catch (err: any) {
      const errMsg = errorToConsole(err);
      console.error('MLB: MlbVrSpaceBrowser.scatterPloOnAfterDrawScene() error:\n' + errMsg);
      //throw err; // Do not throw to prevent disabling event handler
    } finally {
      this.context.restore();
    }
  }

  private treeDfOnMouseOverRowChanged(): void {
    console.debug('MLB: MlbVrSpaceBrowser.treeDfOnMouseOverRowChanged() ');
    if (!this.scatterPlot) return;
    try {
      this.mouseOverTreeIdx = this.treeDf.mouseOverRowIdx;
      //TODO: this.scatterPlot.invalidate();
      // hackish way to force invalidating scatterPlot
      const vp = this.scatterPlot.viewport;
      this.scatterPlot.zoom(vp.left, vp.top, vp.right, vp.bottom); // restore
    } catch (err: any) {
      const errMsg = errorToConsole(err);
      console.error('MLB: MlbVrSpaceBrowser.treeDfOnMouseOverRowChanged() error:\n' + errMsg);
      //throw err; // Do not throw to prevent disabling event handler
    }
  }

  private treeDfOnCurrentRowChanged(): void {
    console.debug('MLB: MlbVrSpaceBrowser.treeDfOnCurrentRowChanged() ');
    if (!this.scatterPlot) return;
    try {
      this.currentTreeIdx = this.treeDf.currentRowIdx;
      if (this.currentTreeIdx != -1) {
        const treeRoot: PointNodeType = this.treeList[this.currentTreeIdx];
        const treeLeafList: PointNodeType[] = this.th.getLeafList(treeRoot);
        const treeLeafWPointList: PointNodeType[] = treeLeafList.filter((n) => n.point != undefined);

        if (treeLeafWPointList.length >= 2) {
          const xList: number[] = treeLeafWPointList.map((n) => n.point!.x);
          const yList: number[] = treeLeafWPointList.map((n) => n.point!.y);
          const [minX, maxX, minY, maxY]: [number, number, number, number] =
            [Math.min(...xList), Math.max(...xList), Math.min(...yList), Math.max(...yList)];
          const vp = this.scatterPlot.viewport;
          const {zoomLeft, zoomRight, zoomTop, zoomBottom} = getZoomCoordinates(
            vp.width, vp.height, minX, minY, maxX, maxY);
          window.setTimeout(() => {
            this.scatterPlot.zoom(zoomLeft, zoomTop, zoomRight, zoomBottom);
          }, 0 /* next event cycle */);
        }
      }
    } catch (err: any) {
      const errMsg = errorToConsole(err);
      console.error('MLB: MlbVrSpaceBrowser.treeDfOnCurrentRowChanged() error:\n' + errMsg);
      //throw err; // Do not throw to prevent disabling event handler
    }
  }

  // -- Trees --

  /** Markup trees once to build trees' points */
  private static buildTreeList(
    th: ITreeHelper, treeCol: DG.Column, vrRows: { [vId: string]: number }, embedCols: { [embed: string]: DG.Column }
  ): PointNodeType[] {
    const embedColNameList = Object.values(EMBED_COL_NAMES);

    function markupNode(node: PointNodeType): void {
      if (isLeaf(node)) {
        node.vId = getVId(node.name);
        const vrRowIdx = vrRows[node.vId];
        node.point = !vrRowIdx ? undefined :
          new DG.Point(embedCols[embedColNameList[0]].get(vrRowIdx), embedCols[embedColNameList[1]].get(vrRowIdx));
      } else {
        for (const childNode of node.children!) markupNode(childNode as PointNodeType);

        const childrenWPoint: PointNodeType[] = (node.children! as PointNodeType[])
          .filter((cn) => cn.point != undefined);

        node.point = childrenWPoint.length == 0 ? undefined :
          new DG.Point(
            childrenWPoint.map((cn) => cn.point!.x)
              .reduce((apx, bpx) => apx + bpx) / childrenWPoint.length,
            childrenWPoint.map((cn) => cn.point!.y)
              .reduce((apy, bpy) => apy + bpy) / childrenWPoint.length);
      }
    }

    const treeRowCount: number = treeCol.length;
    const resTreeDataList: PointNodeType[] = new Array<PointNodeType>(treeRowCount);
    for (let treeRowI: number = 0; treeRowI < treeRowCount; treeRowI++) {
      const treeNewick: string = cleanMlbNewick(treeCol.get(treeRowI));
      const treeRoot: PointNodeType = parseNewick(treeNewick) as PointNodeType;
      markupNode(treeRoot);
      resTreeDataList[treeRowI] = treeRoot;
    }
    return resTreeDataList;
  }

  /** Get list of sequences (per mlbDf row) as concat result of subs derived of {@link regions}  */
  public static getSeqList(df: DG.DataFrame, seqCols: { [chain: string]: DG.Column }, regions: VdRegion[]): string[] {
    const chainPosDict: { [chain: string]: { col: DG.Column, startPos: number, endPos: number } } = Object.assign(
      {}, ...Object.keys(seqCols).map((chain) => {
        const cdr3Reg: VdRegion | undefined = regions.find((r) => r.chain == chain && r.name == 'CDR3')!;
        if (!cdr3Reg)
          throw new Error(`Region 'CDR3' not found for chain '${chain}'.`);
        const chainSeqCol: DG.Column = seqCols[chain];
        const posList: string[] = chainSeqCol.getTag('positionNames').split(', ');
        const k = 11;

        return ({
          [chain]: {
            col: chainSeqCol,
            startPos: posList.indexOf(cdr3Reg.positionStartName) - 1,
            endPos: posList.indexOf(cdr3Reg.positionEndName) + 1,
          }
        });
      }));

    const rowCount: number = df.rowCount;

    const resSeqList: string[] = new Array<string>(rowCount);
    for (let rowI: number = 0; rowI < rowCount; rowI++) {
      let resSeq: string = '';
      for (const chain in seqCols) { // iterate for chain name
        const chainSeqCol: DG.Column = seqCols[chain];
        const chainSeq: string = chainSeqCol.get(rowI);
        const chainPos: { col: DG.Column, startPos: number, endPos: number } = chainPosDict[chain];
        const chainSubSeq = chainSeq.substring(chainPos.startPos, chainPos.endPos);
        resSeq += chainSubSeq;
      }
      resSeqList[rowI] = resSeq;
    }

    return resSeqList;
  }

  // -- Static --

  static readonly emptyTreeDf: TreeDataFrame = TreeDataFrame.Empty;

  static readonly emptyMlbDf: MlbDataFrame = (() => {
    const res: MlbDataFrame = MlbDataFrame.Empty.clone();
    for (const embedColName of Object.values(EMBED_COL_NAMES)) res.columns.addNewFloat(embedColName);
    return res;
  })();
}

/** Updates screen coords and renders {@link node} on {@link ctx} */
function renderAndUpdateNode(
  ctx: CanvasRenderingContext2D, spv: ScatterPlotViewer, node: PointNodeType, styler: EdgeStyler
) {
  const nsp: DG.Point | undefined = node.screenPoint =
    node.point == undefined ? undefined : spv.worldToScreen(node.point.x, node.point.y);

  if (!isLeaf(node) && nsp != undefined) {
    for (const childNode of (node.children! as PointNodeType[])) {
      if (childNode.point != undefined) {
        renderAndUpdateNode(ctx, spv, childNode, styler);
        const cnsp: DG.Point = childNode.screenPoint!;

        ctx.strokeStyle = styler.getStrokeColor(); //toRgba(TreeDefaultPalette[TreeColorNames.Main]);
        ctx.lineWidth = styler.lineWidth;
        ctx.beginPath();
        ctx.moveTo(nsp.x, nsp.y);
        ctx.lineTo(cnsp.x, cnsp.y);
        ctx.stroke();
      }
    }
  }
}

function checkEmbedCols(embedColNameList: string[], df: DG.DataFrame): void {
  const missedList = embedColNameList
    .filter((embedColName) => !df.columns.names().includes(embedColName));
  if (missedList.length > 0)
    throw new Error(`Missed embed columns ${JSON.stringify(missedList)}.`);
}

class EdgeStyler {
  private _lineWidth: number;
  private _strokeColor: string;

  get lineWidth() { return this._lineWidth; }

  constructor(lineWidth: number, strokeColor: string) {
    this._lineWidth = lineWidth;
    this._strokeColor = strokeColor;
  }

  getStrokeColor(): string {
    return this._strokeColor;
  };
}

