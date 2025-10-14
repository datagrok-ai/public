import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {createMCLWorker} from './index';
import {ScatterPlotCurrentLineStyle, ScatterPlotLinesRenderer} from '@datagrok-libraries/utils/src/render-lines-on-sp';
import {KnownMetrics} from '../typed-metrics';
import {DistanceAggregationMethod} from '../distance-matrix/types';
import {PreprocessFunctionReturnType} from '../functionEditors/dimensionality-reduction-editor';
import {Options} from '@datagrok-libraries/utils/src/type-declarations';
import * as rxjs from 'rxjs';


export type MCLClusterViewerResult = {
  sc: DG.ScatterPlotViewer;
  embedXCol: DG.Column;
  embedYCol: DG.Column;
  clusterCol: DG.Column;
  clusterCounterCol: DG.Column;
  connectivityCol: DG.Column;
  i: ArrayLike<number>;
  j: ArrayLike<number>;
}

export async function markovCluster(
  df: DG.DataFrame, cols: DG.Column[], metrics: KnownMetrics[],
  weights: number[], aggregationMethod: DistanceAggregationMethod, preprocessingFuncs: (DG.Func | null | undefined)[],
  preprocessingFuncArgs: any[], threshold: number = 80, maxIterations: number = 10,
  useWebGPU: boolean = false, inflate: number = 2, minClusterSize: number = 5, scp?: DG.ScatterPlotViewer
): Promise<undefined | MCLClusterViewerResult> {
  const scatterPlotProps = {
    showXAxis: false,
    showYAxis: false,
    showXSelector: false,
    showYSelector: false,
  };

  let tv: DG.TableView | null = null;
  let sc: DG.ScatterPlotViewer | undefined = scp;
  if (!sc) {
    tv = grok.shell.tableView(df.name) ?? grok.shell.addTableView(df);
    sc = tv.scatterPlot({...scatterPlotProps, title: 'MCL'});
  }

  ui.setUpdateIndicator(sc.root, true);
  const distanceFnArgs: Options[] = [];
  const encodedColEntries: PreprocessFunctionReturnType[] = [];
  for (let i = 0; i < preprocessingFuncs.length; ++i) {
    const pf = preprocessingFuncs[i];
    if (pf) {
      const colInputName = pf.inputs[0].name;
      const metricInputName = pf.inputs[1].name;
      const {entries, options}: PreprocessFunctionReturnType =
                await pf.apply({[colInputName]: cols[i], [metricInputName]: metrics[i],
                  ...(preprocessingFuncArgs[i] ?? {})});
      encodedColEntries.push({entries, options});
      distanceFnArgs.push(options ?? {});
    } else {
      const entries = cols[i].toList();
      const options = {};
      encodedColEntries.push({entries, options});
      distanceFnArgs.push(options);
    }
  }

  const mclWorker = createMCLWorker(encodedColEntries.map((it) => it.entries),
    threshold, weights, aggregationMethod, metrics, distanceFnArgs, maxIterations, useWebGPU, inflate);

  const terminateSub = grok.events.onViewerClosed.subscribe((args) => {
    if (args.args.viewer?.props?.title === sc.props.title && sc.type === args.args?.viewer?.type) {
      terminateSub.unsubscribe();
      mclWorker.terminate();
    }
  });
  const res = await mclWorker.promise;
  if (!res)
    return;
  const clusterColName = df.columns.getUnusedName('Cluster (MCL)');
  const emberdXColName = df.columns.getUnusedName('EmbedX (MCL)');
  const emberdYColName = df.columns.getUnusedName('EmbedY (MCL)');
  const clustersCounter: {[_: number]: number} = {};
  res.clusters.forEach((c) => {
    if (!clustersCounter[c]) clustersCounter[c] = 0;
    clustersCounter[c]++;
  });
  const connectivity = new Uint32Array(res.embedX.length);
  for (let i = 0; i < res.is.length; i++) {
    connectivity[res.is[i]]++;
    connectivity[res.js[i]]++;
  }
  const clusterCounterColName = df.columns.getUnusedName('Cluster size (MCL)');
  const connectivityColName = df.columns.getUnusedName('Connectivity (MCL)');

  const embedXCol = df.columns.addNewFloat(emberdXColName);
  embedXCol.init((i) => res.embedX[i]);
  const embedYCol = df.columns.addNewFloat(emberdYColName);
  embedYCol.init((i) => res.embedY[i]);
  const clusterCol = df.columns.addNewString(clusterColName);
  clusterCol.init((i) => clustersCounter[res.clusters[i]] >= minClusterSize ? res.clusters[i].toString() : '-1');
  const catColorObj = {'-1': DG.Color.setAlpha(DG.Color.lightBlue, 100)};
  clusterCol.setTag(DG.TAGS.COLOR_CODING_CATEGORICAL,
    JSON.stringify(catColorObj));
  clusterCol.temp[DG.TAGS.COLOR_CODING_CATEGORICAL] = catColorObj;
  // clusterCol.setCategoryOrder(Array.from(new Set(res.clusters)).sort((a, b) => a - b).map((it) => it.toString()));

  const clusterCounterCol = df.columns.addNewInt(clusterCounterColName);
  clusterCounterCol.init((i) => clustersCounter[res.clusters[i]]);
  const connectivityCol = df.columns.addNewInt(connectivityColName);
  connectivityCol.init((i) => connectivity[i]);

  // df.columns.addNewFloat(emberdXColName).init((i) => res.embedX[i]);
  // df.columns.addNewFloat(emberdYColName).init((i) => res.embedY[i]);
  // df.columns.addNewInt(clusterCounterColName).init((i) => clustersCounter[res.clusters[i]]);
  // df.columns.addNewString(clusterColName).init((i) => res.clusters[i].toString());
  sc.props.xColumnName = emberdXColName;
  sc.props.yColumnName = emberdYColName;
  sc.props.colorColumnName = clusterColName;
  sc.props.markerDefaultSize = 6;
  terminateSub.unsubscribe();

  // limit cluster interconnectivity to 20 lines

  const maxInterClusterLinks = 20;
  // I know pushing to array is slow, but it's not a big deal here as its max 100k pushes in most cases
  const filteredIs: number[] = [];
  const filteredJs: number[] = [];
  const clusterInterconnectivity = new Map<number, Map<number, number>>();
  for (let i = 0; i < res.is.length; i++) {
    let from = res.clusters[res.is[i]];
    let to = res.clusters[res.js[i]];
    if (from === to) {
      filteredIs.push(res.is[i]);
      filteredJs.push(res.js[i]);
      continue;
    }
    if (from > to) {
      const tmp = from;
      from = to;
      to = tmp;
    }

    let fromMap = clusterInterconnectivity.get(from);
    if (!fromMap) {
      fromMap = new Map<number, number>();
      clusterInterconnectivity.set(from, fromMap);
    }
    let count = fromMap.get(to);
    if (!count)
      count = 0;

    if (count >= maxInterClusterLinks)
      continue;

    count++;
    fromMap.set(to, count);
    filteredIs.push(res.is[i]);
    filteredJs.push(res.js[i]);
  }


  // const _scLines = new ScatterPlotLinesRenderer(sc, emberdXColName, emberdYColName,
  //   {from: new Uint32Array(filteredIs) as any, to: new Uint32Array(filteredJs) as any,
  //     drawArrows: false, opacity: 0.3, skipMultiLineCalculation: true,
  //     skipShortLines: true, skipMouseOverDetection: true, shortLineThreshold: 6, width: 0.75, color: '128,128,128'},
  //   ScatterPlotCurrentLineStyle.none);

  //const _scLines = new SCLinesRenderer(sc, filteredIs, filteredJs, 6, 0.75, '128,128,128');

  ui.setUpdateIndicator(sc.root, false);
  // sc.close();
  // const scLinesViewer = new ScatterPlotWithLines(sc, res.is, res.js, emberdXColName, emberdYColName);
  // tv.addViewer(scLinesViewer);

  return {sc, embedXCol, embedYCol, clusterCol, clusterCounterCol, connectivityCol, i: filteredIs, j: filteredJs};
}


export class SCLinesRenderer {
  private renderFlag = false;
  renderSub: rxjs.Subscription;
  constructor(public sc: DG.ScatterPlotViewer, public from: ArrayLike<number>,
    public to: ArrayLike<number>, public shortLineThreshold: number, public width: number, public color: string) {
    this.renderSub = DG.debounce(sc.onAfterDrawScene, 200).subscribe(() => {
      if (this.renderFlag) {
        this.renderFlag = false;
        return;
      }
      this.renderFlag = true;
      const tempSub = sc.onBeforeDrawScene.subscribe((_) => {
        this.render();
        tempSub.unsubscribe();
      });
      setTimeout(() => {
        this.sc.invalidateCanvas();
        // this.renderFlag = false;
      });
    });
    sc.subs.push(this.renderSub);
  }

  render() {
    const xCol = this.sc.dataFrame.getCol(this.sc.props.xColumnName);
    const yCol = this.sc.dataFrame.getCol(this.sc.props.yColumnName);
    const filter = this.sc.filter;
    const positions = new Array(this.sc.dataFrame.rowCount).fill(null)
      .map((_, i) => !xCol.isNone(i) && !yCol.isNone(i) && filter.get(i) ? this.sc.pointToScreen(i) : null);
    const canvas = this.sc.canvas;
    const ctx = canvas.getContext('2d');
    if (!ctx)
      return;
    ctx.strokeStyle = `rgba(${this.color}, 0.3)`;
    ctx.lineWidth = this.width;
    const shortLineSquared = this.shortLineThreshold * this.shortLineThreshold;
    for (let i = 0; i < this.from.length; i++) {
      ctx.beginPath();
      const from = this.from[i];
      const to = this.to[i];
      if (positions[from] && positions[to]) {
        const fromPos = positions[from];
        const toPos = positions[to];
        const dx = toPos.x - fromPos.x;
        const dy = toPos.y - fromPos.y;
        const dist = dx * dx + dy * dy;
        if (dist < shortLineSquared)
          continue;
        ctx.moveTo(fromPos.x, fromPos.y);
        ctx.lineTo(toPos.x, toPos.y);
        ctx.stroke();
        ctx.closePath();
        // interestingly, doing stroke for each line turns to be 20 times faster than doing stroke for all lines at once
      }
    }
  }

  destroy() {
    this.renderSub.unsubscribe();
  }
}
