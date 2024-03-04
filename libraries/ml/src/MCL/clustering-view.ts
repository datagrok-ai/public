import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {createMCLWorker} from './index';
import {ScatterPlotCurrentLineStyle, ScatterPlotLinesRenderer} from '@datagrok-libraries/utils/src/render-lines-on-sp';
import {KnownMetrics} from '../typed-metrics';
import {DistanceAggregationMethod} from '../distance-matrix/types';
import {PreprocessFunctionReturnType} from '../functionEditors/dimensionality-reduction-editor';
import {Options} from '@datagrok-libraries/utils/src/type-declarations';


export async function markovCluster(
  df: DG.DataFrame, cols: DG.Column[], metrics: KnownMetrics[],
  weights: number[], aggregationMethod: DistanceAggregationMethod, preprocessingFuncs: (DG.Func | null | undefined)[],
  preprocessingFuncArgs: any[], threshold: number = 80, maxIterations: number = 10
): Promise<void | DG.ScatterPlotViewer> {
  const scatterPlotProps = {
    showXAxis: false,
    showYAxis: false,
    showXSelector: false,
    showYSelector: false,
  };
  const tv = grok.shell.tableView(df.name) ?? grok.shell.addTableView(df);

  const sc = tv.scatterPlot({...scatterPlotProps, title: 'MCL'});

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
    threshold, weights, aggregationMethod, metrics, distanceFnArgs, maxIterations);

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
  const emberdXColName = df.columns.getUnusedName('EmbedX');
  const emberdYColName = df.columns.getUnusedName('EmbedY');
  const clustersCounter: {[_: number]: number} = {};
  res.clusters.forEach((c) => {
    if (!clustersCounter[c]) clustersCounter[c] = 0;
    clustersCounter[c]++;
  });
  const clusterCounterColName = df.columns.getUnusedName('Cluster size (MCL)');
  df.columns.addNewFloat(emberdXColName).init((i) => res.embedX[i]);
  df.columns.addNewFloat(emberdYColName).init((i) => res.embedY[i]);
  df.columns.addNewInt(clusterCounterColName).init((i) => clustersCounter[res.clusters[i]]);
  df.columns.addNewString(clusterColName).init((i) => res.clusters[i].toString());
  sc.props.xColumnName = emberdXColName;
  sc.props.yColumnName = emberdYColName;
  sc.props.colorColumnName = clusterColName;
  sc.props.markerDefaultSize = 5;
  terminateSub.unsubscribe();
  // const sc = tv.scatterPlot({x: emberdXColName, y: emberdYColName});
  // sc.props.colorColumnName = clusterColName;
  // sc.props.markerDefaultSize = 5;
  const _scLines = new ScatterPlotLinesRenderer(sc, emberdXColName, emberdYColName,
    {from: res.is as any, to: res.js as any, drawArrows: false, opacity: 0.3, skipMultiLineCalculation: true},
    ScatterPlotCurrentLineStyle.none);
  ui.setUpdateIndicator(sc.root, false);
  // sc.close();
  // const scLinesViewer = new ScatterPlotWithLines(sc, res.is, res.js, emberdXColName, emberdYColName);
  // tv.addViewer(scLinesViewer);
  return sc;
}

