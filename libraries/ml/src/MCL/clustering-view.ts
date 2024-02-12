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


  const res = await createMCLWorker(encodedColEntries.map((it) => it.entries),
    threshold, weights, aggregationMethod, metrics, distanceFnArgs, maxIterations);
  const clusterColName = df.columns.getUnusedName('Cluster');
  const emberdXColName = df.columns.getUnusedName('EmbedX');
  const emberdYColName = df.columns.getUnusedName('EmbedY');
  const clustersCounter: {[_: number]: number} = {};
  res.clusters.forEach((c) => {
    if (!clustersCounter[c]) clustersCounter[c] = 0;
    clustersCounter[c]++;
  });
  const clusterCounterColName = df.columns.getUnusedName('Cluster size');
  df.columns.addNewFloat(emberdXColName).init((i) => res.embedX[i]);
  df.columns.addNewFloat(emberdYColName).init((i) => res.embedY[i]);
  df.columns.addNewInt(clusterCounterColName).init((i) => clustersCounter[res.clusters[i]]);
  df.columns.addNewString(clusterColName).init((i) => res.clusters[i].toString());
  const tv = grok.shell.tableView(df.name);
  if (tv) {
    const sc = tv.scatterPlot({x: emberdXColName, y: emberdYColName});
    sc.props.colorColumnName = clusterColName;
    sc.props.markerDefaultSize = 5;
    const _scLines = new ScatterPlotLinesRenderer(sc, emberdXColName, emberdYColName,
      {from: res.is as any, to: res.js as any, drawArrows: false, opacity: 0.3, skipMultiLineCalculation: true},
      ScatterPlotCurrentLineStyle.none);
    return sc;
  }
}
