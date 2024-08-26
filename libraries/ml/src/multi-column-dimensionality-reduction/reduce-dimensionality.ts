import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {KnownMetrics} from '../typed-metrics';
import {BYPASS_LARGE_DATA_WARNING, SHOW_SCATTERPLOT_PROGRESS} from '../functionEditors/consts';
import {Matrix, Options} from '@datagrok-libraries/utils/src/type-declarations';
import {IDBScanOptions, getDbscanWorker} from '@datagrok-libraries/math';
import {DIMENSIONALITY_REDUCER_TERMINATE_EVENT} from './consts';
import {PreprocessFunctionReturnType} from '../functionEditors/dimensionality-reduction-editor';
import {getNormalizedEmbeddings} from './embeddings-space';
import {DistanceAggregationMethod} from '../distance-matrix/types';
import {ITSNEOptions, IUMAPOptions} from './multi-column-dim-reducer';
import {DimReductionMethods} from './types';

export type DimRedUiOptions = {
  [BYPASS_LARGE_DATA_WARNING]?: boolean,
  [SHOW_SCATTERPLOT_PROGRESS]?: boolean,
  fastRowCount?: number,
  scatterPlotName?: string,
  tableView?: DG.TableView,
}

export function getEmbeddingColsNames(df: DG.DataFrame) {
  const axes = ['Embed_X', 'Embed_Y'];
  const colNameInd = df.columns.names().filter((it: string) => it.includes(axes[0])).length + 1;
  return axes.map((it) => `${it}_${colNameInd}`);
}

export function getEmbeddingViewerName(columns: DG.Column[], method: DimReductionMethods) {
  const colNames = columns.length > 4 ? `${columns.length} columns` : columns.map((it) => it.name).join(', ');
  return `${method} (${colNames})`;
}

export async function multiColReduceDimensionality(table: DG.DataFrame, columns: DG.Column[],
  method: DimReductionMethods, metrics: KnownMetrics[], weights: number[],
  preprocessingFunctions: (DG.Func | null | undefined)[],
  aggregationMethod: DistanceAggregationMethod, plotEmbeddings: boolean = true, clusterEmbeddings: boolean = false,
  dimRedOptions: (IUMAPOptions | ITSNEOptions) & Partial<IDBScanOptions> & {preprocessingFuncArgs: Options[]} &
    Options = {preprocessingFuncArgs: []},
  uiOptions: DimRedUiOptions = {}, postProcessingFunc: DG.Func | null = null, postProcFuncArgs: Options = {}
): Promise<DG.ScatterPlotViewer | undefined> {
  const scatterPlotProps = {
    showXAxis: false,
    showYAxis: false,
    showXSelector: false,
    showYSelector: false,
  };
  if (columns.length !== metrics.length || columns.length !== preprocessingFunctions.length ||
    columns.length !== weights.length || columns.length !== dimRedOptions.preprocessingFuncArgs.length) {
    throw new Error('columns, metrics and preprocessing functions, weights and function arguments' +
      'must have the same length');
  }

  const tv = plotEmbeddings ? (uiOptions.tableView ?? grok.shell.tableView(table.name) ?? grok.shell.addTableView(table)) : null;

  const doReduce = async () => {
    const pg = DG.TaskBarProgressIndicator.create(
      `Initializing ${uiOptions.scatterPlotName ?? 'dimensionality reduction'} ...`);
    let scatterPlot: DG.ScatterPlotViewer | undefined = undefined;
    try {
      const embedColsNames = getEmbeddingColsNames(table);
      function progressFunc(_nEpoch: number, epochsLength: number, embeddings: number[][]) {
        let embedXCol: DG.Column | null = null;
        let embedYCol: DG.Column | null = null;
        if (!table.columns.names().includes(embedColsNames[0])) {
          embedXCol = table.columns.add(DG.Column.float(embedColsNames[0], table.rowCount));
          embedYCol = table.columns.add(DG.Column.float(embedColsNames[1], table.rowCount));
          if (plotEmbeddings && !scatterPlot) {
            scatterPlot = tv!
              .scatterPlot({...scatterPlotProps, x: embedColsNames[0], y: embedColsNames[1],
                title: uiOptions.scatterPlotName ?? getEmbeddingViewerName(columns, method)});
          }
        } else {
          embedXCol = table.columns.byName(embedColsNames[0]);
          embedYCol = table.columns.byName(embedColsNames[1]);
        }

        if (uiOptions[SHOW_SCATTERPLOT_PROGRESS]) {
          scatterPlot?.root && ui.setUpdateIndicator(scatterPlot!.root, false);
          embedXCol.init((i) => embeddings[0] ? embeddings[0][i] : undefined);
          embedYCol.init((i) => embeddings[1] ? embeddings[1][i] : undefined);
        }
        const progress = (_nEpoch / epochsLength * 100);
        pg.update(progress,
          `Running ${uiOptions.scatterPlotName ?? 'dimensionality reduction'}... ${progress.toFixed(0)}%`);
      }
      async function getDimRed() {
        table.columns.add(DG.Column.float(embedColsNames[0], table.rowCount));
        table.columns.add(DG.Column.float(embedColsNames[1], table.rowCount));
        let resolveF: Function | null = null;
        if (plotEmbeddings) {
          scatterPlot = tv!
            .scatterPlot({...scatterPlotProps, x: embedColsNames[0], y: embedColsNames[1],
              title: uiOptions.scatterPlotName ?? getEmbeddingViewerName(columns, method)});
          ui.setUpdateIndicator(scatterPlot.root, true);
        }

        const sub = grok.events.onViewerClosed.subscribe((args) => {
          const v = args.args.viewer as unknown as DG.Viewer<any>;
          if (v?.getOptions()?.look?.title && scatterPlot?.getOptions()?.look?.title &&
            v?.getOptions()?.look?.title === scatterPlot?.getOptions()?.look?.title) {
            grok.events.fireCustomEvent(DIMENSIONALITY_REDUCER_TERMINATE_EVENT, {});
            sub.unsubscribe();
            resolveF?.();
            pg.close();
          }
        });

        const dimRedResPromise = new Promise<Matrix | undefined>(async (resolve, reject) => {
          try {
            resolveF = resolve;
            const encodedColEntries: PreprocessFunctionReturnType[] = [];
            for (let i = 0; i < preprocessingFunctions.length; ++i) {
              const pf = preprocessingFunctions[i];
              if (!dimRedOptions.distanceFnArgs)
                dimRedOptions.distanceFnArgs = [];
              if (pf) {
                const colInputName = pf.inputs[0].name;
                const metricInputName = pf.inputs[1].name;
                const {entries, options}: PreprocessFunctionReturnType =
                await pf.apply({[colInputName]: columns[i], [metricInputName]: metrics[i],
                  ...(dimRedOptions.preprocessingFuncArgs[i] ?? {})});
                encodedColEntries.push({entries, options});
                dimRedOptions.distanceFnArgs.push(options);
              } else {
                const entries = columns[i].toList();
                const options = {};
                encodedColEntries.push({entries, options});
                dimRedOptions.distanceFnArgs.push(options);
              }
            }
            const res = await getNormalizedEmbeddings(encodedColEntries.map((it) => it.entries), method,
              metrics, weights, aggregationMethod, dimRedOptions,
              uiOptions[BYPASS_LARGE_DATA_WARNING] ? undefined : progressFunc);
            resolve(res);
          } catch (e) {
            reject(e);
          }
        });
        const res = await dimRedResPromise;
        pg.close();
        sub.unsubscribe();
        return res;
      }
      const res = await getDimRed();

      if (clusterEmbeddings && res) {
        const clusterPg = DG.TaskBarProgressIndicator.create(`Clustering embeddings ...`);
        try {
          const clusterRes = await getDbscanWorker(res[0], res[1],
            dimRedOptions.dbScanEpsilon ?? 0.01, dimRedOptions.dbScanMinPts ?? 4);
          const clusterColName = table.columns.getUnusedName('Cluster (DBSCAN)');
          const clusterCol = table.columns.addNewString(clusterColName);
          clusterCol.init((i) => clusterRes[i].toString());
          if (scatterPlot)
            (scatterPlot as DG.ScatterPlotViewer).props.colorColumnName = clusterColName;
        } catch (e) {
          grok.shell.error('Clustering embeddings failed');
          console.error(e);
        } finally {
          clusterPg.close();
        }
      }
      if (res) {
        const embedXCol = table.columns.byName(embedColsNames[0]);
        const embedYCol = table.columns.byName(embedColsNames[1]);
        embedXCol.init((i) => res[0][i]);
        embedYCol.init((i) => res[1][i]);
        if (postProcessingFunc) {
          try {
            const col1InputName = postProcessingFunc.inputs[0].name;
            const col2InputName = postProcessingFunc.inputs[1].name;
            await postProcessingFunc
              .prepare({[col1InputName]: embedXCol, [col2InputName]: embedYCol, ...postProcFuncArgs})
              .call(true);
          } catch (e) {
            grok.shell.error('Post-processing failed');
            console.error(e);
          }
        }
        if (scatterPlot) {
          ui.setUpdateIndicator((scatterPlot as DG.ScatterPlotViewer).root, false);
          (scatterPlot as DG.ScatterPlotViewer).helpUrl = '/help/compute/sequence-space';
          return scatterPlot as DG.ScatterPlotViewer;
        }
      }
    } catch (e) {
      grok.shell.error('Dimensionality reduction failed');
      console.error(e);
      pg.close();
      if (scatterPlot)
        ui.setUpdateIndicator((scatterPlot as DG.ScatterPlotViewer).root, false);
    }
  };
  return new Promise<DG.ScatterPlotViewer | undefined>(async (resolve, reject) => {
    try {
      if (uiOptions.fastRowCount && table.rowCount > uiOptions.fastRowCount && !uiOptions[BYPASS_LARGE_DATA_WARNING]) {
        ui.dialog()
          .add(ui.divText('Analysis might take several minutes. Do you want to continue?'))
          .onOK(async () => {
            try {
              const res = await doReduce();
              resolve(res);
            } catch (e) {
              reject(e);
            }
          })
          .onCancel(() => resolve(undefined))
          .show();
      } else {
        const res = await doReduce();
        resolve(res);
      }
    } catch (e) {
      reject(e);
    }
  });
}
