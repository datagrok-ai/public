/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {debounceTime} from 'rxjs/operators';
import {getSimilarityFromDistance} from '../distance-metrics-methods';
import {Subject} from 'rxjs';
import {MmDistanceFunctionsNames} from '../macromolecule-distance-functions';
import '../../css/styles.css';
import {BitArrayMetrics} from '../typed-metrics/typed-metrics';
import {dmLinearIndex} from '../distance-matrix';
import {SparseMatrixResult, SparseMatrixService} from '../distance-matrix/sparse-matrix-service';
import {ILineSeries, MouseOverLineEvent, ScatterPlotCurrentLineStyle, ScatterPlotLinesRenderer} from '@datagrok-libraries/utils/src/render-lines-on-sp';
import {PreprocessFunctionReturnType} from '../functionEditors/dimensionality-reduction-editor';
import {getNormalizedEmbeddings} from '../multi-column-dimensionality-reduction/embeddings-space';
import {DimReductionMethods} from '../multi-column-dimensionality-reduction/types';
import {Matrix} from '@datagrok-libraries/utils/src/type-declarations';
import {MCLMethodName, createMCLWorker} from '../MCL';
import {multiColWebGPUSparseMatrix} from '@datagrok-libraries/math/src/webGPU/sparse-matrix/webGPU-sparse-matrix';
import {WEBGSLAGGREGATION} from '@datagrok-libraries/math/src/webGPU/multi-col-distances/webGPU-aggregation';
import wu from 'wu';

export let activityCliffsIdx = 0;

export const CLIFFS_DF_NAME = 'cliffsDf';

export const CLIFFS_COL_ENCODE_FN = 'cliffs_col_encode_fn';

export const enum TAGS {
  activityCliffs = '.activityCliffs',
}

export const enum TEMPS{
  cliffsDfGrid = '.cliffsDfGrid',
}

export interface ISequenceSpaceParams {
  seqCol: DG.Column,
  methodName: DimReductionMethods,
  similarityMetric: BitArrayMetrics | MmDistanceFunctionsNames,
  embedAxesNames: string[],
  options?: any
}

export interface ISequenceSpaceResult {
  distance?: Float32Array;
  coordinates: DG.ColumnList;
}

export type SequenceSpaceFunc = (params: ISequenceSpaceParams) => Promise<ISequenceSpaceResult>;

export interface ITooltipAndPanelParams {
  lineId: number,
  points: number[],
  df: DG.DataFrame,
  seqCol: DG.Column,
  activityCol: DG.Column,
  sali?: number
}

export interface IActivityCliffsMetrics {
  simVals: Float32Array;
  saliVals: Float32Array;
  n1: Uint32Array;
  n2: Uint32Array;
  cliffsBitSet: DG.BitSet;
}

interface ISaliLims {
  min: number;
  max: number;
}


export type ActivityCliffsLines = {
  lines: ILineSeries;
  linesDf: DG.DataFrame;
}

const filterCliffsSubj = new Subject<string>();
const LINES_DF_ACT_DIFF_COL_NAME = 'activity_difference';
const LINES_DF_SALI_COL_NAME = 'SALI_index';
const LINES_DF_SIM_COL_NAME = 'similarity';
const LINES_DF_LINE_IND_COL_NAME = 'line_index';
const LINES_DF_MOL_COLS_NAMES = ['1_molecule', '2_molecule'];
const CLIFFS_FILTER_APPLIED = 'filterCliffs';

// Searches for activity cliffs in a chemical dataset by selected cutoff
export async function getActivityCliffs(df: DG.DataFrame, seqCol: DG.Column,
  axesNames: string[], scatterTitle: string, activities: DG.Column, similarity: number,
  similarityMetric: BitArrayMetrics | MmDistanceFunctionsNames,
  methodName: DimReductionMethods, seqSpaceOptions: any, semType: string, tags: {[index: string]: string},
  encodingFunc: DG.Func, tooltipFunc: (params: ITooltipAndPanelParams) => HTMLElement,
  propertyPanelFunc: (params: ITooltipAndPanelParams) => HTMLElement,
  linesGridFunc?: (df: DG.DataFrame, pairColNames: string[]) => DG.Grid,
  cliffsDockRatio?: number, demo?: boolean) : Promise<DG.Viewer> {
  activityCliffsIdx++;
  const similarityLimit = similarity / 100;
  // eslint-disable-next-line prefer-const
  let acc: DG.Accordion;
  let clickedSp = false;

  const encodingFuncInputs = encodingFunc.inputs;
  const encodedColWithOptions: PreprocessFunctionReturnType =
    await encodingFunc.apply({[encodingFuncInputs[0].name]: seqCol,
      [encodingFuncInputs[1].name]: similarityMetric, ...(seqSpaceOptions.preprocessingFuncArgs ?? {})});

  let embeddingsMatrix: Matrix = [];
  if ((methodName as any) === MCLMethodName) {
    const mclRes = await createMCLWorker([encodedColWithOptions.entries], similarity, [1], 'MANHATTAN',
      [similarityMetric], [encodedColWithOptions.options??{}],
      seqSpaceOptions?.maxIterations ?? 5, seqSpaceOptions.useWebGPU ?? false).promise;
    df.columns.addNewInt(df.columns.getUnusedName('MCL Cluster')).init((i) => mclRes.clusters[i]);
    embeddingsMatrix = [mclRes.embedX, mclRes.embedY];
  } else {
    embeddingsMatrix = await getNormalizedEmbeddings([encodedColWithOptions.entries], methodName,
      [similarityMetric], [1], 'MANHATTAN', {...seqSpaceOptions, distanceFnArgs: [encodedColWithOptions.options??{}]});
  }
  if (embeddingsMatrix.length !== axesNames.length)
    throw new Error('Number of axes names should be equal to number of embedding dimensions');
  for (let i = 0; i < embeddingsMatrix.length; ++i)
    df.columns.addNewFloat(axesNames[i]).init((idx) => embeddingsMatrix[i][idx]);

  let sparseMatrixRes: SparseMatrixResult | null = null;
  if (seqSpaceOptions.useWebGPU) {
    try {
      sparseMatrixRes = await multiColWebGPUSparseMatrix(
        [encodedColWithOptions.entries], similarityLimit,
        [similarityMetric as any], WEBGSLAGGREGATION.MANHATTAN, [1],
        [encodedColWithOptions.options ?? {}]);
    } catch (e) {
      console.error(e);
    }
  }
  if (!sparseMatrixRes) {
    if (seqSpaceOptions.useWebGPU)
      console.error('WebGPU sparse matrix calculation failed, falling back to CPU implementation');
    sparseMatrixRes = await new SparseMatrixService()
      .calc(encodedColWithOptions.entries, similarityMetric, similarityLimit, encodedColWithOptions.options);
  }

  const cliffsMetrics: IActivityCliffsMetrics = await getSparseActivityCliffsMetrics(sparseMatrixRes, activities);

  const sali: DG.Column = getSaliCountCol(seqCol.length, `sali_${axesNames[0].substring(axesNames[0].lastIndexOf('_'))}`,
    cliffsMetrics.saliVals, cliffsMetrics.n1, cliffsMetrics.n2, axesNames, activities);
  df.columns.add(sali);

  const cliffsBitSet = cliffsMetrics.cliffsBitSet;

  const saliMinMax = getSaliMinMax(cliffsMetrics.saliVals);
  const saliOpacityCoef = 0.8 / (saliMinMax.max - saliMinMax.min);

  const view = grok.shell.getTableView(df.name);
  const sp = view.addViewer(DG.VIEWER.SCATTER_PLOT, {
    xColumnName: axesNames[0],
    yColumnName: axesNames[1],
    size: sali.name,
    color: activities.name,
    showXSelector: false,
    showYSelector: false,
    showSizeSelector: false,
    showColorSelector: false,
    markerMinSize: 5,
    markerMaxSize: 25,
    title: scatterTitle,
  }) as DG.ScatterPlotViewer;

  const linesRes = createLines(df, cliffsMetrics, seqCol, activities, semType, tags, saliMinMax, saliOpacityCoef);
  linesRes.lines.skipMultiLineCalculation = true;
  linesRes.linesDf.col(LINES_DF_SALI_COL_NAME)!.setTag('description', 'Structure−Activity Landscape Index (activity difference divided by 1 minus similarity)');
  //creating scatter plot lines renderer
  const spEditor = new ScatterPlotLinesRenderer(sp as DG.ScatterPlotViewer,
    axesNames[0], axesNames[1], linesRes.lines, ScatterPlotCurrentLineStyle.none);

  const linesDfGrid = linesGridFunc ?
    linesGridFunc(linesRes.linesDf, LINES_DF_MOL_COLS_NAMES).sort([LINES_DF_SALI_COL_NAME], [false]) :
    linesRes.linesDf.plot.grid().sort([LINES_DF_SALI_COL_NAME], [false]);

  if (linesDfGrid.col(LINES_DF_LINE_IND_COL_NAME))
    linesDfGrid.col(LINES_DF_LINE_IND_COL_NAME)!.visible = false;
  df.temp[TEMPS.cliffsDfGrid] = linesDfGrid;

  const listCliffsLink = ui.button(`${linesRes.linesDf.rowCount} cliffs`, () => {
    const viewerExists = wu(view.viewers).some((v) => v.dataFrame.name === `${CLIFFS_DF_NAME}${activityCliffsIdx}`);
    if (demo && !viewerExists) // Ensure the grid viewer is added only once if not already present in the demo app
      view.addViewer(linesDfGrid);
    view.dockManager.dock(linesDfGrid, 'down', null, 'Activity cliffs', cliffsDockRatio ?? 0.2);
  });
  listCliffsLink.classList.add('scatter_plot_link', 'cliffs_grid');

  /* in case several activity cliffs viewers are opened cliffs filtering can
  be applyed only to one of the viewers. When 'Show only cliffs' is switched on one of the viewers
  switch inputs on other viewers are disabled */
  const filterCliffsButton = ui.input.toggle(`Show only cliffs`, {value: false, onValueChanged: (value) => {
    if (value) {
      sp.dataFrame.setTag(CLIFFS_FILTER_APPLIED, axesNames[0]);
      df.filter.and(cliffsBitSet);
      filterCliffsSubj.next(axesNames[0]);
    } else {
      sp.dataFrame.setTag(CLIFFS_FILTER_APPLIED, '');
      df.filter.setAll(true);
      filterCliffsSubj.next('');
    }
  }});

  sp.root.prepend(ui.divH([
    listCliffsLink,
    filterCliffsButton.root
  ], 'cliffs_div'));

  filterCliffsSubj.subscribe((s: string) => {
    filterCliffsButton.enabled = s !== '' ? s !== axesNames[0] ? false : true : true;
  });

  //need to apply filtering by cliffs on each d4-before-draw-scene since redrawing scatter plot on zoom or move resets the filter
  let cliffsFilterApplied = false;
  sp.onEvent('d4-before-draw-scene')
    .subscribe((_: any) => {
      if (!cliffsFilterApplied) {
        if (filterCliffsButton.value) {
          setTimeout(() => { df.filter.and(cliffsBitSet); }, 100);
          cliffsFilterApplied = true;
        }
      } else {
        cliffsFilterApplied = false;
      }
    });


  //closing docked grid when viewer is closed
  const viewerClosedSub = grok.events.onViewerClosed.subscribe((v) => {
    //@ts-ignore
    if (v.args.viewer === sp) {
      view.dockManager.close(linesDfGrid.root);
      viewerClosedSub.unsubscribe();
      view.subs = view.subs.filter((sub) => sub !== viewerClosedSub);
    }
  });
  view.subs.push(viewerClosedSub);

  linesRes.linesDf.onCurrentCellChanged.subscribe(() => {
    for (let i = 0; i < linesRes.linesDf.rowCount; i++)
      linesRes.lines.widths![i] = i === linesRes.linesDf.currentRowIdx ? 3 : 1;
    spEditor.linesToRender = linesRes.lines;
    const molIdsArray = linesRes.linesDf.currentCol && linesRes.linesDf.currentCol.name === LINES_DF_MOL_COLS_NAMES[1] ?
      linesRes.lines.to : linesRes.lines.from;
    const lineIdx = linesRes.linesDf.currentRowIdx !== -1 ? linesRes.linesDf.currentRowIdx : null;
    sp.dataFrame.currentRowIdx = lineIdx ? molIdsArray[lineIdx] : -1;
    if (lineIdx !== null) {
      const currentLineId = linesRes.linesDf.currentRowIdx;
      spEditor.currentLineId = currentLineId;
      const {zoomLeft, zoomRight, zoomTop, zoomBottom} = getZoomCoordinates(
        sp.viewport.width,
        sp.viewport.height,
        sp.dataFrame.get(axesNames[0], linesRes.lines.from[currentLineId]),
        sp.dataFrame.get(axesNames[1], linesRes.lines.from[currentLineId]),
        sp.dataFrame.get(axesNames[0], linesRes.lines.to[currentLineId]),
        sp.dataFrame.get(axesNames[1], linesRes.lines.to[currentLineId])
      );
      sp.zoom(zoomLeft,
        zoomTop,
        zoomRight,
        zoomBottom);
      if (filterCliffsButton.value)
        df.filter.and(cliffsBitSet);
      else
      if (filterCliffsButton.enabled === true)
        df.filter.setAll(true);
      setTimeout(() => {
        updatePropertyPanel(df, acc, linesRes.lines.from[lineIdx], linesRes.lines.to[lineIdx], lineIdx,
          seqCol, activities, linesRes.linesDf.get(LINES_DF_SALI_COL_NAME, lineIdx), propertyPanelFunc);
        const order = sp.dataFrame.getSortedOrder(view.grid.sortByColumns, view.grid.sortTypes);
        view.grid.scrollToCell(seqCol.name, order.indexOf(sp.dataFrame.currentRowIdx));
      }, 1000);
    }
  });

  const updateParentDfSelectionAndLineColors = () => {
    const selection = DG.BitSet.create(df.rowCount);
    for (let i = 0; i < linesRes.linesDf.rowCount; i++) {
      const selected = linesRes.linesDf.selection.get(i);
      if (selected) {
        selection.set(linesRes.lines.from[i], true);
        selection.set(linesRes.lines.to[i], true);
      }
      linesRes.lines.colors![i] = selected ? '255,255,0' : '0,128,0';
    }
    df.selection.copyFrom(selection);
    spEditor.linesToRender = linesRes.lines;
  };

  const resetSelectionOnEsc = (dataframe: DG.DataFrame) => {
    dataframe.selection.setAll(false);
    for (let i = 0; i < linesRes.lines.colors!.length; i++)
      linesRes.lines.colors![i] = '0,128,0'; //reset color to default for all lines
    spEditor.linesToRender = linesRes.lines;
  };

  linesRes.linesDf.onSelectionChanged.subscribe((_e: any) => {
    setTimeout(() => updateParentDfSelectionAndLineColors(), 100);
  });

  df.onSelectionChanged.subscribe((e: any) => {
    if (df.selection.anyTrue === false && typeof e === 'number') { //catching event when initial df selection is reset by pushing Esc
      if (!clickedSp)
        resetSelectionOnEsc(linesDfGrid.dataFrame);
      else
        clickedSp = false;
    }
  });

  spEditor.lineClicked.subscribe((event: MouseOverLineEvent) => {
    clickedSp = true;
    spEditor.currentLineId = event.id;
    if (event.id !== -1) {
      const savedSelection = linesRes.linesDf.selection.clone();
      setTimeout(() => { //clicking on a scatter plot leads to selection reset, thus we need to set selection after a little timeout
        if (event.event.ctrlKey) {
          savedSelection.set(event.id, !savedSelection.get(event.id));
          linesRes.linesDf.selection.copyFrom(savedSelection);
        } else {
          if (linesRes.linesDf.currentRowIdx !== event.id) {
            linesRes.linesDf.currentRowIdx = event.id;
            df.currentRowIdx = linesRes.lines.from[event.id];
          }
          linesRes.linesDf.selection.copyFrom(savedSelection);
        }
        const order = linesRes.linesDf.getSortedOrder(linesDfGrid.sortByColumns, linesDfGrid.sortTypes);
        linesDfGrid.scrollToCell(LINES_DF_MOL_COLS_NAMES[0], order.indexOf(event.id));
      }, 500);
    }
  });


  spEditor.lineHover.pipe(debounceTime(500)).subscribe((event: MouseOverLineEvent) => {
    if (event.id !== -1 && df.mouseOverRowIdx === -1) {
      ui.tooltip.show(tooltipFunc({lineId: event.id,
        points: [linesRes.lines.from[event.id], linesRes.lines.to[event.id]], df: df, seqCol: seqCol,
        activityCol: activities}), event.x, event.y);
    }
  });

  sp.addProperty('similarityLimit', 'double', similarityLimit);
  acc = createPopertyPanel();
  return sp;
}


export async function getActivityCliffsEmbeddings(df: DG.DataFrame, seqCol: DG.Column,
  axesNames: string[], similarity: number, similarityMetric: BitArrayMetrics | MmDistanceFunctionsNames,
  methodName: DimReductionMethods, seqSpaceOptions: any, encodingFunc: DG.Func) : Promise<void> {
  // eslint-disable-next-line prefer-const

  const encodingFuncInputs = encodingFunc.inputs;
  const encodedColWithOptions: PreprocessFunctionReturnType =
    await encodingFunc.apply({[encodingFuncInputs[0].name]: seqCol,
      [encodingFuncInputs[1].name]: similarityMetric, ...(seqSpaceOptions.preprocessingFuncArgs ?? {})});

  let embeddingsMatrix: Matrix = [];
  if ((methodName as any) === MCLMethodName) {
    const mclRes = await createMCLWorker([encodedColWithOptions.entries], similarity, [1], 'MANHATTAN',
      [similarityMetric], [encodedColWithOptions.options??{}],
      seqSpaceOptions?.maxIterations ?? 5, seqSpaceOptions.useWebGPU ?? false).promise;
    df.columns.addNewInt(df.columns.getUnusedName('MCL Cluster')).init((i) => mclRes.clusters[i]);
    embeddingsMatrix = [mclRes.embedX, mclRes.embedY];
  } else {
    embeddingsMatrix = await getNormalizedEmbeddings([encodedColWithOptions.entries], methodName,
      [similarityMetric], [1], 'MANHATTAN', {...seqSpaceOptions, distanceFnArgs: [encodedColWithOptions.options??{}]});
  }
  if (embeddingsMatrix.length !== axesNames.length)
    throw new Error('Number of axes names should be equal to number of embedding dimensions');
  for (let i = 0; i < embeddingsMatrix.length; ++i)
    df.columns.addNewFloat(axesNames[i]).init((idx) => embeddingsMatrix[i][idx]);
}


export async function runActivityCliffs(sp: DG.ScatterPlotViewer, df: DG.DataFrame, seqCol: DG.Column,
  encodedColWithOptions: PreprocessFunctionReturnType, activities: DG.Column, axesNames: string[],
  similarity: number, similarityMetric: BitArrayMetrics | MmDistanceFunctionsNames,
  seqSpaceOptions: any, semType: string, tags: {[index: string]: string},
  tooltipFunc: (params: ITooltipAndPanelParams) => HTMLElement,
  propertyPanelFunc: (params: ITooltipAndPanelParams) => HTMLElement,
  linesGridFunc?: (df: DG.DataFrame, pairColNames: string[]) => DG.Grid,
  cliffsDockRatio?: number, demo?: boolean) : Promise<void> {
  activityCliffsIdx++;
  const similarityLimit = similarity / 100;
  // eslint-disable-next-line prefer-const
  let acc: DG.Accordion;
  let clickedSp = false;
  const view = grok.shell.getTableView(df.name);

  let sparseMatrixRes: SparseMatrixResult | null = null;
  if (seqSpaceOptions.useWebGPU) {
    try {
      sparseMatrixRes = await multiColWebGPUSparseMatrix(
        [encodedColWithOptions.entries], similarityLimit,
        [similarityMetric as any], WEBGSLAGGREGATION.MANHATTAN, [1],
        [encodedColWithOptions.options ?? {}]);
    } catch (e) {
      console.error(e);
    }
  }
  if (!sparseMatrixRes) {
    if (seqSpaceOptions.useWebGPU)
      console.error('WebGPU sparse matrix calculation failed, falling back to CPU implementation');
    sparseMatrixRes = await new SparseMatrixService()
      .calc(encodedColWithOptions.entries, similarityMetric, similarityLimit, encodedColWithOptions.options);
  }

  const cliffsMetrics: IActivityCliffsMetrics = await getSparseActivityCliffsMetrics(sparseMatrixRes, activities);
  const saliColName = `sali_${axesNames[0].substring(axesNames[0].lastIndexOf('_'))}`;
  if (!df.columns.names().includes(saliColName)) {
    const sali: DG.Column = getSaliCountCol(seqCol.length, saliColName, cliffsMetrics.saliVals, cliffsMetrics.n1,
      cliffsMetrics.n2, axesNames, activities);
    df.columns.add(sali);
    sp.setOptions({size: sali.name});
  }

  const saliMinMax = getSaliMinMax(cliffsMetrics.saliVals);
  const saliOpacityCoef = 0.8 / (saliMinMax.max - saliMinMax.min);

  const linesRes = createLines(df, cliffsMetrics, seqCol, activities, semType, tags, saliMinMax, saliOpacityCoef);
  linesRes.lines.skipMultiLineCalculation = true;
  linesRes.linesDf.col(LINES_DF_SALI_COL_NAME)!.setTag('description', 'Structure−Activity Landscape Index (activity difference divided by 1 minus similarity)');
  //creating scatter plot lines renderer
  const spEditor = new ScatterPlotLinesRenderer(sp as DG.ScatterPlotViewer,
    axesNames[0], axesNames[1], linesRes.lines, ScatterPlotCurrentLineStyle.none);

  const linesDfGrid = linesGridFunc ?
    linesGridFunc(linesRes.linesDf, LINES_DF_MOL_COLS_NAMES).sort([LINES_DF_SALI_COL_NAME], [false]) :
    linesRes.linesDf.plot.grid().sort([LINES_DF_SALI_COL_NAME], [false]);

  if (linesDfGrid.col(LINES_DF_LINE_IND_COL_NAME))
    linesDfGrid.col(LINES_DF_LINE_IND_COL_NAME)!.visible = false;
  df.temp[TEMPS.cliffsDfGrid] = linesDfGrid;

  const listCliffsLink = ui.button(`${linesRes.linesDf.rowCount} cliffs`, () => {
    const viewerExists = wu(view.viewers).some((v) => v.dataFrame.name === `${CLIFFS_DF_NAME}${activityCliffsIdx}`);
    if (demo && !viewerExists) // Ensure the grid viewer is added only once if not already present in the demo app
      view.addViewer(linesDfGrid);
    view.dockManager.dock(linesDfGrid, 'down', undefined, 'Activity cliffs', cliffsDockRatio ?? 0.2);
  });
  listCliffsLink.classList.add('scatter_plot_link', 'cliffs_grid');

  /* in case several activity cliffs viewers are opened cliffs filtering can
  be applyed only to one of the viewers. When 'Show only cliffs' is switched on one of the viewers
  switch inputs on other viewers are disabled */
  const filterCliffsButton = ui.input.toggle(`Show only cliffs`, {value: false, onValueChanged: (value) => {
    if (value) {
      sp.dataFrame.setTag(CLIFFS_FILTER_APPLIED, axesNames[0]);
      df.filter.and(cliffsMetrics.cliffsBitSet);
      filterCliffsSubj.next(axesNames[0]);
    } else {
      sp.dataFrame.setTag(CLIFFS_FILTER_APPLIED, '');
      df.filter.setAll(true);
      filterCliffsSubj.next('');
    }
  }});

  sp.root.prepend(ui.divH([
    listCliffsLink,
    filterCliffsButton.root
  ], 'cliffs_div'));

  filterCliffsSubj.subscribe((s: string) => {
    filterCliffsButton.enabled = s !== '' ? s !== axesNames[0] ? false : true : true;
  });

  //need to apply filtering by cliffs on each d4-before-draw-scene since redrawing scatter plot on zoom or move resets the filter
  let cliffsFilterApplied = false;
  sp.onEvent('d4-before-draw-scene')
    .subscribe((_: any) => {
      if (!cliffsFilterApplied) {
        if (filterCliffsButton.value) {
          setTimeout(() => { df.filter.and(cliffsMetrics.cliffsBitSet); }, 100);
          cliffsFilterApplied = true;
        }
      } else {
        cliffsFilterApplied = false;
      }
    });


  //closing docked grid when viewer is closed
  const viewerClosedSub = grok.events.onViewerClosed.subscribe((v) => {
    //@ts-ignore
    if (v.args.viewer === sp) {
      view.dockManager.close(linesDfGrid.root);
      viewerClosedSub.unsubscribe();
      view.subs = view.subs.filter((sub) => sub !== viewerClosedSub);
    }
  });
  view.subs.push(viewerClosedSub);

  linesRes.linesDf.onCurrentCellChanged.subscribe(() => {
    for (let i = 0; i < linesRes.linesDf.rowCount; i++)
      linesRes.lines.widths![i] = i === linesRes.linesDf.currentRowIdx ? 3 : 1;
    spEditor.linesToRender = linesRes.lines;
    const molIdsArray = linesRes.linesDf.currentCol && linesRes.linesDf.currentCol.name === LINES_DF_MOL_COLS_NAMES[1] ?
      linesRes.lines.to : linesRes.lines.from;
    const lineIdx = linesRes.linesDf.currentRowIdx !== -1 ? linesRes.linesDf.currentRowIdx : null;
    sp.dataFrame.currentRowIdx = lineIdx ? molIdsArray[lineIdx] : -1;
    if (lineIdx !== null) {
      const currentLineId = linesRes.linesDf.currentRowIdx;
      spEditor.currentLineId = currentLineId;
      const {zoomLeft, zoomRight, zoomTop, zoomBottom} = getZoomCoordinates(
        sp.viewport.width,
        sp.viewport.height,
        sp.dataFrame.get(axesNames[0], linesRes.lines.from[currentLineId]),
        sp.dataFrame.get(axesNames[1], linesRes.lines.from[currentLineId]),
        sp.dataFrame.get(axesNames[0], linesRes.lines.to[currentLineId]),
        sp.dataFrame.get(axesNames[1], linesRes.lines.to[currentLineId])
      );
      sp.zoom(zoomLeft,
        zoomTop,
        zoomRight,
        zoomBottom);
      if (filterCliffsButton.value)
        df.filter.and(cliffsMetrics.cliffsBitSet);
      else
      if (filterCliffsButton.enabled === true)
        df.filter.setAll(true);
      setTimeout(() => {
        updatePropertyPanel(df, acc, linesRes.lines.from[lineIdx], linesRes.lines.to[lineIdx], lineIdx,
          seqCol, activities, linesRes.linesDf.get(LINES_DF_SALI_COL_NAME, lineIdx), propertyPanelFunc);
        const order = sp.dataFrame.getSortedOrder(view.grid.sortByColumns, view.grid.sortTypes);
        view.grid.scrollToCell(seqCol.name, order.indexOf(sp.dataFrame.currentRowIdx));
      }, 1000);
    }
  });

  const updateParentDfSelectionAndLineColors = () => {
    const selection = DG.BitSet.create(df.rowCount);
    for (let i = 0; i < linesRes.linesDf.rowCount; i++) {
      const selected = linesRes.linesDf.selection.get(i);
      if (selected) {
        selection.set(linesRes.lines.from[i], true);
        selection.set(linesRes.lines.to[i], true);
      }
      linesRes.lines.colors![i] = selected ? '255,255,0' : '0,128,0';
    }
    df.selection.copyFrom(selection);
    spEditor.linesToRender = linesRes.lines;
  };

  const resetSelectionOnEsc = (dataframe: DG.DataFrame) => {
    dataframe.selection.setAll(false);
    for (let i = 0; i < linesRes.lines.colors!.length; i++)
      linesRes.lines.colors![i] = '0,128,0'; //reset color to default for all lines
    spEditor.linesToRender = linesRes.lines;
  };

  linesRes.linesDf.onSelectionChanged.subscribe((_e: any) => {
    setTimeout(() => updateParentDfSelectionAndLineColors(), 100);
  });

  df.onSelectionChanged.subscribe((e: any) => {
    if (df.selection.anyTrue === false && typeof e === 'number') { //catching event when initial df selection is reset by pushing Esc
      if (!clickedSp)
        resetSelectionOnEsc(linesDfGrid.dataFrame);
      else
        clickedSp = false;
    }
  });

  spEditor.lineClicked.subscribe((event: MouseOverLineEvent) => {
    clickedSp = true;
    spEditor.currentLineId = event.id;
    if (event.id !== -1) {
      const savedSelection = linesRes.linesDf.selection.clone();
      setTimeout(() => { //clicking on a scatter plot leads to selection reset, thus we need to set selection after a little timeout
        if (event.event.ctrlKey) {
          savedSelection.set(event.id, !savedSelection.get(event.id));
          linesRes.linesDf.selection.copyFrom(savedSelection);
        } else {
          if (linesRes.linesDf.currentRowIdx !== event.id) {
            linesRes.linesDf.currentRowIdx = event.id;
            df.currentRowIdx = linesRes.lines.from[event.id];
          }
          linesRes.linesDf.selection.copyFrom(savedSelection);
        }
        const order = linesRes.linesDf.getSortedOrder(linesDfGrid.sortByColumns, linesDfGrid.sortTypes);
        linesDfGrid.scrollToCell(LINES_DF_MOL_COLS_NAMES[0], order.indexOf(event.id));
      }, 500);
    }
  });


  spEditor.lineHover.pipe(debounceTime(500)).subscribe((event: MouseOverLineEvent) => {
    if (event.id !== -1 && df.mouseOverRowIdx === -1) {
      ui.tooltip.show(tooltipFunc({lineId: event.id,
        points: [linesRes.lines.from[event.id], linesRes.lines.to[event.id]], df: df, seqCol: seqCol,
        activityCol: activities}), event.x, event.y);
    }
  });

  sp.addProperty('similarityLimit', 'double', similarityLimit);
  acc = createPopertyPanel();
}


async function getSparseActivityCliffsMetrics(
  sparseMatrix: SparseMatrixResult, activities: DG.Column
): Promise<IActivityCliffsMetrics> {
  const saliVals = sparseMatrix.distance.map((d, idx) => {
    const diff = Math.abs(activities.get(sparseMatrix.i[idx]) - activities.get(sparseMatrix.j[idx]));
    return d != 0 ? diff / d : Infinity;
  });
  const simVals = sparseMatrix.distance.map((d) => 1 - d);
  const n1 = sparseMatrix.i;
  const n2 = sparseMatrix.j;
  const bitset = DG.BitSet.create(activities.length);
  sparseMatrix.distance.forEach((_, idx) => {
    bitset.set(sparseMatrix.i[idx], true);
    bitset.set(sparseMatrix.j[idx], true);
  });
  return {simVals: simVals, saliVals: saliVals, n1: n1, n2: n2, cliffsBitSet: bitset} as IActivityCliffsMetrics;
}


function getSaliMinMax(saliVals: Float32Array): ISaliLims {
  const saliValsWithoutInfinity = saliVals.filter((it) => it !== Infinity);
  const saliMin = saliValsWithoutInfinity.reduce((acc, val) => Math.min(acc, val), Number.MAX_VALUE); //Math.min(...saliValsWithoutInfinity);
  const saliMax = saliValsWithoutInfinity.reduce((acc, val) => Math.max(acc, val), saliMin); //Math.max(...saliValsWithoutInfinity);
  return {max: saliMax, min: saliMin};
}

function getSaliCountCol(length: number, saliColName: string, saliVals: Float32Array, n1: Uint32Array, n2: Uint32Array,
  axesNames: string[], activities: DG.Column): DG.Column {
  const saliCount = new Array(length).fill(0);

  for (let i = 0; i != n1.length; ++i) {
    if (saliVals[i] != Infinity) {
      if (activities.get(n1[i]) > activities.get(n2[i]))
        saliCount[n1[i]] += saliVals[i];
      else
        saliCount[n2[i]] += saliVals[i];
    }
  }

  return DG.Column.fromList('double', saliColName, saliCount);
}


function createPopertyPanel(): DG.Accordion {
  const acc = ui.accordion();
  const accIcon = ui.element('i');
  accIcon.className = 'grok-icon svg-icon svg-view-layout';
  acc.addTitle(ui.span([accIcon, ui.label(`Activity cliffs`)]));
  acc.addPane('Cliff Details', () => ui.divText('Cliff has not been selected'), true);
  grok.shell.o = acc.root;
  return acc;
}

function updatePropertyPanel(df: DG.DataFrame, acc: DG.Accordion, pointFrom: number, pointTo: number, lineId: number,
  seqCol: DG.Column, activities: DG.Column, sali: number, propPanelFunc: (params: ITooltipAndPanelParams) => HTMLElement) {
  const panel = acc.getPane('Cliff Details');
  ui.empty(panel.root);
  const panelElement = propPanelFunc({points: [pointFrom, pointTo], lineId, df: df, seqCol: seqCol,
    activityCol: activities, sali: sali});
  panel.root.append(panelElement);
  setTimeout(() => {
    grok.shell.o = acc.root;
  }, 500);
}

function getZoomCoordinates(W0: number, H0: number, x1: number, y1: number, x2: number, y2: number) {
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


function createLines(df: DG.DataFrame, params: IActivityCliffsMetrics, seq: DG.Column, activities: DG.Column, semType: string,
  tags: {[index: string]: string}, saliMinMax: ISaliLims, saliOpacityCoef: number) : ActivityCliffsLines {
  const lines: ILineSeries = {
    from: new Uint32Array(params.n1.length),
    to: new Uint32Array(params.n1.length),
    opacities: new Float32Array(params.n1.length),
    colors: new Array<string>(params.n1.length),
    widths: new Float32Array(params.n1.length)
  };
  for (let i = 0; i < params.n1.length; i++) {
    lines.from[i] = params.n1[i];
    lines.to[i] = params.n2[i];
    lines.opacities![i] = params.saliVals[i] === Infinity ? 1 : 0.2 + (params.saliVals[i] - saliMinMax.min) * saliOpacityCoef;
    lines.colors![i] = df.selection.get(lines.from[i]) && df.selection.get(lines.to[i]) ? '255,255,0' : '0,128,0';
    lines.widths![i] = 1;
  }
  const linesDf = DG.DataFrame.create(lines.from.length);
  LINES_DF_MOL_COLS_NAMES.forEach((it, idx) => {
    linesDf.columns.addNewString(it).init((i: number) => seq.get(idx === 0 ? lines.from[i] : lines.to[i]));
    setTags(linesDf.col(it)!, tags);
    linesDf.col(it)!.semType = semType;
  });
  linesDf.columns.addNewFloat(LINES_DF_ACT_DIFF_COL_NAME)
    .init((i: number) => Math.abs(activities.get(lines.from[i]) - activities.get(lines.to[i])));
  linesDf.columns.addNewInt(LINES_DF_LINE_IND_COL_NAME).init((i: number) => i);
  linesDf.columns.addNewFloat(LINES_DF_SALI_COL_NAME).init((i: number) => params.saliVals[i]);
  linesDf.columns.addNewFloat(LINES_DF_SIM_COL_NAME).init((i: number) => params.simVals[i]);
  linesDf.name = `${CLIFFS_DF_NAME}${activityCliffsIdx}`;
  return {lines, linesDf};
}

function setTags(col: DG.Column, tags: {[index: string]: string}) {
  Object.keys(tags).forEach((tag) => {
    col.tags[tag] = tags[tag];
  });
}

export async function getSimilaritiesMatrix( dim: number, seqCol: DG.Column, dfSeq: DG.DataFrame, simArr: DG.Column[],
  simFunc: (col: DG.Column, mol: string) => Promise<DG.Column | null>): Promise<DG.Column[]> {
  for (let i = 0; i != dim - 1; ++i) {
    const mol = seqCol.get(i);
    dfSeq.rows.removeAt(0, 1, false);
    simArr[i] = (await simFunc(dfSeq.col('seq')!, mol))!;
  }
  return simArr;
}

export function getSimilaritiesFromDistances(dim: number, distances: Float32Array,
  simArr: (DG.Column | null)[]): (DG.Column | null)[] {
  const linearIndexFunc = dmLinearIndex(dim);
  for (let i = 0; i < dim - 1; ++i) {
    const similarityArr = new Float32Array(dim - i - 1).fill(0);
    for (let j = i + 1; j < dim; ++j) {
      const index = linearIndexFunc(i, j);
      similarityArr[j - i - 1] = distances[index] === DG.FLOAT_NULL ? 0 : getSimilarityFromDistance(distances[index]);
    }

    simArr[i] = DG.Column.fromFloat32Array('similarity', similarityArr);
  }
  return simArr;
}

export function getCliffsBooleanCol(df: DG.DataFrame, ids: Set<number>): DG.Column {
  const colname = 'containsCliff';
  const colNameInd = df.columns.names().filter((it: string) => it.includes(colname)).length + 1;
  const newColName = `${colname}_${colNameInd}`;
  return DG.Column.bool(newColName, df.rowCount).init((i) => ids.has(i));
}

export function getCliffsBitset(df: DG.DataFrame, ids: Set<number>): DG.BitSet {
  const bitset = DG.BitSet.create(df.rowCount);
  for (let i = 0; i < df.rowCount; i++)
    bitset.set(i, ids.has(i));
  return bitset;
}
