import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Matrix} from '@datagrok-libraries/utils/src/type-declarations';
import {getSimilarityFromDistance} from '../distance-metrics-methods';
import {removeEmptyStringRows} from '@datagrok-libraries/utils/src/dataframe-utils';
import {Subject} from 'rxjs';
import '../../css/styles.css';

export let activityCliffsIdx = 0;

export const CLIFFS_DF_NAME = 'cliffsDf';

export const enum TAGS {
  activityCliffs = '.activityCliffs',
}

export const enum TEMPS{
  cliffsDfGrid = '.cliffsDfGrid',
}

export interface ILine {
  id: number;
  mols: number[];
  selected: boolean;
  a: number[]; // [x, y]
  b: number[]; // [x, y]
}

interface IRenderedLines {
  lines: ILine[];
  linesDf: DG.DataFrame;
}

export interface ISequenceSpaceParams {
  seqCol: DG.Column,
  methodName: string,
  similarityMetric: string,
  embedAxesNames: string[],
  options?: any
}

export interface ISequenceSpaceResult {
  distance: Matrix;
  coordinates: DG.ColumnList;
}

export interface ITooltipAndPanelParams {
  cashedData: any,
  line: ILine,
  df: DG.DataFrame,
  seqCol: DG.Column,
  activityCol: DG.Column,
  sali?: number
}

export interface IActivityCliffsMetrics {
  simVals: number[];
  saliVals: number[];
  n1: number[];
  n2: number[];
  cliffsMolIds: Set<number>;
}

interface ISaliLims {
  min: number;
  max: number;
}

const filterCliffsSubj = new Subject<string>();
const nonNormalizedDistances = ['Levenshtein'];
const LINES_DF_ACT_DIFF_COL_NAME = 'act_diff';
const LINES_DF_SALI_COL_NAME = 'sali';
const LINES_DF_SIM_COL_NAME = 'sim';
const LINES_DF_LINE_IND_COL_NAME = 'line_index';
const LINES_DF_MOL_COLS_NAMES = ['1_seq', '2_seq'];
const CLIFFS_FILTER_APPLIED = 'filterCliffs';

// Searches for activity cliffs in a chemical dataset by selected cutoff
export async function getActivityCliffs(df: DG.DataFrame, seqCol: DG.Column, encodedCol: DG.Column | null,
  axesNames: string[], scatterTitle: string, activities: DG.Column, similarity: number, similarityMetric: string,
  methodName: string, semType: string, tags: {[index: string]: string},
  seqSpaceFunc: (params: ISequenceSpaceParams) => Promise<ISequenceSpaceResult>,
  simMatrixFunc: (dim: number, seqCol: DG.Column, df: DG.DataFrame, colName: string,
    simArr: DG.Column[]) => Promise<DG.Column[]>,
  tooltipFunc: (params: ITooltipAndPanelParams) => HTMLElement,
  propertyPanelFunc: (params: ITooltipAndPanelParams) => HTMLElement,
  linesGridFunc?: (df: DG.DataFrame, pairColNames: string[]) => DG.Grid,
  seqSpaceOptions?: any, cliffsDockRatio?: number) : Promise<DG.Viewer> {

  activityCliffsIdx++;
  const similarityLimit = similarity / 100;
  const dimensionalityReduceCol = encodedCol ?? seqCol;
  let zoom = false;
  let ignoreSelectionChange = false;
  const cashedLinesData: any = {};
  let acc: DG.Accordion;
  let timer: NodeJS.Timeout;

  const seqSpaceParams = {
    seqCol: dimensionalityReduceCol,
    methodName: methodName,
    similarityMetric: similarityMetric,
    embedAxesNames: axesNames,
    options: seqSpaceOptions,
  };

  const { distance, coordinates } = await seqSpaceFunc(seqSpaceParams);
  for (const col of coordinates)
    df.columns.add(col);

  const simArr = await createSimilaritiesMatrix(dimensionalityReduceCol, distance,
    !distance || nonNormalizedDistances.includes(similarityMetric), simMatrixFunc);

  const cliffsMetrics: IActivityCliffsMetrics = getActivityCliffsMetrics(simArr, similarityLimit, activities);

  const sali: DG.Column = getSaliCountCol(dimensionalityReduceCol.length, cliffsMetrics.saliVals, cliffsMetrics.n1,
    cliffsMetrics.n2, axesNames, activities);
  df.columns.add(sali);

  const cliffCol: DG.Column = getCliffsBooleanCol(df, cliffsMetrics.cliffsMolIds);
  df.columns.add(cliffCol);

  const saliMinMax = getSaliMinMax(cliffsMetrics.saliVals);
  const saliOpacityCoef = 0.8/(saliMinMax.max - saliMinMax.min);;

  const view = grok.shell.getTableView(df.name);
  view.grid.columns.byName(cliffCol.name)!.visible = false;
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

  const canvas = (sp.getInfo() as any)['canvas'];
  const linesRes = createLines(cliffsMetrics, seqCol, activities, semType, tags);

  const linesDfGrid = linesGridFunc ?
    linesGridFunc(linesRes.linesDf, LINES_DF_MOL_COLS_NAMES).sort([LINES_DF_SALI_COL_NAME], [false]) :
    linesRes.linesDf.plot.grid().sort([LINES_DF_SALI_COL_NAME], [false]);
  df.temp[TEMPS.cliffsDfGrid] = linesDfGrid;

  const listCliffsLink = ui.button(`${linesRes.linesDf.rowCount} cliffs`, () => {
    view.dockManager.dock(linesDfGrid, 'down', null, 'Activity cliffs', cliffsDockRatio ?? 0.2);
  });
  listCliffsLink.classList.add('scatter_plot_link', 'cliffs_grid');
  sp.root.append(listCliffsLink);

  /* in case several activity cliffs viewers are opened cliffs filtering can
  be applyed only to one of the viewers. When 'Show only cliffs' is switched on one of the viewers
  switch inputs on other viewers are disabled */
  const filterCliffsButton = ui.switchInput(`Show only cliffs`, false, () => {
    if (filterCliffsButton.value) {
      sp.dataFrame.setTag(CLIFFS_FILTER_APPLIED, cliffCol.name);
      df.filter.copyFrom(createCliffsOnlyFilter(df, cliffCol.name));
      filterCliffsSubj.next(cliffCol.name);
    } else {
      sp.dataFrame.setTag(CLIFFS_FILTER_APPLIED, '');
      df.filter.setAll(true, true);
      filterCliffsSubj.next('');
    }
  });
  filterCliffsButton.root.classList.add('scatter_plot_link', 'show_only_cliffs');
  sp.root.append(filterCliffsButton.root);

  filterCliffsSubj.subscribe((s: string) => {
    if (s !== '') {
      if (s !== cliffCol.name)
        filterCliffsButton.enabled = false;
    } else {
      filterCliffsButton.enabled = true;
    }
  })

  //closing docked grid when viewer is closed
  const viewerClosedSub = grok.events.onViewerClosed.subscribe((v) => {
    //@ts-ignore
    if(v.args.viewer === sp) {
      view.dockManager.close(linesDfGrid.root);
      viewerClosedSub.unsubscribe();
      view.subs = view.subs.filter((sub) => sub !== viewerClosedSub);
    }
  })
  view.subs.push(viewerClosedSub);


  linesRes.linesDf.onCurrentCellChanged.subscribe(() => {
    zoom = true;
    const currentMolIdx = linesRes.linesDf.currentCol && linesRes.linesDf.currentCol.name === LINES_DF_MOL_COLS_NAMES[1] ? 1 : 0;
    const line = linesRes.linesDf.currentRowIdx !== -1 ? linesRes.lines[linesRes.linesDf.currentRowIdx] : null;
    sp.dataFrame.currentRowIdx = line ? line.mols[currentMolIdx] : -1;
    sp.dataFrame.filter.set(0, !linesRes.lines[0].selected); //calling filter to force scatter plot re-rendering
    sp.dataFrame.filter.set(0, linesRes.lines[0].selected);
    if (line) {
      setTimeout(() => {
        updatePropertyPanel(df, acc, cashedLinesData, line, seqCol, activities,
          linesRes.linesDf.get(LINES_DF_SALI_COL_NAME, line.id), propertyPanelFunc);
        const order = sp.dataFrame.getSortedOrder(view.grid.sortByColumns, view.grid.sortTypes);
        view.grid.scrollToCell(seqCol.name, order.indexOf(sp.dataFrame.currentRowIdx));
      }, 1000);
    }
  });

  linesRes.linesDf.onSelectionChanged.subscribe((_: any) => {
    if (linesRes.linesDf.selection.anyTrue === false) { //case when selection is reset by pushing Esc
      linesRes.lines.forEach((l) => { l.selected = false; });
    } else {
      if (linesRes.linesDf.mouseOverRowIdx !== -1) {
        const line = linesRes.lines[linesRes.linesDf.mouseOverRowIdx];
        line.selected = !line.selected;
      }
    }
    setTimeout(() => {
      const selection = DG.BitSet.create(df.rowCount);
      linesRes.lines.forEach((l) => {
        if (l.selected) {
          l.mols.forEach((m) => {
            selection.set(m, l.selected, true);
          });
        }
      });
      df.selection.copyFrom(selection);
      view.grid.invalidate();
    }, 300); //timeout is required because of resetting selection by clicking on scatter plot. First selection is reset and after it is set again using by this method
  });

  df.onSelectionChanged.subscribe((e: any) => {
    if (ignoreSelectionChange) {
      ignoreSelectionChange = false;
    } else {
      if (df.selection.anyTrue === false && typeof e === 'number') { //catching event when initial df selection is reset by pushing Esc
        linesRes.lines.forEach((l) => {l.selected = false; });
        linesDfGrid.dataFrame.selection.setAll(false, false);
        linesDfGrid.invalidate();
      }
    }
  });

  canvas.addEventListener('mousemove', function(event: MouseEvent) {
    clearTimeout(timer);
    timer = global.setTimeout(function() {
      const line = checkCursorOnLine(event, canvas, linesRes.lines);
      if (line && df.mouseOverRowIdx === -1) {
        ui.tooltip.show(tooltipFunc({cashedData: cashedLinesData, line: line, df: df, seqCol: seqCol,
          activityCol: activities}), event.clientX, event.clientY);
      }
    }, 500);
  });

  canvas.addEventListener('mousedown', function(event: MouseEvent) {
    ignoreSelectionChange = true; //clicking on a scatter plot leads to selection reset. Variable is required to prevent from resetting (is used in df.onSelectionChanged.subscribe)
    const line = checkCursorOnLine(event, canvas, linesRes.lines);
    if (line && df.mouseOverRowIdx === -1) {
      if (event.ctrlKey) {
        line.selected = !line.selected;
        linesRes.linesDf.selection.set(line.id, line.selected);
      } else {
        if (linesRes.linesDf.currentRowIdx !== line.id) {
          linesRes.linesDf.currentRowIdx = line.id;
          df.currentRowIdx = line.mols[0];
          df.filter.set(0, !linesRes.lines[0].selected); //calling filter to force scatter plot re-rendering
          df.filter.set(0, linesRes.lines[0].selected);
        }
      }
      const order = linesRes.linesDf.getSortedOrder(linesDfGrid.sortByColumns, linesDfGrid.sortTypes);
      linesDfGrid.scrollToCell(LINES_DF_MOL_COLS_NAMES[0], order.indexOf(line.id));
    }
  });

  sp.onEvent('d4-before-draw-scene')
    .subscribe((_: any) => {
      const lines = renderLines(sp,
        axesNames[0], axesNames[1], linesRes, cliffsMetrics.saliVals, saliOpacityCoef, saliMinMax.min);
      if (zoom) {
        const currentLine = lines[linesRes.linesDf.currentRowIdx];
        setTimeout(()=> {
          const {zoomLeft, zoomRight, zoomTop, zoomBottom} = getZoomCoordinates(
            sp.viewport.width,
            sp.viewport.height,
            sp.dataFrame.get(axesNames[0], currentLine.mols[0]),
            sp.dataFrame.get(axesNames[1], currentLine.mols[0]),
            sp.dataFrame.get(axesNames[0], currentLine.mols[1]),
            sp.dataFrame.get(axesNames[1], currentLine.mols[1]),
          );
          sp.zoom(zoomLeft,
            zoomTop,
            zoomRight,
            zoomBottom);
        }, 300);
        zoom = false;
      }
      if (filterCliffsButton.value) {
        df.filter.copyFrom(createCliffsOnlyFilter(df, cliffCol.name));
      }
      else
        if (filterCliffsButton.enabled === true)
          df.filter.setAll(true, false);
    });

  sp.addProperty('similarityLimit', 'double', similarityLimit);
  acc = createPopertyPanel();
  return sp;
}

function createCliffsOnlyFilter(df: DG.DataFrame, colName: string): DG.BitSet {
  const filter = DG.BitSet.create(df.rowCount);
  const raw = df.col(colName)!.getRawData();
  for (let i = 0; i < raw.length; i++) {
    filter.set(i, !!raw[i], false);
  }
  return filter;
}

async function createSimilaritiesMatrix(col: DG.Column, distance: Matrix, countFromDistance: boolean,
  simMatrixFunc: (dim: number, seqCol: DG.Column, df: DG.DataFrame, colName: string, simArr: DG.Column[]) => Promise<DG.Column[]>):Promise<DG.Column[]> {
  const cats = col.categories;
  const raw = col.getRawData();
  const newCol = DG.Column.string('seq', col.length).init((i) => cats[raw[i]]);
  const dfSeq = DG.DataFrame.fromColumns([newCol]);
  const dim = col.length;
  let simArr: DG.Column[] = Array(dim - 1);

  if (countFromDistance) {
    simArr = await simMatrixFunc(dim, col, dfSeq, 'seq', simArr);
  } else {
    getSimilaritiesFromDistances(dim, distance, simArr);
  }

  return simArr;
}

function getActivityCliffsMetrics(simArr: DG.Column[], similarityLimit: number, activities: DG.Column): IActivityCliffsMetrics {
  const simVals: number[] = [];
  const saliVals: number[] = [];
  const n1: number[] = [];
  const n2: number[] = [];
  const cliffsMolIds = new Set<number>();

  for (let i = 0; i != simArr.length; ++i) {
    for (let j = 0; j != simArr.length - i; ++j) {
      const sim: number = simArr[i] ? simArr[i].get(j) : 0;

      if (sim >= similarityLimit) {
        n1.push(i);
        n2.push(i + j + 1);
        cliffsMolIds.add(i);
        cliffsMolIds.add(i + j + 1);
        simVals.push(sim);
        const diff = Math.abs(activities.get(i) - activities.get(i + j + 1));
        if (sim != 1) {
          saliVals.push(diff / (1 - sim));
        } else {
          saliVals.push(Infinity);
        }
      }
    }
  }

  return {simVals: simVals, saliVals: saliVals, n1: n1, n2: n2, cliffsMolIds: cliffsMolIds};
}

function getSaliMinMax(saliVals: number[]): ISaliLims {
  const saliValsWithoutInfinity = saliVals.filter((it) => it !== Infinity);
  const saliMin = Math.min(...saliValsWithoutInfinity);
  const saliMax = Math.max(...saliValsWithoutInfinity);
  return {max: saliMax, min: saliMin};
}

function getSaliCountCol(length: number, saliVals: number[], n1: number[], n2: number[],
  axesNames: string[], activities: DG.Column): DG.Column {
  const saliCount = new Array(length).fill(0);

  for (let i = 0; i != n1.length; ++i) {
    if (saliVals[i] != Infinity) {
      if (activities.get(n1[i]) > activities.get(n2[i])) {
        saliCount[n1[i]] += saliVals[i];
      } else {
        saliCount[n2[i]] += saliVals[i];
      }
    }
  }

  return DG.Column.fromList('double', `sali_${axesNames[0].substring(axesNames[0].lastIndexOf('_'))}`, saliCount);
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

function updatePropertyPanel(df: DG.DataFrame, acc: DG.Accordion, cashedData: any, line: ILine, seqCol: DG.Column,
  activities: DG.Column, sali: number, propPanelFunc: (params: ITooltipAndPanelParams) => HTMLElement) {
  const panel = acc.getPane('Cliff Details');
  ui.empty(panel.root);
  const panelElement = propPanelFunc({cashedData: cashedData, line: line, df: df, seqCol: seqCol,
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

function checkCursorOnLine(event: any, canvas: any, lines: ILine[]): ILine | null {
  const rect = canvas.getBoundingClientRect();
  const x = event.clientX - rect.left;
  const y = event.clientY - rect.top;
  let closestLine = null;
  let minDist = 0;
  for (const line of lines) {
    const dist =
      Math.abs(Math.hypot(line.a[0] - x, line.a[1] - y) +
      Math.hypot(line.b[0] - x, line.b[1] - y) - Math.hypot(line.a[0] - line.b[0], line.a[1] - line.b[1]));
    if ((!minDist && dist < 2) || dist < minDist) {
      minDist = dist;
      closestLine = line;
    }
  }
  return closestLine;
}

function renderLines(sp: DG.ScatterPlotViewer, xAxis: string, yAxis: string, linesRes: IRenderedLines,
  saliVals: number[], saliOpacityCoef: number, saliMin: number): ILine [] {
  const lines = linesRes.lines;
  const canvas = (sp.getInfo() as {[index: string] : any})['canvas'];
  const ctx = canvas.getContext('2d') as CanvasRenderingContext2D;
  const x = sp.dataFrame!.columns.byName(xAxis);
  const y = sp.dataFrame!.columns.byName(yAxis);
  for (let i = 0; i < lines.length; i++) {
    const pointFrom = sp.worldToScreen(x.get(lines[i].mols[0]), y.get(lines[i].mols[0]));
    const pointTo = sp.worldToScreen(x.get(lines[i].mols[1]), y.get(lines[i].mols[1]));
    lines[i].a = [pointFrom.x, pointFrom.y];
    lines[i].b = [pointTo.x, pointTo.y];
    const line = new Path2D();
    line.moveTo(lines[i].a[0], lines[i].a[1]);
    const color = lines[i].selected ? '255,255,0' : '0,128,0';
    const opacity = saliVals[i] === Infinity ? 1 : 0.2 + (saliVals[i] - saliMin)*saliOpacityCoef;
    ctx.strokeStyle = `rgba(${color},${opacity})`;
    ctx.lineWidth = lines[i].id === linesRes.linesDf.currentRowIdx ? 3 : 1;
    line.lineTo(lines[i].b[0], lines[i].b[1]);
    ctx.stroke(line);
  }
  return lines;
}

function createLines(params: IActivityCliffsMetrics, seq: DG.Column, activities: DG.Column, semType: string,
  tags: {[index: string]: string}) : IRenderedLines {
  const lines: ILine[] = new Array(params.n1.length).fill(null);
  for (let i = 0; i < params.n1.length; i++) {
    const num1 = params.n1[i];
    const num2 = params.n2[i];
    lines[i] = ({id: i, mols: [num1, num2], selected: false, a: [], b: []});
  }
  const linesDf = DG.DataFrame.create(lines.length);
  LINES_DF_MOL_COLS_NAMES.forEach((it, idx) => {
    linesDf.columns.addNewString(it).init((i: number) => seq.get(lines[i].mols[idx]));
    setTags(linesDf.col(it)!, tags);
    linesDf.col(it)!.semType = semType;
  });
  linesDf.columns.addNewFloat(LINES_DF_ACT_DIFF_COL_NAME)
    .init((i: number) => Math.abs(activities.get(lines[i].mols[0]) - activities.get(lines[i].mols[1])));
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

export function getSimilaritiesFromDistances(dim: number, distances: Matrix, simArr: DG.Column[])
  : DG.Column[] {
  for (let i = 0; i < dim - 1; ++i) {
    const similarityArr = new Float32Array(dim - i - 1).fill(0);
    for (let j = i + 1; j < dim; ++j) {
      similarityArr[j - i - 1] = distances[i][j] === DG.FLOAT_NULL ? 0 : getSimilarityFromDistance(distances[i][j]);
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
