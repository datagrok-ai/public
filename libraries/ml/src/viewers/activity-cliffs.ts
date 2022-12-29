import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Matrix} from '@datagrok-libraries/utils/src/type-declarations';
import {getSimilarityFromDistance} from '../distance-metrics-methods';
import {removeEmptyStringRows} from '@datagrok-libraries/utils/src/dataframe-utils';

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

let zoom = false;

const molColumnNames = ['1_seq', '2_seq'];

const nonNormalizedDistances = ['Levenshtein'];

// Searches for activity cliffs in a chemical dataset by selected cutoff
export async function getActivityCliffs( df: DG.DataFrame, seqCol: DG.Column, encodedCol: DG.Column | null,
  axesNames: string[], scatterTitle: string, activities: DG.Column, similarity: number, similarityMetric: string,
  methodName: string, semType: string, tags: {[index: string]: string},
  seqSpaceFunc: (params: ISequenceSpaceParams) => Promise<ISequenceSpaceResult>,
  simMatrixFunc: (dim: number, seqCol: DG.Column, df: DG.DataFrame, colName: string,
    simArr: DG.Column[]) => Promise<DG.Column[]>,
  tooltipFunc: (params: ITooltipAndPanelParams) => HTMLElement,
  propertyPanelFunc: (params: ITooltipAndPanelParams) => HTMLElement,
  linesGridFunc?: (df: DG.DataFrame, pairColNames: string[]) => DG.Grid,
  seqSpaceOptions?: any) : Promise<DG.Viewer> {
  const automaticSimilarityLimit = false;
  const MIN_SIMILARITY = 80;

  const initialSimilarityLimit = automaticSimilarityLimit ? MIN_SIMILARITY : similarity / 100;

  const dimensionalityReduceCol = encodedCol ?? seqCol;
  const withoutEmptyValues = DG.DataFrame.fromColumns([dimensionalityReduceCol]).clone();
  //@ts-ignore
  const emptyValsIdxs = removeEmptyStringRows(withoutEmptyValues, dimensionalityReduceCol);

  const {distance, coordinates} = await seqSpaceFunc({
    seqCol: withoutEmptyValues.col(dimensionalityReduceCol.name)!,
    methodName: methodName,
    similarityMetric: similarityMetric,
    embedAxesNames: axesNames,
    options: seqSpaceOptions,
  });

  for (const col of coordinates) {
    const listValues = col.toList();
    emptyValsIdxs.forEach((ind: number) => listValues.splice(ind, 0, null));
    df.columns.add(DG.Column.float(col.name, df.rowCount).init((i)=> listValues[i]));
  }

  const dfSeq = DG.DataFrame.fromColumns([DG.Column.fromList('string', 'seq', dimensionalityReduceCol.toList())]);
  const dim = dimensionalityReduceCol.length;
  let simArr: DG.Column[] = Array(dim - 1);

  if (!distance || emptyValsIdxs.length !== 0 || nonNormalizedDistances.includes(similarityMetric)) {
    simArr = await simMatrixFunc(dim, dimensionalityReduceCol, dfSeq, 'seq', simArr);
  } else {
    getSimilaritiesFromDistances(dim, distance, simArr);
  }

  const optSimilarityLimit = initialSimilarityLimit;

  const simVals: number[] = [];
  const saliVals: number[] = [];
  const n1: number[] = [];
  const n2: number[] = [];
  const cliffsMolIds = new Set<number>();

  for (let i = 0; i != dim - 1; ++i) {
    for (let j = 0; j != dim - 1 - i; ++j) {
      const sim: number = simArr[i] ? simArr[i].get(j) : 0;

      if (sim >= optSimilarityLimit) {
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

  const saliValsWithoutInfinity = saliVals.filter((it) => it !== Infinity);
  const saliMin = Math.min(...saliValsWithoutInfinity);
  const saliMax = Math.max(...saliValsWithoutInfinity);
  const saliOpacityCoef = 0.8/(saliMax - saliMin);


  const neighboursCount = new Array(dim).fill(0);
  const similarityCount = new Array(dim).fill(0);
  const saliCount = new Array(dim).fill(0);

  for (let i = 0; i != n1.length; ++i) {
    neighboursCount[n1[i]] += 1;
    neighboursCount[n2[i]] += 1;
    similarityCount[n1[i]] += simVals[i];
    similarityCount[n2[i]] += simVals[i];
    if (saliVals[i] != Infinity) {
      if (activities.get(n1[i]) > activities.get(n2[i])) {
        saliCount[n1[i]] += saliVals[i];
      } else {
        saliCount[n2[i]] += saliVals[i];
      }
    }
  }

  const sali: DG.Column = DG.Column.fromList('double',
    `sali_${axesNames[0].substring(axesNames[0].lastIndexOf('_'))}`, saliCount);

  df.columns.add(sali);
  const cliffsColName = addCliffsBooleanCol(df, cliffsMolIds);

  const view = grok.shell.getTableView(df.name);
  view.grid.columns.byName(cliffsColName)!.visible = false;
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
  const linesRes = createLines(n1, n2, seqCol, activities, saliVals, simVals, semType, tags);

  const cashedLinesData: any = {};
  let acc: DG.Accordion;

  linesRes.linesDf.onCurrentCellChanged.subscribe(() => {
    const currentMolIdx = linesRes.linesDf.currentCol && linesRes.linesDf.currentCol.name === '2_seq' ? 1 : 0;
    const line = linesRes.linesDf.currentRowIdx !== -1 ? linesRes.lines[linesRes.linesDf.currentRowIdx] : null;
    sp.dataFrame.currentRowIdx = line ? line.mols[currentMolIdx] : -1;
    sp.dataFrame.selection.set(0, !linesRes.lines[0].selected);
    sp.dataFrame.selection.set(0, linesRes.lines[0].selected);
    linesDfGrid.invalidate();
    if (line) {
      setTimeout(() => {
        updatePropertyPanel(df, acc, cashedLinesData, line, seqCol, activities,
          linesRes.linesDf.get('sali', line.id), propertyPanelFunc);
        const order = sp.dataFrame.getSortedOrder(view.grid.sortByColumns, view.grid.sortTypes);
        view.grid.scrollToCell(seqCol.name, order.indexOf(sp.dataFrame.currentRowIdx));
      }, 1000);
    }
  });

  linesRes.linesDf.onSelectionChanged.subscribe((_: any) => {
    if (linesRes.linesDf.mouseOverRowIdx !== -1) {
      const line = linesRes.lines[linesRes.linesDf.mouseOverRowIdx];
      line.selected = !line.selected;
      if (!line.selected) {
        df.selection.setAll(false);
      }
    }
    linesRes.lines.forEach((l) => {
      if (l.selected) {
        l.mols.forEach((m) => df.selection.set(m, true));
      }
    });
    linesDfGrid.invalidate();
  });

  const linesDfGrid = linesGridFunc ?
    linesGridFunc(linesRes.linesDf, molColumnNames).sort(['sali'], [false]) :
    linesRes.linesDf.plot.grid().sort(['sali'], [false]);

  linesDfGrid.onCellClick.subscribe(() => {
    zoom = true;
  });

  const listCliffsLink = ui.button(`${linesRes.linesDf.rowCount} cliffs`, () => {
    view.dockManager.dock(linesDfGrid.root, 'down', null, 'Activity cliffs', 0.2);
  });

  const viewerClosedSub = grok.events.onViewerClosed.subscribe((v) => {
    //@ts-ignore
    if(v.args.viewer === sp) {
      view.dockManager.close(linesDfGrid.root);
      viewerClosedSub.unsubscribe();
      view.subs = view.subs.filter((sub) => sub !== viewerClosedSub);
    }  
  })
  view.subs.push(viewerClosedSub);

  editSpLinkStyle(listCliffsLink, sp, 10);

  const filterCliffsButton = ui.switchInput(`Show only cliffs`, false);
  filterCliffsButton.onChanged(()=> {
    filterCliffsButton.value ? df.rows.filter((row) => row[cliffsColName] === true) : df.filter.setAll(true);
  });
  editSpLinkStyle(filterCliffsButton.root, sp, 30);

  let timer: NodeJS.Timeout;
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
    const line = checkCursorOnLine(event, canvas, linesRes.lines);
    if (line && df.mouseOverRowIdx === -1) {
      if (event.ctrlKey) {
        line.selected = !line.selected;
        linesRes.linesDf.selection.set(line.id, line.selected);

        if (!line.selected) {
          df.selection.setAll(false);
        }

        linesRes.lines.forEach((l) => {
          if (l.selected) {
            l.mols.forEach((m) => df.selection.set(m, true));
          }
        });
      } else {
        if (linesRes.linesDf.currentRowIdx !== line.id) {
          linesRes.linesDf.currentRowIdx = line.id;
          df.currentRowIdx = line.mols[0];
          df.selection.set(0, !linesRes.lines[0].selected);
          df.selection.set(0, linesRes.lines[0].selected);
        }
      }
      const order = linesRes.linesDf.getSortedOrder(linesDfGrid.sortByColumns, linesDfGrid.sortTypes);
      linesDfGrid.scrollToCell(molColumnNames[0], order.indexOf(line.id));
    }
  });

  sp.onEvent('d4-before-draw-scene')
    .subscribe((_: any) => {
      const lines = renderLines(sp,
        axesNames[0], axesNames[1], linesRes, saliVals, saliOpacityCoef, saliMin);
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
      filterCliffsButton.value ? df.rows.filter((row) => row[cliffsColName] === true) : df.filter.setAll(true);
    });

  sp.addProperty('similarityLimit', 'double', optSimilarityLimit);
  acc = createPopertyPanel();
  return sp;
}

function editSpLinkStyle(el: HTMLElement, sp: DG.ScatterPlotViewer, topMargin: number) {
  el.style.position = 'absolute';
  el.style.top = `${topMargin}px`;
  el.style.right = '10px';
  sp.root.append(el);
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

function createLines(n1: number[], n2: number[], seq: DG.Column, activities: DG.Column, saliVals: number[],
  simVals: number[], semType: string, tags: {[index: string]: string}) : IRenderedLines {
  const lines: ILine[] = new Array(n1.length).fill(null);
  for (let i = 0; i < n1.length; i++) {
    const num1 = n1[i];
    const num2 = n2[i];
    lines[i] = ({id: i, mols: [num1, num2], selected: false, a: [], b: []});
  }
  const linesDf = DG.DataFrame.create(lines.length);
  molColumnNames.forEach((it, idx) => {
    linesDf.columns.addNewString(it).init((i: number) => seq.get(lines[i].mols[idx]));
    setTags(linesDf.col(it)!, tags);
    linesDf.col(it)!.semType = semType;
  });
  linesDf.columns.addNewFloat('act_diff')
    .init((i: number) => Math.abs(activities.get(lines[i].mols[0]) - activities.get(lines[i].mols[1])));
  linesDf.columns.addNewInt('line_index').init((i: number) => i);
  linesDf.columns.addNewFloat('sali').init((i: number) => saliVals[i]);
  linesDf.columns.addNewFloat('sim').init((i: number) => simVals[i]);
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
      similarityArr[j - i - 1] = getSimilarityFromDistance(distances[i][j]);
    }
    simArr[i] = DG.Column.fromFloat32Array('similarity', similarityArr);
  }
  return simArr;
}

export function addCliffsBooleanCol(df: DG.DataFrame, ids: Set<number>): string {
  const colname = 'containsCliff';
  const colNameInd = df.columns.names().filter((it: string) => it.includes(colname)).length + 1;
  const newColName = `${colname}_${colNameInd}`;
  df.columns.addNewBool(newColName).init((i) => ids.has(i));
  return newColName;
}
