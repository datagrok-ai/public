import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import { Matrix } from '@datagrok-libraries/utils/src/type-declarations';
import { getSimilarityFromDistance } from '@datagrok-libraries/utils/src/similarity-metrics';
import {removeEmptyStringRows} from '@datagrok-libraries/utils/src/dataframe-utils'

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
  hosts: HTMLDivElement[],
  line: ILine,
  seqCol: DG.Column,
}

let zoom = false;

const molColumnNames = ['1_seq', '2_seq'];

// Searches for activity cliffs in a chemical dataset by selected cutoff
export async function getActivityCliffs(
    df: DG.DataFrame, 
    seqCol: DG.Column,
    encodedCol: DG.Column, 
    axesNames: string[],
    scatterTitle: string,
    activities: DG.Column, 
    similarity: number,
    similarityMetric: string, 
    methodName: string, 
    semType: string,
    tags: {[index: string]: string},
    seqSpaceFunc: (params: ISequenceSpaceParams) => Promise<ISequenceSpaceResult>,
    simFunc: (col: DG.Column, mol: string) => Promise<DG.Column | null>,
    renderSeqFunction: (params: ITooltipAndPanelParams) => void,
    options?: any) : Promise<DG.Viewer> {
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
    options: options
  });

  for (const col of coordinates) {
      const listValues = col.toList();
      emptyValsIdxs.forEach((ind: number) => listValues.splice(ind, 0, null));
      df.columns.add(DG.Column.fromList('double', col.name, listValues));
  }

  const dfSeq = DG.DataFrame.fromColumns([DG.Column.fromList('string', 'seq', dimensionalityReduceCol.toList())]);
  const dim = dimensionalityReduceCol.length;
  const simArr: DG.Column[] = Array(dim - 1);

  if (!distance || emptyValsIdxs.length !== 0)
    await getSimilaritiesMarix(dim, dimensionalityReduceCol, dfSeq, simArr, simFunc);
  else
    getSimilaritiesMarixFromDistances(dim, distance, simArr);

  const optSimilarityLimit = initialSimilarityLimit;

  const simVals: number[] = [];
  const saliVals: number[] = [];
  const n1: number[] = [];
  const n2: number[] = [];

  for (let i = 0; i != dim - 1; ++i) {
    for (let j = 0; j != dim - 1 - i; ++j) {
      const sim: number = simArr[i] ? simArr[i].get(j) : 0;

      if (sim >= optSimilarityLimit) {
        n1.push(i);
        n2.push(i + j + 1);
        simVals.push(sim);
        const diff = Math.abs(activities.get(i) - activities.get(i + j + 1));
        if (sim != 1)
          saliVals.push(diff / (1 - sim));
        else
          saliVals.push(Infinity);
      }
    }
  }

  const saliValsWithoutInfinity = saliVals.filter(it => it !== Infinity);
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
      if (activities.get(n1[i]) > activities.get(n2[i]))
        saliCount[n1[i]] += saliVals[i];
      else
        saliCount[n2[i]] += saliVals[i];
    }
  }

  const sali: DG.Column = DG.Column.fromList('double', `sali_${axesNames[0].substring(axesNames[0].lastIndexOf('_'))}`, saliCount);

  df.columns.add(sali);

  const view = grok.shell.getTableView(df.name);
  const sp = view.addViewer(DG.Viewer.scatterPlot(df, {
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
    title: scatterTitle
  })) as DG.ScatterPlotViewer;

  const canvas = (sp.getInfo() as any)['canvas'];
  const linesRes = createLines(n1, n2, seqCol, activities, saliVals, semType, tags);
  const cashedLinesData: any = {};
  let acc: DG.Accordion;

  linesRes.linesDf.onCurrentCellChanged.subscribe(() => {
    const currentMolIdx = linesRes.linesDf.currentCol && linesRes.linesDf.currentCol.name === '2_seq' ? 1 : 0;
    sp.dataFrame.currentRowIdx =
      linesRes.linesDf.currentRowIdx !== -1 ? linesRes.lines[linesRes.linesDf.currentRowIdx].mols[currentMolIdx] : -1;
    sp.dataFrame.selection.set(0, !linesRes.lines[0].selected);
    sp.dataFrame.selection.set(0, linesRes.lines[0].selected);
    linesDfGrid.invalidate();
  });

  linesRes.linesDf.onSelectionChanged.subscribe((_) => {
    if (linesRes.linesDf.mouseOverRowIdx !== -1) {
      const line = linesRes.lines[linesRes.linesDf.mouseOverRowIdx];
      line.selected = !line.selected;
      if (!line.selected)
        df.selection.setAll(false);
    }
    linesRes.lines.forEach((l) => {
      if (l.selected)
        l.mols.forEach((m) => df.selection.set(m, true));
    });
    linesDfGrid.invalidate();
  });

  const linesDfGrid = linesRes.linesDf.plot.grid().sort(['sali'], [false]);
  linesDfGrid.root.style.minWidth = '650px';

  linesDfGrid.onCellClick.subscribe(() => {
    zoom = true;
  });

  const listCliffsLink = ui.button(`${linesRes.linesDf.rowCount} cliffs`, () => {
    const cliffsDialog = ui.dialog({title: 'Activity cliffs'})
      .add(linesDfGrid.root)
      .show({resizable: true});
    ui.tools.waitForElementInDom(linesDfGrid.root).then(() => {
      molColumnNames.forEach(it => { linesDfGrid.col(it)!.width = 200; });
      linesDfGrid.invalidate();
    });
    cliffsDialog.root.id = 'cliffs_dialog';
  });
  listCliffsLink.style.position = 'absolute';
  listCliffsLink.style.top = '10px';
  listCliffsLink.style.right = '10px';
  sp.root.append(listCliffsLink);

  let timer: NodeJS.Timeout;
  canvas.addEventListener('mousemove', function (event: MouseEvent) {
    clearTimeout(timer);
    timer = global.setTimeout(function () {
      const line = checkCursorOnLine(event, canvas, linesRes.lines);
      if (line && df.mouseOverRowIdx === -1) {
          ui.tooltip.show(createTooltipElement(cashedLinesData, line, seqCol, activities, renderSeqFunction), event.clientX, event.clientY);
      }
    }, 500);
  });

  canvas.addEventListener('mousedown', function(event: MouseEvent) {
    const line = checkCursorOnLine(event, canvas, linesRes.lines);
    if (line && df.mouseOverRowIdx === -1) {
      if (event.ctrlKey) {
        line.selected = !line.selected;
        linesRes.linesDf.selection.set(line.id, line.selected);

        if (!line.selected)
          df.selection.setAll(false);

        linesRes.lines.forEach((l) => {
          if (l.selected)
            l.mols.forEach((m) => df.selection.set(m, true));
        });
      } else {
        if (linesRes.linesDf.currentRowIdx !== line.id) {
          linesRes.linesDf.currentRowIdx = line.id;
          df.currentRowIdx = line.mols[0];
          df.selection.set(0, !linesRes.lines[0].selected);
          df.selection.set(0, linesRes.lines[0].selected);
        }
      }
      updatePropertyPanel(df, acc, cashedLinesData, line, seqCol, activities, linesRes.linesDf.get('sali', line.id), renderSeqFunction);
    }
  });

  sp.onEvent('d4-before-draw-scene')
    .subscribe((_) => {
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
            sp.dataFrame.get(axesNames[1], currentLine.mols[1])
          )        
          sp.zoom(zoomLeft,
            zoomTop,
            zoomRight,
            zoomBottom);
        }, 300);
        zoom = false;
      }
    });

  sp.addProperty('similarityLimit', 'double', optSimilarityLimit);
  acc = createPopertyPanel();
  return sp;
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

function updatePropertyPanel(
  df: DG.DataFrame,
  acc: DG.Accordion,
  cashedData: any,
  line: ILine, 
  seqCol: DG.Column, 
  activities: DG.Column,
  sali: number,
  renderSeqFunction: (params: ITooltipAndPanelParams) => void){
  const panel = acc.getPane('Cliff Details');
  ui.empty(panel.root);
  const panelElement = createPropPanelElement(df, cashedData, line, seqCol, activities, sali, renderSeqFunction);
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

function renderLines(sp: DG.ScatterPlotViewer,
  xAxis: string, yAxis: string, linesRes: IRenderedLines, saliVals: number[], saliOpacityCoef: number, saliMin: number): ILine [] {
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

function createLines(
  n1: number[], 
  n2: number[], 
  seq: DG.Column, 
  activities: DG.Column, 
  saliVals: number[],
  semType: string,
  tags: {[index: string]: string}) : IRenderedLines {
  const lines: ILine[] = [];
  for (let i = 0; i < n1.length; i++) {
    const num1 = n1[i];
    const num2 = n2[i];
    lines.push({id: i, mols: [num1, num2], selected: false, a: [], b: []});
  }
  const linesDf = DG.DataFrame.create(lines.length);
  molColumnNames.forEach((it, idx) => {
    linesDf.columns.addNewString(it).init((i: number) => seq.get(lines[i].mols[idx]));
    setTags(linesDf.col(it)!, tags);
    linesDf.col(it)!.semType = semType;
  })
  linesDf.columns.addNewFloat('act_diff')
    .init((i: number) => Math.abs(activities.get(lines[i].mols[0]) - activities.get(lines[i].mols[1])));
  linesDf.columns.addNewInt('line_index').init((i: number) => i);
  linesDf.columns.addNewFloat('sali').init((i: number) => saliVals[i]);
  return {lines, linesDf};
}

function setTags(col: DG.Column, tags: {[index: string]: string}) {
  Object.keys(tags).forEach(tag => {
    col.tags[tag] = tags[tag];
  })
}

export async function getSimilaritiesMarix(
    dim: number, 
    seqCol: DG.Column, 
    dfSeq: DG.DataFrame, 
    simArr: DG.Column[],
    simFunc: (col: DG.Column, mol: string) => Promise<DG.Column | null>
    )
  : Promise<DG.Column[]> {
  for (let i = 0; i != dim - 1; ++i) {
    const mol = seqCol.get(i);
    dfSeq.rows.removeAt(0, 1, false);
    simArr[i] = (await simFunc(dfSeq.col('seq')!, mol))!;
  }
  return simArr;
}

export function getSimilaritiesMarixFromDistances(dim: number, distances: Matrix, simArr: DG.Column[])
  : DG.Column[] {
  for (let i = 0; i < dim - 1; ++i) {
    const similarityArr = [];
    for (let j = i + 1; j < dim; ++j)
      similarityArr.push(getSimilarityFromDistance(distances[i][j]));
    simArr[i] = DG.Column.fromFloat32Array('similarity', Float32Array.from(similarityArr));
  }
  return simArr;
}


export function createPropPanelElement(
  df: DG.DataFrame,
  cashedData: any,
  line: ILine,
  seqCol: DG.Column,
  activities: DG.Column,
  sali: number,
  renderSeqFunction: (params: ITooltipAndPanelParams) => void): HTMLDivElement {
  const propPanel = ui.divV([]);
  const columnNames = ui.divH([
    ui.divText(seqCol.name),
    ui.divText(activities.name),
  ]);
  columnNames.style.fontWeight = 'bold';
  columnNames.style.justifyContent = 'space-between';
  propPanel.append(columnNames);
  const hosts: HTMLDivElement[] = [];
  line.mols.forEach((molIdx: number, hostIdx: number) => {
    const activity = ui.divText(activities.get(molIdx).toFixed(2));
    activity.style.paddingLeft = '15px';
    activity.style.paddingLeft = '10px';
    const molHost = ui.div();
    if (df.currentRowIdx === molIdx) {
      molHost.style.border = 'solid 1px lightgrey';
    }
    ui.tooltip.bind(molHost, () => moleculeInfo(df, molIdx, seqCol.name));
    molHost.onclick = () => {
      const obj = grok.shell.o;
      molHost.style.border = 'solid 1px lightgrey';
      df.currentRowIdx = molIdx;
      hosts.forEach((h, i) => {
        if (i !== hostIdx) {
          h.style.border = '';
        }
      })
      setTimeout(() => {
        grok.shell.o = obj
      }, 1000);
    };
    propPanel.append(ui.divH([
      molHost,
      activity,
    ]));
    hosts.push(molHost);
  });
  renderSeqFunction({ cashedData: cashedData, hosts: hosts, line: line, seqCol: seqCol});
  propPanel.append(ui.divH([
    ui.divText(`SALI: `, {style: {fontWeight: 'bold', paddingRight: '5px'}}),
    ui.divText(sali.toFixed(2))
  ]));
  return propPanel;
}

function moleculeInfo(df: DG.DataFrame, idx: number, seqColName: string): HTMLElement {
  let dict: {[key: string]: string} = {};
  for (let col of df.columns) {
    if(col.name !== seqColName) {
      dict[col.name] = df.get(col.name, idx);
    }
  }
  return ui.tableFromMap(dict);
}

export function createTooltipElement(
  cashedData: any,
  line: ILine,
  seqCol: DG.Column,
  activities: DG.Column,
  renderSeqFunction: (params: ITooltipAndPanelParams) => void): HTMLDivElement {
  const tooltipElement = ui.divH([]);
  const columnNames = ui.divV([
    ui.divText(seqCol.name),
    ui.divText(activities.name),
  ]);
  columnNames.style.fontWeight = 'bold';
  columnNames.style.display = 'flex';
  columnNames.style.justifyContent = 'space-between';
  tooltipElement.append(columnNames);
  const hosts: HTMLDivElement[] = [];
  line.mols.forEach((molIdx: number, idx: number) => {
    const activity = ui.divText(activities.get(molIdx).toFixed(2));
    activity.style.display = 'flex';
    activity.style.justifyContent = 'left';
    activity.style.paddingLeft = '30px';
    const molHost = ui.div();
    tooltipElement.append(ui.divV([
      molHost,
      activity,
    ]));
    hosts.push(molHost);
  });
  renderSeqFunction({ cashedData: cashedData, hosts: hosts, line: line, seqCol: seqCol });
  return tooltipElement;
}