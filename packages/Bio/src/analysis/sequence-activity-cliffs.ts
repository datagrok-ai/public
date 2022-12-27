import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ITooltipAndPanelParams} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {getSimilarityFromDistance} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import {AvailableMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {TAGS} from '../utils/constants';
import {drawMoleculeDifferenceOnCanvas} from '../utils/cell-renderer';
import * as C from '../utils/constants';
import {GridColumn} from 'datagrok-api/dg';
import {invalidateMols, MONOMERIC_COL_TAGS} from '../substructure-search/substructure-search';
import {getSplitter} from '@datagrok-libraries/bio';

export async function getDistances(col: DG.Column, seq: string): Promise<Array<number>> {
  const stringArray = col.toList();
  const distances = new Array(stringArray.length).fill(0);
  for (let i = 0; i < stringArray.length; ++i) {
    const distance = stringArray[i] ? AvailableMetrics['String']['Levenshtein'](stringArray[i], seq) : null;
    distances[i] = distance ? distance / Math.max((stringArray[i] as string).length, seq.length) : null;
  }
  return distances;
}

export async function getSimilaritiesMarix(dim: number, seqCol: DG.Column, df: DG.DataFrame, colName: string, simArr: DG.Column[])
  : Promise<DG.Column[]> {
  const distances = new Array(simArr.length).fill(null);
  for (let i = 0; i != dim - 1; ++i) {
    const seq: string = seqCol.get(i);
    df.rows.removeAt(0, 1, false);
    distances[i] = (await getDistances(df.col(colName)!, seq))!;
  }

  for (let i = 0; i < distances.length; i++) {
    for (let j = 0; j < distances[i].length; j++)
      distances[i][j] = getSimilarityFromDistance(distances[i][j]);

    simArr[i] = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'distances', distances[i]);
  }
  return simArr;
}

export async function getChemSimilaritiesMarix(dim: number, seqCol: DG.Column,
  df: DG.DataFrame, colName: string, simArr: DG.Column[])
  : Promise<DG.Column[]> {
  if (seqCol.version !== seqCol.temp[MONOMERIC_COL_TAGS.LAST_INVALIDATED_VERSION])
    await invalidateMols(seqCol, false);
  const fpDf = DG.DataFrame.create(seqCol.length);
  fpDf.columns.addNewString(colName).init((i) => seqCol.temp[MONOMERIC_COL_TAGS.MONOMERIC_MOLS].get(i));
  const res = await grok.functions.call('Chem:getChemSimilaritiesMatrix', {
    dim: dim,
    col: seqCol.temp[MONOMERIC_COL_TAGS.MONOMERIC_MOLS],
    df: fpDf,
    colName: colName,
    simArr: simArr
  });
  return res;
}

export function createTooltipElement(params: ITooltipAndPanelParams): HTMLDivElement {
  const tooltipElement = ui.divH([]);
  const columnNames = ui.divV([
    ui.divText(params.seqCol.name),
    ui.divText(params.activityCol.name),
  ]);
  columnNames.style.fontWeight = 'bold';
  columnNames.style.display = 'flex';
  columnNames.style.justifyContent = 'space-between';
  tooltipElement.append(columnNames);
  params.line.mols.forEach((molIdx: number, idx: number) => {
    const activity = ui.divText(params.activityCol.get(molIdx).toFixed(2));
    activity.style.display = 'flex';
    activity.style.justifyContent = 'left';
    activity.style.paddingLeft = '30px';
    tooltipElement.append(ui.divV([
      ui.divText(params.seqCol.get(molIdx)),
      activity,
    ]));
  });
  return tooltipElement;
}

function moleculeInfo(df: DG.DataFrame, idx: number, seqColName: string): HTMLElement {
  const dict: { [key: string]: string } = {};
  for (const col of df.columns) {
    if (col.name !== seqColName)
      dict[col.name] = df.get(col.name, idx);
  }
  return ui.tableFromMap(dict);
}


export function createPropPanelElement(params: ITooltipAndPanelParams): HTMLDivElement {
  const propPanel = ui.div();

  propPanel.append(ui.divText(params.seqCol.name, {style: {fontWeight: 'bold'}}));

  const sequencesArray = new Array<string>(2);
  const activitiesArray = new Array<number>(2);
  params.line.mols.forEach((molIdx, idx) => {
    sequencesArray[idx] = params.seqCol.get(molIdx);
    activitiesArray[idx] = params.activityCol.get(molIdx);
  });

  const molDifferences: { [key: number]: HTMLCanvasElement } = {};
  const units = params.seqCol.getTag(DG.TAGS.UNITS);
  const separator = params.seqCol.getTag(TAGS.SEPARATOR);
  const splitter = getSplitter(units, separator);
  const subParts1 = splitter(sequencesArray[0]);
  const subParts2 = splitter(sequencesArray[1]);
  const canvas = createDifferenceCanvas(subParts1, subParts2, units, molDifferences);
  propPanel.append(ui.div(canvas, {style: {width: '300px', overflow: 'scroll'}}));

  propPanel.append(createDifferencesWithPositions(molDifferences));

  propPanel.append(createPropPanelField('Activity delta', Math.abs(activitiesArray[0] - activitiesArray[1])));
  propPanel.append(createPropPanelField('Cliff', params.sali!));

  return propPanel;
}

function createPropPanelField(name: string, value: number): HTMLDivElement {
  return ui.divH([
    ui.divText(`${name}: `, {style: {fontWeight: 'bold', paddingRight: '5px'}}),
    ui.divText(value.toFixed(2))
  ], {style: {paddingTop: '10px'}});
}

export function createDifferenceCanvas(
  subParts1: string[],
  subParts2: string[],
  units: string,
  molDifferences: { [key: number]: HTMLCanvasElement }): HTMLCanvasElement {
  const canvas = document.createElement('canvas');
  const context = canvas.getContext('2d');
  canvas.height = 30;
  drawMoleculeDifferenceOnCanvas(context!, 0, 0, 0, 30, subParts1, subParts2, units, true, molDifferences);
  return canvas;
}

export function createDifferencesWithPositions(
  molDifferences: { [key: number]: HTMLCanvasElement }): HTMLDivElement {
  const div = ui.div();
  if (Object.keys(molDifferences).length > 0) {
    const diffsPanel = ui.divV([]);
    diffsPanel.append(ui.divH([
      ui.divText('Pos', {style: {fontWeight: 'bold', width: '30px', borderBottom: '1px solid'}}),
      ui.divText('Difference', {style: {fontWeight: 'bold', borderBottom: '1px solid'}})
    ]));
    for (const key of Object.keys(molDifferences)) {
      molDifferences[key as any].style.borderBottom = '1px solid lightgray';
      diffsPanel.append(ui.divH([
        ui.divText((parseInt(key) + 1).toString(), {style: {width: '30px', borderBottom: '1px solid lightgray'}}),
        molDifferences[key as any]
      ]));
    }
    div.append(diffsPanel);
  }
  return div;
}

export function createLinesGrid(df: DG.DataFrame, colNames: string[]): DG.Grid {
  const seqDiffCol = DG.Column.string('seq_diff', df.rowCount)
    .init((i) => `${df.get(colNames[0], i)}#${df.get(colNames[1], i)}`);
  seqDiffCol.semType = 'MacromoleculeDifference';
  seqDiffCol.setTag(DG.TAGS.UNITS, df.col(colNames[0])!.getTag(DG.TAGS.UNITS));
  seqDiffCol.setTag(C.TAGS.SEPARATOR, df.col(colNames[0])!.getTag(C.TAGS.SEPARATOR));
  df.columns.add(seqDiffCol);
  const grid = df.plot.grid();
  grid.col(colNames[0])!.visible = false;
  grid.col(colNames[1])!.visible = false;
  return grid;
}
