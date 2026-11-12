import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {ITooltipAndPanelParams} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {ISeqSplitted} from '@datagrok-libraries/bio/src/utils/macromolecule/types';
import {HelmType} from '@datagrok-libraries/bio/src/helm/types';

import {drawMoleculeDifferenceOnCanvas} from '../utils/cell-renderer';

import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {BitArrayMetrics} from '@datagrok-libraries/ml/src/typed-metrics';

import {_package} from '../package';

export type SeqActivityCliffsParams = {
  seqColName: string,
  activityColName: string,
  similarityMetric: MmDistanceFunctionsNames | BitArrayMetrics,
  similarity: number,
  options: any,
  isDemo?: boolean,
}

export function createTooltipElement(params: ITooltipAndPanelParams): HTMLDivElement {
  const tooltipElement = ui.divH([], {style: {gap: '10px'}});
  const columnNames = ui.divV([
    ui.divText(params.seqCol.name),
    ui.divText(params.activityCol.name),
  ]);
  columnNames.style.fontWeight = 'bold';
  columnNames.style.display = 'flex';
  columnNames.style.justifyContent = 'space-between';
  tooltipElement.append(columnNames);
  params.points.forEach((molIdx: number) => {
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

export function createPropPanelElement(params: ITooltipAndPanelParams): HTMLDivElement {
  const propPanel = ui.div();

  propPanel.append(ui.divText(params.seqCol.name, {style: {fontWeight: 'bold'}}));

  const activitiesArray = new Array<number>(2);
  params.points.forEach((molIdx, idx) => {
    activitiesArray[idx] = params.activityCol.get(molIdx);
  });

  const molDifferences: { [key: number]: HTMLCanvasElement } = {};
  const sh = _package.seqHelper.getSeqHandler(params.seqCol);
  const biotype = sh.defaultBiotype;
  const subParts1 = sh.getSplitted(params.points[0]);
  const subParts2 = sh.getSplitted(params.points[1]);
  const canvas = createDifferenceCanvas(subParts1, subParts2, biotype, molDifferences);
  propPanel.append(ui.div(canvas, {style: {width: '300px', overflow: 'scroll'}}));

  propPanel.append(createDifferencesWithPositions(molDifferences));

  propPanel.append(createPropPanelField('Activity delta', Math.abs(activitiesArray[0] - activitiesArray[1])));
  propPanel.append(createPropPanelField('Cliff', params.sali!));

  return propPanel;
}

function createPropPanelField(name: string, value: number): HTMLDivElement {
  return ui.divH([
    ui.divText(`${name}: `, {style: {fontWeight: 'bold', paddingRight: '5px'}}),
    ui.divText(value.toFixed(2)),
  ], {style: {paddingTop: '10px'}});
}

export function createDifferenceCanvas(
  subParts1: ISeqSplitted,
  subParts2: ISeqSplitted,
  biotype: HelmType,
  molDifferences: { [key: number]: HTMLCanvasElement }): HTMLCanvasElement {
  const canvas = document.createElement('canvas');
  const context = canvas.getContext('2d');
  canvas.height = 30;
  const monomerLib = _package.monomerLib;
  drawMoleculeDifferenceOnCanvas(context!, 0, 0, 0, 30,
    wu.count(0).take(subParts1.length).map((posIdx) => subParts1.getCanonical(posIdx)).toArray(),
    wu.count(0).take(subParts2.length).map((posIdx) => subParts2.getCanonical(posIdx)).toArray(),
    biotype, monomerLib, true, molDifferences);
  return canvas;
}

export function createDifferencesWithPositions(
  molDifferences: { [key: number]: HTMLCanvasElement }): HTMLDivElement {
  const div = ui.div();
  if (Object.keys(molDifferences).length > 0) {
    const diffsPanel = ui.divV([]);
    diffsPanel.append(ui.divH([
      ui.divText('Pos', {style: {fontWeight: 'bold', width: '30px', borderBottom: '1px solid'}}),
      ui.divText('Difference', {style: {fontWeight: 'bold', borderBottom: '1px solid'}}),
    ]));
    for (const key of Object.keys(molDifferences)) {
      molDifferences[key as any].style.borderBottom = '1px solid lightgray';
      diffsPanel.append(ui.divH([
        ui.divText((parseInt(key) + 1).toString(), {style: {width: '30px', borderBottom: '1px solid lightgray'}}),
        molDifferences[key as any],
      ]));
    }
    div.append(diffsPanel);
  }
  return div;
}

export function createLinesGrid(df: DG.DataFrame, colNames: string[]): DG.Grid {
  if (!df.col('seq_diff')) {
    const seqDiffCol = DG.Column.string('seq_diff', df.rowCount)
      .init((i) => `${df.get(colNames[0], i)}#${df.get(colNames[1], i)}`);
    seqDiffCol.semType = 'MacromoleculeDifference';
    seqDiffCol.meta.units = df.col(colNames[0])!.meta.units;
    seqDiffCol.setTag(bioTAGS.separator, df.col(colNames[0])!.getTag(bioTAGS.separator));
    df.columns.add(seqDiffCol);
  }
  const grid = df.plot.grid();
  grid.col(colNames[0])!.visible = false;
  grid.col(colNames[1])!.visible = false;
  return grid;
}
