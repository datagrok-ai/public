import {ITooltipAndPanelParams} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {getSimilarityFromDistance} from '@datagrok-libraries/utils/src/similarity-metrics';
import {AvailableMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import * as grok from 'datagrok-api/grok';
import { SplitterFunc, WebLogo } from '@datagrok-libraries/bio/src/viewers/web-logo';
import { UnitsHandler } from '@datagrok-libraries/bio/src/utils/units-handler';

export async function getDistances(col: DG.Column, seq: string): Promise<Array<number>> {
  const stringArray = col.toList();
  const distances = new Array(stringArray.length).fill(0);
  for (let i = 0; i < stringArray.length; ++i) {
    const distance = stringArray[i] ? AvailableMetrics['String']['Levenshtein'](stringArray[i], seq) : null;
    distances[i] = distance ? distance/Math.max((stringArray[i] as string).length, seq.length) : null;
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
    for (let j = 0; j < distances[i].length; j++) {
      distances[i][j] = getSimilarityFromDistance(distances[i][j]);
    }
    simArr[i] = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'distances', distances[i]);
  }
  return simArr;
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
  let dict: {[key: string]: string} = {};
  for (let col of df.columns) {
    if(col.name !== seqColName) {
      dict[col.name] = df.get(col.name, idx);
    }
  }
  return ui.tableFromMap(dict);
}


export function createPropPanelElement(params: ITooltipAndPanelParams): HTMLDivElement {
  const propPanel = ui.divV([]);
  const columnNames = ui.divH([
    ui.divText(params.seqCol.name),
    ui.divText(params.activityCol.name),
  ]);
  columnNames.style.fontWeight = 'bold';
  columnNames.style.justifyContent = 'space-between';
  propPanel.append(columnNames);
  const hosts: HTMLDivElement[] = [];
  params.line.mols.forEach((molIdx: number, hostIdx: number) => {
    const activity = ui.divText(params.activityCol.get(molIdx).toFixed(2));
    activity.style.paddingLeft = '15px';
    activity.style.paddingLeft = '10px';
    const molHost = ui.divText(params.seqCol.get(molIdx));
    if (params.df.currentRowIdx === molIdx) {
      molHost.style.border = 'solid 1px lightgrey';
    }
    //@ts-ignore
    ui.tooltip.bind(molHost, () => moleculeInfo(params.df, molIdx, params.seqCol.name));
    molHost.onclick = () => {
      const obj = grok.shell.o;
      molHost.style.border = 'solid 1px lightgrey';
      params.df.currentRowIdx = molIdx;
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
  propPanel.append(ui.divH([
    ui.divText(`Cliff: `, {style: {fontWeight: 'bold', paddingRight: '5px'}}),
    ui.divText(params.sali!.toFixed(2))
  ]));
  return propPanel;
}