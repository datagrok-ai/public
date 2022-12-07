import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {findMCS} from '../scripts-api';
import {drawMoleculeToCanvas} from '../utils/chem-common-rdkit';
import {ITooltipAndPanelParams} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';

const canvasWidth = 200;
const canvasHeight = 100;

export function findMcsAndUpdateDrawings(params: ITooltipAndPanelParams, hosts: HTMLElement[]) {
  if (!params.cashedData[params.line.id])
    drawMoleculesWithMcsAsync(params, hosts);

  drawMolecules(params, hosts);
}

async function drawMoleculesWithMcsAsync(params: ITooltipAndPanelParams, hosts: HTMLElement[]) {
  const mcsDf = DG.DataFrame.create(2);
  mcsDf.columns.addNewString('smiles').init((i) => params.seqCol.get(params.line.mols[i]));
  const mcs = await findMCS('smiles', mcsDf, true);
  params.cashedData[params.line.id] = mcs;
  drawMolecules(params, hosts);
}

function drawMolecules(params: ITooltipAndPanelParams, hosts: HTMLElement[]) {
  params.line.mols.forEach((mol: number, index: number) => {
    const imageHost = ui.canvas(canvasWidth, canvasHeight);
    const r = window.devicePixelRatio;
    imageHost.width = canvasWidth * r;
    imageHost.height = canvasHeight * r;
    imageHost.style.width = (canvasWidth).toString() + 'px';
    imageHost.style.height = (canvasHeight).toString() + 'px';
    drawMoleculeToCanvas(0, 0, canvasWidth, canvasHeight, imageHost, params.seqCol.get(mol), params.cashedData[params.line.id],
      {normalizeDepiction: false, straightenDepiction: false});
    ui.empty(hosts[index]);
    if (!params.cashedData[params.line.id])
      hosts[index].append(ui.divText('MCS loading...'));
    hosts[index].append(imageHost);
  });
}

export function createTooltipElement(params: ITooltipAndPanelParams): HTMLDivElement {
  return createElementTemplate(params, drawTooltipElement, true, true);
}

function drawTooltipElement(params: ITooltipAndPanelParams, element: HTMLDivElement,
  hosts: HTMLDivElement[], molIdx: number, idx: number) {
  const activity = ui.divText(params.activityCol.get(molIdx).toFixed(2));
  activity.style.display = 'flex';
  activity.style.justifyContent = 'left';
  activity.style.paddingLeft = '30px';
  const molHost = ui.div();
  element.append(ui.divV([
    molHost,
    activity,
  ]));
  hosts.push(molHost);
}

function moleculeInfo(df: DG.DataFrame, idx: number, seqColName: string): HTMLElement {
  const dict: {[key: string]: string} = {};
  for (const col of df.columns) {
    if (col.name !== seqColName) 
      dict[col.name] = df.get(col.name, idx);
  }
  return ui.tableFromMap(dict);
}

function createElementTemplate(params: ITooltipAndPanelParams,
  drawFunc: (params: ITooltipAndPanelParams, element: HTMLDivElement, hosts: HTMLDivElement[], molIdx: number, idx: number) => void,
  vertical: boolean, flexColNames: boolean, colNameStyle?: any) {
  const element = vertical ? ui.divH([]) : ui.divV([]);
  const columnNames = vertical ? ui.divV([]) : ui.divH([]);
  columnNames.append(colNameStyle ? ui.divText(params.seqCol.name, colNameStyle) : ui.divText(params.seqCol.name));
  columnNames.append(ui.divText(params.activityCol.name));
  columnNames.style.fontWeight = 'bold';
  columnNames.style.justifyContent = 'space-between';
  if (flexColNames) columnNames.style.display = 'flex';
  element.append(columnNames);
  const hosts: HTMLDivElement[] = [];
  params.line.mols.forEach((molIdx: number, idx: number) => {
    drawFunc(params, element, hosts, molIdx, idx);
  });
  findMcsAndUpdateDrawings(params, hosts);
  return element;
}


export function createPropPanelElement(params: ITooltipAndPanelParams): HTMLDivElement {
  const propPanel = createElementTemplate(params, drawPropPanelElement, false, false, {style: {minWigth: `${canvasWidth}px`}});
  propPanel.append(ui.divH([
    ui.divText(`Cliff: `, {style: {fontWeight: 'bold', paddingRight: '5px'}}),
    ui.divText(params.sali!.toFixed(2)),
  ]));
  return propPanel;
}

function drawPropPanelElement(params: ITooltipAndPanelParams, element: HTMLDivElement,
  hosts: HTMLDivElement[], molIdx: number, idx: number) {
  const activity = ui.divText(params.activityCol.get(molIdx).toFixed(2));
  activity.style.paddingLeft = '15px';
  activity.style.paddingLeft = '10px';
  const molHost = ui.div();
  if (params.df.currentRowIdx === molIdx) 
    molHost.style.border = 'solid 1px lightgrey';
    
  ui.tooltip.bind(molHost, () => moleculeInfo(params.df, molIdx, params.seqCol.name));
  molHost.onclick = () => {
    const obj = grok.shell.o;
    molHost.style.border = 'solid 1px lightgrey';
    params.df.currentRowIdx = molIdx;
    hosts.forEach((h, i) => {
      if (i !== idx) 
        h.style.border = '';
    });
    setTimeout(() => {
      grok.shell.o = obj;
    }, 1000);
  };
  element.append(ui.divH([
    molHost,
    activity,
  ]));
  hosts.push(molHost);
}
