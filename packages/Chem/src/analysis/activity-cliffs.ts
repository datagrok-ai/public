import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {drawMoleculeToCanvas} from '../utils/chem-common-rdkit';
import {ITooltipAndPanelParams} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {convertMolNotation, getRdKitModule} from '../package';
import { RDMol } from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getMCS} from '../utils/most-common-subs';
import {BitArrayMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import { getUncommonAtomsAndBonds } from '../utils/chem-common';

const canvasWidth = 200;
const canvasHeight = 100;

const cashedData: DG.LruCache<String, any> = new DG.LruCache<string, string>();

export type ActivityCliffsParams = {
  molColName: string,
  activityColName: string,
  similarityMetric: BitArrayMetrics,
  similarity: number,
  options: any,
  isDemo?: boolean,
}

export function findMcsAndUpdateDrawings(params: ITooltipAndPanelParams, hosts: HTMLElement[]) {
  const molecule1 = params.seqCol.get(params.points[0]);
  const molecule2 = params.seqCol.get(params.points[1]);
  if (!cashedData.has(`${molecule1}_${molecule2}`))
    drawMoleculesWithMcsAsync(params, hosts, [molecule1, molecule2]);
  drawMolecules(params, hosts, [molecule1, molecule2]);
}

async function drawMoleculesWithMcsAsync(params: ITooltipAndPanelParams, hosts: HTMLElement[], molecules: [string, string]) {
  const mcsDf = DG.DataFrame.create(2);
  mcsDf.columns.addNewString('smiles').init((i) => molecules[i]);
  const mcs = await getMCS(mcsDf.col('smiles')!, true, true);
  cashedData.set(`${molecules[0]}_${molecules[1]}`, mcs);
  drawMolecules(params, hosts, molecules);
}

function drawMolecules(params: ITooltipAndPanelParams, hosts: HTMLElement[], molecules: [string, string]) {
  const rdkit = getRdKitModule();
  let mcsMol: RDMol | null = null;
  const mcsGenerated = cashedData.has(`${molecules[0]}_${molecules[1]}`);
  try {
    if (mcsGenerated)
      mcsMol = rdkit.get_qmol(cashedData.get(`${molecules[0]}_${molecules[1]}`));
    molecules.forEach((molecule: string, index: number) => {
      const imageHost = ui.canvas(canvasWidth, canvasHeight);
      if (params.seqCol.meta.units === DG.chem.Notation.Smiles) {
        //convert to molFile to draw in coordinates similar to dataframe cell
        molecule = convertMolNotation(molecule, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
      }
      const substruct = mcsGenerated ? getUncommonAtomsAndBonds(molecule!, mcsMol, rdkit) : null;
      drawMoleculeToCanvas(0, 0, canvasWidth, canvasHeight, imageHost, molecule, '',
        { normalizeDepiction: true, straightenDepiction: true }, substruct);
      ui.empty(hosts[index]);
      if (!cashedData.has(`${molecules[0]}_${molecules[1]}`))
        hosts[index].append(ui.divText('MCS loading...'));
      hosts[index].append(imageHost);
    });
  } finally {
    mcsMol?.delete();
  }
}

export function createTooltipElement(params: ITooltipAndPanelParams): HTMLDivElement {
  return createElementTemplate(params, drawTooltipElement, true, true);
}

function drawTooltipElement(params: ITooltipAndPanelParams, element: HTMLDivElement,
  hosts: HTMLDivElement[], molIdx: number) {
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
  let tooltipColsNum = Math.min(df.columns.length, 7);
  const dict: {[key: string]: string} = {};
  for (let i = 0; i < tooltipColsNum; i++) {
    const colName = df.columns.byIndex(i).name;
    colName !== seqColName ? dict[colName] = df.get(colName, idx) : tooltipColsNum++;
  }
  return ui.tableFromMap(dict);
}

function createElementTemplate(params: ITooltipAndPanelParams,
  drawFunc: (params: ITooltipAndPanelParams, element: HTMLDivElement, hosts: HTMLDivElement[],
  molIdx: number, idx: number) => void,
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
  params.points.forEach((molIdx: number, idx: number) => {
    drawFunc(params, element, hosts, molIdx, idx);
  });
  findMcsAndUpdateDrawings(params, hosts);
  return element;
}


export function createPropPanelElement(params: ITooltipAndPanelParams): HTMLDivElement {
  const propPanel = createElementTemplate(params, drawPropPanelElement, false, false,
    {style: {minWigth: `${canvasWidth}px`}});
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
