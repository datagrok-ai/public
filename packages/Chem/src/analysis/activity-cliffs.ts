import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {drawMoleculeToCanvas} from '../utils/chem-common-rdkit';
import {ITooltipAndPanelParams} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {PackageFunctions} from '../package';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {ISubstruct} from '@datagrok-libraries/chem-meta/src/types';
import {getMCS} from '../utils/most-common-subs';
import {BitArrayMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {getUncommonAtomsAndBonds} from '../utils/chem-common';
import {byEngineAsync, isCruxRequested} from '../engine/chem-engine';
import {getCruxService} from '../crux/crux-searches';

const canvasWidth = 200;
const canvasHeight = 100;

/** Each pair's {@link cliffMcs}, by {@link pairKey}. */
const cachedData = new DG.LruCache<string, string | ISubstruct[] | null>();
/** A pair's key in {@link cachedData}: per engine (Crux mode's MCES, RDKit's FMCS), so that a switch computes it
 * again. */
const pairKey = (molecule1: string, molecule2: string) =>
  `${isCruxRequested() ? 'Crux' : 'RDKit'}:${molecule1}_${molecule2}`;

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
  if (!cachedData.has(pairKey(molecule1, molecule2)))
    drawMoleculesWithMcsAsync(params, hosts, [molecule1, molecule2]);
  drawMolecules(params, hosts, [molecule1, molecule2]);
}

async function drawMoleculesWithMcsAsync(params: ITooltipAndPanelParams, hosts: HTMLElement[],
  molecules: [string, string]) {
  const key = pairKey(molecules[0], molecules[1]);
  const mcs = await cliffMcs(molecules);
  cachedData.set(key, mcs);
  drawMolecules(params, hosts, molecules);
}

/**
 * The pair's MCS, which the tooltip and the panel highlight the rest of: the SMARTS of RDKit's FMCS from the RDKit
 * workers; in Crux mode ("Chemical engine"), each molecule's atoms and bonds outside RDKit's Rascal MCES from the Crux
 * workers (the most bonds the two share, possibly in several fragments; null when they share too little), with FMCS
 * behind it. Rascal's own atom matches, not its SMARTS matched back: the fragments of a disconnected MCES can match
 * elsewhere.
 */
export async function cliffMcs(molecules: [string, string]): Promise<string | ISubstruct[] | null> {
  const fmcs = () => {
    const mcsDf = DG.DataFrame.create(2);
    mcsDf.columns.addNewString('smiles').init((i) => molecules[i]);
    return getMCS(mcsDf.col('smiles')!, true, true);
  };
  return byEngineAsync<string | ISubstruct[] | null>('activity cliffs MCS',
    async () => (await getCruxService()).pairUncommonAtomsAndBonds(molecules[0], molecules[1]), fmcs);
}

function drawMolecules(params: ITooltipAndPanelParams, hosts: HTMLElement[], molecules: [string, string]) {
  // a SMARTS is FMCS's (cliffMcs' RDKit side), highlighted on RDKit's module; Crux's MCES gives each molecule's atoms
  const rdkit = PackageFunctions.getRdKitModule();
  let mcsMol: RDMol | null = null;
  const mcsGenerated = cachedData.has(pairKey(molecules[0], molecules[1]));
  const mcs = cachedData.get(pairKey(molecules[0], molecules[1])) ?? null;
  try {
    if (mcsGenerated && typeof mcs === 'string')
      mcsMol = rdkit.get_qmol(mcs);
    molecules.forEach((molecule: string, index: number) => {
      const imageHost = ui.canvas(canvasWidth, canvasHeight);
      if (params.seqCol.meta.units === DG.chem.Notation.Smiles) {
        //convert to molFile to draw in coordinates similar to dataframe cell
        molecule = PackageFunctions.convertMolNotation(molecule, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
      }
      const substruct = !mcsGenerated ? null :
        typeof mcs === 'string' ? getUncommonAtomsAndBonds(molecule!, mcsMol, rdkit) : mcs?.[index] ?? null;
      drawMoleculeToCanvas(0, 0, canvasWidth, canvasHeight, imageHost, molecule, '',
        {normalizeDepiction: true, straightenDepiction: true}, substruct);
      ui.empty(hosts[index]);
      if (!cachedData.has(pairKey(molecules[0], molecules[1])))
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
  molIdx: number, idx: number, vertical?: boolean) => void,
  vertical: boolean, flexColNames: boolean, colNameStyle?: any) {
  const element = vertical ? ui.divH([]) : ui.divV([]);
  const columnNames = vertical ? ui.divV([]) : ui.divH([]);
  columnNames.append(colNameStyle ? ui.divText(params.seqCol.name, colNameStyle) : ui.divText(params.seqCol.name));
  columnNames.append(ui.divText(params.activityCol.name));
  columnNames.classList.add(vertical ? 'chem-activity-cliffs-tooltip' : 'chem-activity-cliffs-context-panel');
  if (flexColNames) columnNames.style.display = 'flex';
  element.append(columnNames);
  const hosts: HTMLDivElement[] = [];
  params.points.forEach((molIdx: number, idx: number) => {
    drawFunc(params, element, hosts, molIdx, idx, vertical);
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
  hosts: HTMLDivElement[], molIdx: number, idx: number, vertical?: boolean) {
  const activity = ui.divText(params.activityCol.get(molIdx).toFixed(2));
  activity.style.paddingLeft = '10px';
  const molHost = ui.div([], 'chem-activity-cliffs-molecule');
  molHost.setAttribute('role', 'button');
  molHost.setAttribute('aria-label', `molecule of row ${molIdx + 1}`);
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
  const content = ui.divH([
    molHost,
    activity,
  ]);
  if (!vertical)
    content.classList.add('chem-activity-cliffs-context-panel-content');
  element.append(content);
  hosts.push(molHost);
}
