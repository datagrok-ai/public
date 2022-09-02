import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {findMCS} from '../scripts-api';
import {drawMoleculeToCanvas} from '../utils/chem-common-rdkit';
import { ITooltipAndPanelParams } from '@datagrok-libraries/ml/src/viewers/activity-cliffs';

const canvasWidth = 200;
const canvasHeight = 100;

export async function findMcsAndUpdateDrawings(params: ITooltipAndPanelParams) {
  if (!params.cashedData[params.line.id]) {
    drawMoleculesWithMcsAsync(params); 
  }
  drawMolecules(params);
}

async function drawMoleculesWithMcsAsync(params: ITooltipAndPanelParams){
  const mcsDf = DG.DataFrame.create(2);
  mcsDf.columns.addNewString('smiles').init((i) => params.seqCol.get(params.line.mols[i]));
  const mcs = await findMCS('smiles', mcsDf, true);
  params.cashedData[params.line.id] = mcs;
  drawMolecules(params);
}

function drawMolecules(params: ITooltipAndPanelParams) {
  params.line.mols.forEach((mol: number, index: number) => {
    const imageHost = ui.canvas(canvasWidth, canvasHeight);
    const r = window.devicePixelRatio;
    imageHost.width = canvasWidth * r;
    imageHost.height = canvasHeight * r;
    imageHost.style.width = (canvasWidth).toString() + 'px';
    imageHost.style.height = (canvasHeight).toString() + 'px';
    drawMoleculeToCanvas(0, 0, canvasWidth, canvasHeight, imageHost, params.seqCol.get(mol), params.cashedData[params.line.id]);
    ui.empty(params.hosts[index]);
    if (!params.cashedData[params.line.id]) {
      params.hosts[index].append(ui.divText('MCS loading...'));
    }
    params.hosts[index].append(imageHost);
  });
}
