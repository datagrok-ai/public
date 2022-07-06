import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {findMCS} from '../scripts-api';
import {drawMoleculeToCanvas} from '../utils/chem-common-rdkit';
import { IDrawTooltipParams } from '@datagrok-libraries/ml/src/viewers/activity-cliffs';

export function drawTooltip(params: IDrawTooltipParams) {
  params.tooltips[params.line.id] = ui.divText('Loading...');
  const mcsDf = DG.DataFrame.create(2);
  mcsDf.columns.addNewString('smiles').init((i) => params.df.get(params.seqCol.name, params.line.mols[i]));
  findMCS('smiles', mcsDf, true).then((mcs) => {
    params.tooltips[params.line.id] = ui.divH([]);
    const columnNames = ui.divV([
      ui.divText('smiles'),
      ui.divText(params.activity.name),
    ]);
    columnNames.style.fontWeight = 'bold';
    columnNames.style.display = 'flex';
    columnNames.style.justifyContent = 'space-between';
    params.tooltips[params.line.id].append(columnNames);
    params.line.mols.forEach((mol: any) => {
      const imageHost = ui.canvas(200, 100);
      drawMoleculeToCanvas(0, 0, 200, 100, imageHost, params.df.get('smiles', mol), mcs);
      const activity = ui.divText(params.df.get(params.activity.name, mol).toFixed(2));
      activity.style.display = 'flex';
      activity.style.justifyContent = 'left';
      activity.style.paddingLeft = '30px';
      params.tooltips[params.line.id].append(ui.divV([
        imageHost,
        activity,
      ]));
    });
    ui.tooltip.show(params.tooltips[params.line.id], params.x, params.y);
  });
}