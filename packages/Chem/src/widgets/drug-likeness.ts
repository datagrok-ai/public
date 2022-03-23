import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';
import {oclMol, renderDescription} from '../utils/chem-common-ocl';

const _dlp = new OCL.DruglikenessPredictor();

export function drugLikenessWidget(molString: string): DG.Widget {
  const mol = oclMol(molString);
  return new DG.Widget(ui.divV([
    ui.label(`Score: ${_dlp.assessDruglikeness(mol)}`),
    renderDescription(_dlp.getDetail(), molString),
  ]));
}
