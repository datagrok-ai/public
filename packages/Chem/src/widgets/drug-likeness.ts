import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full';
import {oclMol, renderDescription} from '../utils/chem-common-ocl';

const _dlp = new OCL.DruglikenessPredictor();

export function assessDruglikeness(molString: string): [number, OCL.IParameterizedString[]] {
  const mol = oclMol(molString);
  return [_dlp.assessDruglikeness(mol), _dlp.getDetail()];
}

export function drugLikenessWidget(molString: string): DG.Widget {
  const [score, description] = assessDruglikeness(molString);
  return new DG.Widget(ui.divV([ui.label(`Score: ${score}`), renderDescription(description)]));
}
