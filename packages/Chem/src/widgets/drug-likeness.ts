import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full';
import {oclMol, renderDescription} from '../utils/chem-common-ocl';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {getRdKitModule} from '../package';

const _dlp = new OCL.DruglikenessPredictor();

export function assessDruglikeness(molString: string): [number, OCL.IParameterizedString[]] {
  const mol = oclMol(molString);
  return [_dlp.assessDruglikeness(mol), _dlp.getDetail()];
}

export function drugLikenessWidget(molString: string): DG.Widget {
  const rdKitModule = getRdKitModule();
  try {
    molString = _convertMolNotation(molString, 'unknown', 'smiles', rdKitModule);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule is possible malformed'));
  }
  let score: number;
  let description: OCL.IParameterizedString[];
  try {
    [score, description] = assessDruglikeness(molString);
  } catch (e) {
    return new DG.Widget(ui.divText('Could not asses drug likeness'));
  }
  return new DG.Widget(ui.divV([ui.label(`Score: ${score}`), renderDescription(description)]));
}
