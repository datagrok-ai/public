import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full';
import {oclMol, renderDescription} from '../utils/chem-common-ocl';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {getRdKitModule} from '../package';
import {OCLService} from '../open-chem/ocl-service';
import '../../css/chem.css';

const _dlp = new OCL.DruglikenessPredictor();

export function assessDruglikeness(molString: string): [number, OCL.IParameterizedString[]] {
  const mol = oclMol(molString);
  return [_dlp.assessDruglikeness(mol), _dlp.getDetail()];
}

export function drugLikenessWidget(molString: DG.SemanticValue): DG.Widget {
  const rdKitModule = getRdKitModule();
  let convertedMolstring: string = '';
  try {
    convertedMolstring = _convertMolNotation(molString.value,
      DG.chem.Notation.Unknown, DG.chem.Notation.Smiles, rdKitModule);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  }
  let score: number;
  let description: OCL.IParameterizedString[];
  try {
    [score, description] = assessDruglikeness(convertedMolstring);
  } catch (e) {
    return new DG.Widget(ui.divText('Could not asses drug likeness'));
  }
  const descriptionHost = renderDescription(description);
  descriptionHost.style.overflow = 'hidden';
  descriptionHost.style.maxHeight = '400px';
  const calcFowWholeButton = ui.icons.add(async () => {
    const pi = DG.TaskBarProgressIndicator.create('Calculating drug likeness scores');
    try {
      const oclService = new OCLService();
      const res = await oclService.getDrugLikeness(molString.cell.column);
      oclService.terminate();
      Object.keys(res).forEach((k) => {
        molString.cell.dataFrame.columns.add(DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, k, res[k]));
      });
    } catch (e) {
      console.error(e);
    } finally {
      pi.close();
    }
  }, 'Calculate drug likeness for the whole table');
  const scoreDiv = ui.divH([ui.label(`Score: ${score.toFixed(2)}`), calcFowWholeButton],
    {classes: 'chem-drug-likeness-widget-score-div'});
  calcFowWholeButton.classList.add('chem-drug-likeness-calc-all-button');

  const host = ui.divV([scoreDiv,
    ui.label(` ${description[0].value}`), descriptionHost], {classes: 'ui-box'});
  ui.tools.setHoverVisibility(host, [calcFowWholeButton]);

  return new DG.Widget(host);
}
