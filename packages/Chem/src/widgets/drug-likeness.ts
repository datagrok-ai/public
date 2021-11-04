import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';
import {renderDescription} from '../chem_common';

export function drugLikenessWidget(smiles: string) {
  const mol = OCL.Molecule.fromSmiles(smiles);
  const dlp = new OCL.DruglikenessPredictor();
  return new DG.Widget(ui.divV([
    ui.label(`Score: ${dlp.assessDruglikeness(mol)}`),
    renderDescription(dlp.getDetail()),
  ]));
}
