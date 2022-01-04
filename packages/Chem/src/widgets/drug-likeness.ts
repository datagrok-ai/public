import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';
import {renderDescription} from '../chem_common_ocl';

export function oclMol(mol: string): OCL.Molecule {
  if (mol.includes('M  END'))
    return OCL.Molecule.fromMolfile(mol);
  else
    return OCL.Molecule.fromSmiles(mol);
}

export function drugLikenessWidget(molString: string) {
  const mol = oclMol(molString);
  const dlp = new OCL.DruglikenessPredictor();
  return new DG.Widget(ui.divV([
    ui.label(`Score: ${dlp.assessDruglikeness(mol)}`),
    renderDescription(dlp.getDetail()),
  ]));
}
