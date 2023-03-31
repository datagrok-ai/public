import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {renderMolecule} from '../rendering/render-molecule';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import { structure3d } from './structure3d';
import { identifiers } from './identifiers';

const WIDTH = 200;
const HEIGHT = 100;
const WIDTH3D = 300;
const HEIGHT3D = 300;
const MORE_ICON_FONT_WEIGHT = '500';

export function moleculeOverviewWidget(molecule: string): DG.Widget {
  const rdKitModule = getRdKitModule();
  let molblock: string;
  try {
    molblock = _convertMolNotation(molecule, DG.chem.Notation.Unknown, DG.chem.Notation.MolBlock, rdKitModule);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  }

  const acc = ui.accordion('chem-molecule-overview');
  acc.addPane(`Identifiers`, () => ui.wait(async () => await identifiers(molblock)));
  acc.addPane(`3D Structure`, () => ui.wait(async () => await structure3d(molecule, WIDTH3D, HEIGHT3D)));

  return new DG.Widget(ui.divV([
    get2dMolecule(molecule),
    acc.root
  ]));
}

function get2dMolecule(molecule: string): HTMLElement{
  const molecule2d = renderMolecule(molecule, { renderer: 'RDKit', width: WIDTH, height: HEIGHT });
  const moreIcon = molecule2d.querySelector('.chem-mol-view-icon.pep-more-icon') as HTMLElement;
  moreIcon.style.left = `${WIDTH}px`;
  moreIcon.style.fontWeight = MORE_ICON_FONT_WEIGHT;
  moreIcon.classList.remove('pep-more-icon');
  return molecule2d;
}
