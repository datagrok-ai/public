import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {renderMolecule} from '../rendering/render-molecule';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import { getMolSafe } from '../utils/mol-creation_rdkit';

const WIDTH = 200;
const HEIGHT = 100;
const MORE_ICON_FONT_WEIGHT = '500';

export function structure2dWidget(molecule: string): DG.Widget {
  const rdKitModule = getRdKitModule();
  let mol = getMolSafe(molecule, {}, rdKitModule).mol;
  const resultWidget = mol && mol.is_valid() ? new DG.Widget(get2dMolecule(molecule)) : new DG.Widget(ui.divText('Molecule is possibly malformed'));
  mol?.delete();
  return resultWidget;
}

function get2dMolecule(molecule: string): HTMLElement{
  const molecule2d = renderMolecule(molecule, { renderer: 'RDKit', width: WIDTH, height: HEIGHT });
  const moreIcon = molecule2d.querySelector('.chem-mol-view-icon.pep-more-icon') as HTMLElement;
  moreIcon.style.left = `${WIDTH}px`;
  moreIcon.style.fontWeight = MORE_ICON_FONT_WEIGHT;
  moreIcon.classList.remove('pep-more-icon');
  return molecule2d;
}
