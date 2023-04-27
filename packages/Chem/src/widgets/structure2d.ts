import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {renderMolecule} from '../rendering/render-molecule';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {_convertMolNotation} from '../utils/convert-notation-utils';

const WIDTH = 200;
const HEIGHT = 100;
const MORE_ICON_FONT_WEIGHT = '500';

export function structure2dWidget(molecule: string): DG.Widget {
  const rdKitModule = getRdKitModule();
  let molblock: string;
  try {
    molblock = _convertMolNotation(molecule, DG.chem.Notation.Unknown, DG.chem.Notation.MolBlock, rdKitModule);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  }

  return new DG.Widget(get2dMolecule(molecule));
}

function get2dMolecule(molecule: string): HTMLElement{
  const molecule2d = renderMolecule(molecule, { renderer: 'RDKit', width: WIDTH, height: HEIGHT });
  const moreIcon = molecule2d.querySelector('.chem-mol-view-icon.pep-more-icon') as HTMLElement;
  moreIcon.style.left = `${WIDTH}px`;
  moreIcon.style.fontWeight = MORE_ICON_FONT_WEIGHT;
  moreIcon.classList.remove('pep-more-icon');
  return molecule2d;
}
