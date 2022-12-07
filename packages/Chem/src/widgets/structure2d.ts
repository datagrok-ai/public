import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {renderMolecule} from '../rendering/render-molecule';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {_convertMolNotation} from '../utils/convert-notation-utils';

export function structure2dWidget(smiles: string): DG.Widget {
  const rdKitModule = getRdKitModule();
  try {
    smiles = _convertMolNotation(smiles, 'unknown', 'smiles', rdKitModule);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule is possible malformed'));
  }
  const width = 200;
  const height = 100;
  // const host = ui.canvas(width, height);
  // drawMoleculeToCanvas(0, 0, width, height, host, smiles);
  const host = renderMolecule(smiles, {renderer: 'RDKit', width: width, height: height});
  return new DG.Widget(host);
}
