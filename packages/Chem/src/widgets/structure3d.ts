import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {smilesTo3DCoordinates} from '../scripts-api';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {_convertMolNotation} from '../utils/convert-notation-utils';

export async function structure3dWidget(molecule: string): Promise<DG.Widget> {
  const rdKitModule = getRdKitModule();
  try {
    molecule = _convertMolNotation(molecule, 'unknown', 'molblock', rdKitModule);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  }
  let sdf: string;
  try {
    sdf = (await smilesTo3DCoordinates(molecule)).replaceAll('\\n', '\n');
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule has no atoms'));
  }
  const stringBlob = new Blob([sdf], {type: 'text/plain'});

  const nglHost = ui.div([], {classes: 'd4-ngl-viewer', id: 'ngl-3d-host'});
  nglHost.style.width = '300px';
  nglHost.style.height = '300px';
  nglHost.style.backgroundColor = 'white';

  //@ts-ignore
  const stage = new NGL.Stage(nglHost, {backgroundColor: 'white'});
  //@ts-ignore
  stage.loadFile(stringBlob, {ext: 'sdf'}).then(function(comp: NGL.StructureComponent) {
    stage.setSize(300, 300);
    comp.addRepresentation('ball+stick');
    comp.autoView();
  });

  return new DG.Widget(nglHost);
}
