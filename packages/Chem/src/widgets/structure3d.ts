import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {smilesTo3DCoordinates} from '../scripts-api';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {_convertMolNotation} from '../utils/convert-notation-utils';


export async function structure3d(molecule: string, width: number, height: number): Promise<HTMLElement> {
  let sdf: string;
  try {
    sdf = (await smilesTo3DCoordinates(molecule)).replaceAll('\\n', '\n');
  } catch (e) {
    return ui.divText('Molecule has no atoms');
  }
  const stringBlob = new Blob([sdf], {type: 'text/plain'});

  const nglHost = ui.div([], {classes: 'd4-ngl-viewer', id: 'ngl-3d-host'});
  nglHost.style.width = `${width}px`;
  nglHost.style.height = `${height}px`;
  nglHost.style.backgroundColor = 'white';

  //@ts-ignore
  const stage = new NGL.Stage(nglHost, {backgroundColor: 'white'});
  //@ts-ignore
  stage.loadFile(stringBlob, {ext: 'sdf'}).then(function(comp: NGL.StructureComponent) {
    stage.setSize(width, height);
    comp.addRepresentation('ball+stick');
    comp.autoView();
  });
  return nglHost;
}