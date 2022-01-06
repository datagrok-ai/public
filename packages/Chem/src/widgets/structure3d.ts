import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {smilesTo3DCoordinates} from '../scripts-api';

export async function structure3dWidget(smiles: string) {
  const sdf = (await smilesTo3DCoordinates(smiles)).replaceAll('\\n', '\n');
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
