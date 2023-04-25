import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {_package} from '../package';
import {NglViewerApp} from '../apps/ngl-viewer-app';
import {TaskBarProgressIndicator} from 'datagrok-api/dg';

export async function dockingDemoApp(appName: string, pi: TaskBarProgressIndicator): Promise<void> {
  const piMsg = pi.description;
  const pdbStr: string = await _package.files.readAsText('samples/protease.pdb').then((value: string) => {
    pi.update(33, piMsg);
    return value;
  });
  const sdfBytes: Uint8Array = await _package.files.readAsBytes('samples/1bdq.sdf').then((value: Uint8Array) => {
    pi.update(66, piMsg);
    return value;
  });
  const ligands: DG.DataFrame = (await grok.functions.call(
    'Chem:importSdf', {bytes: sdfBytes}))[0];

  const app = new NglViewerApp(appName);
  await app.init({ligands: ligands, macromolecule: pdbStr}).then(() => {
    pi.update(100, piMsg);
  });
}
