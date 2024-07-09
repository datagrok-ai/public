import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {getPolyToolEnumerationDialog} from './pt-dialog';

export function polyToolEnumerateUI(cell?: DG.Cell): void {
  getPolyToolEnumerationDialog(cell)
    .then((dialog) => {
      dialog.show({resizable: true});
    })
    .catch((_err: any) => {
      grok.shell.warning('To run PolyTool Enumeration, sketch the macromolecule and select monomers to vary');
    });
}
