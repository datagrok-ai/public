import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {defaultErrorHandler} from './err-info';
import {openMonomerLibrary} from './open-monomer-library';

export type HelmWindowType = Window & {
  $helm?: {
    contextMenuError?: any;
  }
}
declare const window: HelmWindowType;

export function addContextMenuForFileInfoJsonUI(event: DG.EventData): void {
  try {
    const fi: DG.FileInfo = event.args.item;
    const menu: DG.Menu = event.args.menu;
    menu.item('Open monomer library', async () => {
      openMonomerLibrary(fi)
        .catch(defaultErrorHandler);
    });
  } catch (err: any) {
    defaultErrorHandler(err);
    if (!window.$helm) window.$helm = {};
    window.$helm.contextMenuError = err;
  }
}
