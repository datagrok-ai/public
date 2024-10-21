import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {defaultErrorHandler} from './err-info';
import {polyToolEnumerateHelmUI} from '../polytool/pt-enumeration-helm-dialog';
import {polyToolEnumerateChemUI} from '../polytool/pt-dialog';

import {_package} from '../package';

export type SequenceTranslatorWindowType = Window & {
  $sequenceTranslator?: {
    contextMenuError?: any
  },
};
declare const window: SequenceTranslatorWindowType;

export function addContextMenuUI(event: DG.EventData): void {
  try {
    const item = event.args.item;
    if (item) {
      const menu: DG.Menu = event.args.menu;
      if (item instanceof DG.GridCell || item.constructor.name == 'GridCell') {
        const gridCell: DG.GridCell = item;
        if (addContextMenuForCell(gridCell, menu))
          event.preventDefault();
      }
    }
  } catch (err: any) {
    defaultErrorHandler(err);
    if (!window.$sequenceTranslator) window.$sequenceTranslator = {};
    window.$sequenceTranslator.contextMenuError = err;
  }
}

function addContextMenuForCell(gridCell: DG.GridCell, menu: DG.Menu): boolean {
  const logPrefix: string = `ST: addContextMenuForCell()`;
  _package.logger.debug(`${logPrefix}, start`);

  if (gridCell && gridCell.tableColumn) {
    switch (gridCell.tableColumn.semType) {
    case DG.SEMTYPE.MACROMOLECULE: {
      menu.item('PolyTool-Enumerate', () => { polyToolEnumerateHelmUI(gridCell.cell); });
      return true;
    }
    case DG.SEMTYPE.MOLECULE: {
      menu.item('PolyTool-Enumerate', () => { polyToolEnumerateChemUI(gridCell.cell); });
      return true;
    }
    }
  }
  return false;
}
