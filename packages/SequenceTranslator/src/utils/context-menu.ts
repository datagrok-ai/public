import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {defaultErrorHandler} from './err-info';
import {_package, addContextMenu} from '../package';
import {SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {polyToolEnumerateUI} from '../polytool/pt-ui';

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

  const polyToolEnumerate = () => {
    polyToolEnumerateUI(gridCell.cell);
  };

  if (gridCell && gridCell.tableColumn) {
    switch (gridCell.tableColumn.semType) {
    case DG.SEMTYPE.MACROMOLECULE: {
      const sh = SeqHandler.forColumn(gridCell.tableColumn);
      if (sh.notation === NOTATION.HELM) {
        menu
          .item('PolyTool-Enumerate', polyToolEnumerate);
        true;
      }
    }
    }
  }
  return false;
}
