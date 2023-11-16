import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {Molecule3DUnitsHandler} from '@datagrok-libraries/bio/src/molecule-3d';

import {defaultErrorHandler} from './err-info';
import {openPdbResidues} from '../package';

import {_package} from '../package';

export function addContextMenuUI(event: DG.EventData): void {
  try {
    if (event.args.item) {
      const item = event.args.item;
      // TODO: TreeViewNode.value is not real DG.FileInfo (no extension property)
      // if ((item instanceof DG.TreeViewNode))
      //   item = item.value;
      const menu: DG.Menu = event.args.menu;
      if (item instanceof DG.FileInfo || item.constructor.name == 'FileInfo') {
        const fi: DG.FileInfo = item;
        if (addContextMenuForFileInfo(fi, menu))
          event.preventDefault();
      } else if (item instanceof DG.GridCell || item.constructor.name == 'GridCell') {
        const gridCell: DG.GridCell = item;
        if (addContextMenuForCell(gridCell, menu))
          event.preventDefault();
      }
    }
  } catch (err: any) {
    defaultErrorHandler(err);
    throw err;
  }
}

export function addContextMenuForFileInfo(fi: DG.FileInfo, menu: DG.Menu): boolean {
  switch (fi.extension.toLowerCase()) {
    case 'pdb': {
      menu.item('Open table residues', async () => {
        openPdbResidues(fi)
          .catch(defaultErrorHandler);
      });
      return true;
    }
  }
  return false;
}

export function addContextMenuForCell(gridCell: DG.GridCell, menu: DG.Menu): boolean {
  function copyRawValue(): void {
    const tgtValue: string = gridCell.cell.value;
    if (!navigator.clipboard) {
      //
      grok.shell.warning('The clipboard functionality requires a secure origin, either HTTPS or localhost.');
    } else {
      navigator.clipboard.writeText(tgtValue)
        .catch(defaultErrorHandler);
      grok.shell.info(`Value copied to clipboard`);
    }
  }

  function downloadRawValue(): void {
    const tableCol = gridCell.tableColumn!;
    const uh = Molecule3DUnitsHandler.getOrCreate(tableCol);
    const fileName: string = uh.getFileName(gridCell.cell);
    const tgtValue: string = gridCell.cell.value;
    DG.Utils.download(fileName, tgtValue);
  }

  switch (gridCell.tableColumn!.semType) {
    case DG.SEMTYPE.MOLECULE3D: {
      menu.item('Copy', copyRawValue);
      menu.item('Download', downloadRawValue);
      return true;
    }
  }
  return false;
}
