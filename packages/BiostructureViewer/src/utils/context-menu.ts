import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {Molecule3DUnitsHandler} from '@datagrok-libraries/bio/src/molecule-3d';
import {BiostructureData, BiostructureDataJson} from '@datagrok-libraries/bio/src/pdb/types';

import {defaultErrorHandler} from './err-info';

import {_package} from '../package';

export type BiostructureViewerWindowType = Window & {
  $biostructureViewer?: {
    contextMenuError?: any,
  },
};

export function copyRawValue(gridCell: DG.GridCell): void {
  const tgtValue: string = gridCell.cell?.value ?? '';
  if (!navigator.clipboard) {
    //
    grok.shell.warning('The clipboard functionality requires a secure origin, either HTTPS or localhost.');
  } else {
    navigator.clipboard.writeText(tgtValue)
      .catch(defaultErrorHandler).then(() => grok.shell.info(`Value copied to clipboard`));
  }
}

export function downloadRawValue(gridCell: DG.GridCell): void {
  try {
    const tableCol = gridCell.tableColumn!;
    const mh = Molecule3DUnitsHandler.getOrCreate(tableCol);
    const fileName: string = mh.getFileName(gridCell.cell);
    const tgtValue: string = gridCell.cell.value;
    DG.Utils.download(fileName, tgtValue);
  } catch (err: any) {
    defaultErrorHandler(err);
  }
}

export async function showBiostructureViewer(gridCell: DG.GridCell): Promise<void> {
  try {
    const tableCol = gridCell.tableColumn!;
    const uh = Molecule3DUnitsHandler.getOrCreate(tableCol);
    const valueData: BiostructureData = {
      binary: false,
      data: gridCell.cell.value,
      ext: uh.fileExt(),
    };
    if (gridCell.grid.view instanceof DG.TableView) {
      const view = gridCell.grid.view as DG.TableView;
      const viewer = await gridCell.tableColumn?.dataFrame.plot.fromType('Biostructure', {
        dataJson: BiostructureDataJson.fromData(valueData),
      });
      if (viewer) view.dockManager.dock(viewer.root);
    }
  } catch (err: any) {
    defaultErrorHandler(err);
  }
}

export async function showNglViewer(gridCell: DG.GridCell): Promise<void> {
  try {
    const tableCol = gridCell.tableColumn!;
    const mh = Molecule3DUnitsHandler.getOrCreate(tableCol);
    const valueData: BiostructureData = {
      binary: false,
      data: gridCell.cell.value,
      ext: mh.fileExt(),
    };
    if (gridCell.grid.view instanceof DG.TableView) {
      const view = gridCell.grid.view as DG.TableView;
      const viewer = await gridCell.tableColumn?.dataFrame.plot.fromType('NGL', {
        dataJson: BiostructureDataJson.fromData(valueData),
      });
      if (viewer) view.dockManager.dock(viewer.root);
    }
  } catch (err: any) {
    defaultErrorHandler(err);
  }
}
