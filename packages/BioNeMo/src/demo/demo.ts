import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { _package } from '../package';
import { delay } from '@datagrok-libraries/utils/src/test';
import { CONSTANTS, DiffDockModel } from '../diffdock/diffdock-model';
import { DIFFDOCK_HELP, ESMFOLD_HELP, openDataset, showHelpPanel } from './utils';

export async function _demoEsmFoldModel(): Promise<void> {
  const tv: DG.TableView = await openDataset('demo/folding-demo.csv');
  _package.files.readAsText('demo/folding-demo.layout').then((layoutString: string) => {
    const layout = DG.ViewLayout.fromJson(layoutString);
    tv.loadLayout(layout);
  });
  tv.dataFrame.currentCell = tv.dataFrame.cell(0, 'Protein');
  showHelpPanel(ESMFOLD_HELP);
}

export async function _demoDiffDockModel(): Promise<void> {
  const tv = await openDataset('demo/docking-demo.csv');
  const grid = tv.grid;
  const table = tv.dataFrame;
  const posesColumnName = CONSTANTS.POSES_COLUMN_NAME;

  await delay(1000);

  const posesColumn = table.columns.byName(posesColumnName);
  posesColumn.setTag(DG.TAGS.SEMTYPE, DG.SEMTYPE.MOLECULE3D);
  posesColumn.setTag(DG.TAGS.CELL_RENDERER, 'xray');
  posesColumn.setTag('docking.role', 'ligand');

  await grok.data.detectSemanticTypes(table);
  grid.invalidate();

  await delay(100);

  const layoutString = await _package.files.readAsText('demo/docking-demo.layout');
  const layout = DG.ViewLayout.fromJson(layoutString);
  tv.loadLayout(layout);

  const ligandsColName = table.columns.bySemType(DG.SEMTYPE.MOLECULE)?.name;
  table.currentCell = table.cell(0, ligandsColName!);
  showHelpPanel(DIFFDOCK_HELP);
}