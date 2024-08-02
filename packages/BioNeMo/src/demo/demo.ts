import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { _package } from '../package';
import { delay } from '@datagrok-libraries/utils/src/test';
import { CONSTANTS } from '../diffdock/diffdock-model';

async function openDataset(name: string): Promise<DG.TableView> {
  const df = await _package.files.readAsText(name);
  const table = DG.DataFrame.fromCsv(df);
  return grok.shell.addTableView(table);
}

export async function _demoEsmFoldModel(): Promise<void> {
  const tv: DG.TableView = await openDataset('demo/folding-demo.csv');
  _package.files.readAsText('demo/folding-demo.layout').then((layoutString: string) => {
    const layout = DG.ViewLayout.fromJson(layoutString);
    tv.loadLayout(layout);
  });
}

export async function _demoDiffDockModel(): Promise<void> {
  const tv: DG.TableView = await openDataset('demo/docking-demo.csv');
  const grid = tv.grid;
  const table = tv.dataFrame;
  await delay(1500);
  table.columns.byName(CONSTANTS.POSES_COLUMN_NAME).setTag(DG.TAGS.SEMTYPE, DG.SEMTYPE.MOLECULE3D);
  table.columns.byName(CONSTANTS.POSES_COLUMN_NAME).setTag(DG.TAGS.CELL_RENDERER, 'xray');
  await grok.data.detectSemanticTypes(table);
  grid.invalidate();
  await delay(2500);
  const layoutString = await _package.files.readAsText('demo/docking-demo.layout');
  const layout = DG.ViewLayout.fromJson(layoutString);
  tv.loadLayout(layout);
}