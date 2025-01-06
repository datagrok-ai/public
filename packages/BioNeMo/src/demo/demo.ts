import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { _package } from '../package';
import { delay } from '@datagrok-libraries/utils/src/test';
import { CONSTANTS, DiffDockModel } from '../diffdock/diffdock-model';
import { DIFFDOCK_HELP, ESMFOLD_HELP, openDataset, showHelpPanel } from './utils';

export async function _demoEsmFoldModel(): Promise<void> {
  const tv: DG.TableView = await openDataset('demo/folding-demo.csv');
  const layoutString = await _package.files.readAsText('demo/folding-demo.layout');
  const layout = DG.ViewLayout.fromJson(layoutString);
  tv.loadLayout(layout);
  tv.dataFrame.currentCell = tv.dataFrame.cell(0, 'Protein');
  showHelpPanel(ESMFOLD_HELP);
}

export async function _demoDiffDockModel(): Promise<void> {
  const tv = await openDataset('demo/docking-demo.csv');
  const layoutString = await _package.files.readAsText('demo/docking-demo.layout');
  const layout = DG.ViewLayout.fromJson(layoutString);
  tv.loadLayout(layout);

  const grid = tv.grid;
  const table = tv.dataFrame;
  const posesColumnName = CONSTANTS.POSES_COLUMN_NAME;
    
  await delay(1000);
    
  const posesColumn = table.columns.byName(posesColumnName);
  posesColumn.setTag(DG.TAGS.SEMTYPE, DG.SEMTYPE.MOLECULE3D);
  posesColumn.setTag(DG.TAGS.CELL_RENDERER, DG.SEMTYPE.MOLECULE3D);
  posesColumn.setTag('docking.role', 'ligand');
    
  await grok.data.detectSemanticTypes(table);
  grid.invalidate();
    
  await delay(100);
    
  const ligandsColName = table.columns.bySemType(DG.SEMTYPE.MOLECULE)?.name;
  table.currentCell = table.cell(0, ligandsColName!);
  showHelpPanel(DIFFDOCK_HELP);
  await delay(1000);
  leaveDiffDockPanelOnly();
}

const waitForButtonAndClick = async (selector: string, timeout: number = 5000) => {
  const startTime = Date.now();
  while (Date.now() - startTime < timeout) {
    const button = document.querySelector(selector) as HTMLElement;
    if (button) {
      button.click();
      return;
    }
    await new Promise(resolve => setTimeout(resolve, 100));
  }
};

const leaveDiffDockPanelOnly = async () => {
  const container = document.querySelectorAll('.panel-content')[1] as HTMLElement;
  const targetPanel = document.querySelector('div.d4-accordion-pane[d4-title="DiffDock"]') as HTMLElement;

  if (!container || !targetPanel) return;
  
  const child = targetPanel.firstChild as HTMLElement;
  if (!child.classList.contains('expanded'))
    child.click();
  
  ui.empty(container);
  container.appendChild(targetPanel);
    
  await waitForButtonAndClick('button.ui-btn.ui-btn-ok[name="button-Run"]');
};
