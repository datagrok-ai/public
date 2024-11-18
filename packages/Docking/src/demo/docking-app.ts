import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {BaseViewApp} from '@datagrok-libraries/tutorials/src/demo-app';
import { addColorCoding, prepareAutoDockData, runAutodock5, runDocking } from '../package';
import { AutoDockDataType } from '../apps/auto-dock-app';
import { BINDING_ENERGY_COL, POSE_COL, setPose } from '../utils/constants';

export class DockingViewApp extends BaseViewApp {
  protected async processFileData(): Promise<void> {
    const {grid} = this.tableView!;
    const table = this.tableView!.dataFrame;
    const desirableHeight = 100;
    const desirableWidth = 100;
    addColorCoding(grid.columns.byName(BINDING_ENERGY_COL)!);
    grid.sort([BINDING_ENERGY_COL]);
    await grok.data.detectSemanticTypes(table);
    grid.onCellRender.subscribe((args: any) => {
      grid.setOptions({ 'rowHeight': desirableHeight });
      grid.col(POSE_COL)!.width = desirableWidth;
      grid.col(BINDING_ENERGY_COL)!.width = desirableWidth + 50;
    });
    
    await DG.delay(100);
    
    table.col(POSE_COL)!.setTag('docking.role', 'ligand');
    grid.invalidate();
    setPose(POSE_COL);
    
    const autodockResults = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText('System:AppData/Docking/demo_files/autodock_results.csv'));
    const data: AutoDockDataType = await prepareAutoDockData('BACE1', table, 'SMILES', 10);
    //@ts-ignore
    CACHED_DOCKING.K.push(data);
    //@ts-ignore
    CACHED_DOCKING.V.push(autodockResults);
  }

  constructor(parentCall: DG.FuncCall) {
    super(parentCall);
    this.setFormGenerator(this.customFormGenerator);
    this.setFunction = () => this.performDocking();
    this.filePath = 'System:AppData/Docking/demo_files/demo_dataset_small.csv';
    this.browseView.path = 'browse/apps/Admetica'
  }

  private async customFormGenerator(): Promise<HTMLElement> {
    const semValue = DG.SemanticValue.fromGridCell(this.tableView!.grid.cell('smiles', 0));
    const molstarWidget = await runDocking(semValue, 'BACE1', 10);
    molstarWidget!.root.style.width = '100%';
    molstarWidget!.root.style.height = '100%';
    return molstarWidget!.root;
  }

  private async performDocking() {   
    await runAutodock5(
      this.tableView!.dataFrame,
      this.tableView!.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE)!,
      'BACE1',
      10
    );
  }
}