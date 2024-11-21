import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {BaseViewApp} from '@datagrok-libraries/tutorials/src/demo-base-view';
import { runAutodock5, runDocking } from '../package';

export class DockingViewApp extends BaseViewApp {
  constructor(parentCall: DG.FuncCall) {
    super(parentCall);
    this.setFormGenerator(this.customFormGenerator);
    this.setFunction = () => this.performDocking();
    this.filePath = 'System:AppData/Docking/demo_files/demo_dataset_small.csv';
    this.browseView.path = 'browse/apps/Admetica';
    this.addTabControl = false;
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