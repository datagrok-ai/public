import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {BaseViewApp} from '@datagrok-libraries/tutorials/src/demo-base-view';

import { getAutodockSingle, runDocking } from '../package';

export class DockingViewApp extends BaseViewApp {
  protected STORAGE_NAME: string = 'docking-sketcher-values';

  constructor(parentCall: DG.FuncCall) {
    super(parentCall);
    this.sketcherValue = {'(5~{S})-5-(1~{H}-indol-2-yl)pyrrolidin-2-one': 'O=C(CC1)N[C@@H]1c1cc2ccccc2[nH]1'};
    this.addTabControl = false;
    this.tableName = 'Docking';
  }

  protected cached(): boolean {
    return this.target?.stringValue === 'kras';
  }

  protected async customInit(): Promise<void> {
    const items = await grok.functions.call('Docking:getConfigFiles');
    const helpIcon = ui.icons.help(() => {
      grok.shell.windows.showHelp = true;
      grok.shell.windows.help.showHelp('/help/develop/domains/chem/docking');
    });
    this.target = ui.input.choice('Target', {value: 'kras', items: items, onValueChanged: async () => {
      await this.onChanged(this.sketcherInstance.getSmiles());
    }});
    this.target.root.classList.add('demo-target-root');
    const container = ui.divH([this.target.root, helpIcon]);
    container.style.cssText = 'overflow: visible !important; gap: 10px; align-items: baseline;';
    this.formContainer.insertBefore(container, this.formContainer.firstChild);
  }

  protected async formGenerator(): Promise<HTMLElement | null> {
    const smilesCell = this.tableView!.grid.cell('smiles', 0);
    const semValue = DG.SemanticValue.fromGridCell(smilesCell);

    const molstarWidget = await runDocking(semValue, this.target!.stringValue, 10);
    this.applyFullSizeStyles(molstarWidget?.root?.firstElementChild as HTMLElement);

    if (molstarWidget) {
      return molstarWidget.root;
    }
    return null;
  }

  protected async uploadCachedData(): Promise<HTMLElement> {
    const datasetPath = 'System:AppData/Docking/demo_files/demo_app_sample.csv';
    const csvText = await grok.dapi.files.readAsText(datasetPath);
    const sampleData = DG.DataFrame.fromCsv(csvText);

    const poseColumn = sampleData.getCol('pose');
    poseColumn.semType = DG.SEMTYPE.MOLECULE3D;
    poseColumn.setTag(DG.TAGS.UNITS, 'pdb');

    this.tableView!.dataFrame = sampleData;
    await grok.data.detectSemanticTypes(sampleData);
    this.tableView!.grid.invalidate();

    const dockingWidget = await this.createDockingWidget(sampleData.cell(0, 'pose'), sampleData);
    return dockingWidget!;
  }

  private async createDockingWidget(cell: DG.Cell, dataFrame: DG.DataFrame): Promise<HTMLElement | null> {
    const semValue = DG.SemanticValue.fromTableCell(cell);
    const widget = (await getAutodockSingle(semValue, false, dataFrame));

    if (widget) {
      this.applyFullSizeStyles(widget.root.firstElementChild as HTMLElement);
      return widget.root;
    }
    return null;
  }

  private applyFullSizeStyles(element: HTMLElement | null): void {
    if (!element) return;

    element.style.width = '100%';
    element.style.height = '100%';
  }
}