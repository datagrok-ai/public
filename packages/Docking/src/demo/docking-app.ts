import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {BaseViewApp} from '@datagrok-libraries/tutorials/src/demo-base-view';
import { getAutodockSingle, prepareAutoDockData, runDocking } from '../package';
import { CACHED_DOCKING } from '../utils/constants';
import { prepareDockingData } from './demo-docking';

export class DockingViewApp extends BaseViewApp {
  constructor(parentCall: DG.FuncCall) {
    super(parentCall);

    this.setFormGenerator(this.generateCustomForm);
    this.filePath = 'System:AppData/Docking/demo_files/demo_dataset_small.csv';
    this.browseView.path = 'browse/apps/Admetica';
    this.uploadCachedData = this.loadCachedAutodockData.bind(this);
    this.addTabControl = false;
  }

  private async generateCustomForm(): Promise<HTMLElement> {
    const smilesCell = this.tableView!.grid.cell('smiles', 0);
    const semValue = DG.SemanticValue.fromGridCell(smilesCell);

    const molstarWidget = await runDocking(semValue, 'BACE1', 10);
    this.applyFullSizeStyles(molstarWidget?.root?.firstElementChild as HTMLElement);

    return molstarWidget!.root;
  }

  private async loadCachedAutodockData(): Promise<HTMLElement> {
    const datasetPath = 'System:AppData/Docking/demo_files/demo_app_sample.csv';
    const csvText = await grok.dapi.files.readAsText(datasetPath);
    const sampleData = DG.DataFrame.fromCsv(csvText);

    const poseColumn = sampleData.getCol('pose');
    poseColumn.semType = DG.SEMTYPE.MOLECULE3D;
    poseColumn.setTag(DG.TAGS.UNITS, 'pdb');

    this.tableView!.dataFrame = sampleData;
    await grok.data.detectSemanticTypes(sampleData);
    this.tableView!.grid.invalidate();

    const resultsPath = 'System:AppData/Docking/demo_files/demo_app_results.csv';
    await prepareDockingData(sampleData, resultsPath, 'BACE1', 'smiles', 10);

    const dockingWidget = await this.createDockingWidget(sampleData.cell(0, 'pose'), sampleData);
    return dockingWidget!;
  }

  private async createDockingWidget(cell: DG.Cell, dataFrame: DG.DataFrame): Promise<HTMLElement> {
    const semValue = DG.SemanticValue.fromTableCell(cell);
    const widgetRoot = (await getAutodockSingle(semValue, false, dataFrame))?.root;

    if (widgetRoot)
      this.applyFullSizeStyles(widgetRoot.firstElementChild as HTMLElement);

    return widgetRoot!;
  }

  private applyFullSizeStyles(element: HTMLElement | null): void {
    if (!element) return;

    element.style.width = '100%';
    element.style.height = '100%';
  }
}
