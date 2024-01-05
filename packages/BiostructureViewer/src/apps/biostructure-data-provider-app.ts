import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';

export type BiostructureDataProviderAppData = { df: DG.DataFrame, idColumnName: string, dataProviderNqName: string };

export class BiostructureDataProviderApp {
  private readonly appName: string;

  constructor(appName: string = 'biostructureDataProviderApp') {
    this.appName = appName;
  }

  async init(data?: BiostructureDataProviderAppData): Promise<void> {
    if (!data)
      data = await BiostructureDataProviderApp.loadData();
    await this.setData(data);
  }

  static async loadData(): Promise<BiostructureDataProviderAppData> {
    const fn = 'pdb_id.csv';
    const df = await _package.files.readCsv(fn);
    return {df: df, idColumnName: 'pdb_id', dataProviderNqName: `${_package.name}:getBiostructureRcsbPdb`};
  }

  private data: BiostructureDataProviderAppData;

  async setData(data: BiostructureDataProviderAppData): Promise<void> {
    this.data = data;
    await this.buildView();
  }

  // -- View --
  private view?: DG.TableView;

  async buildView(): Promise<void> {
    this.view = grok.shell.addTableView(this.data.df);
    this.view.path = this.view.basePath = `/func/${_package.name}.${this.appName}`;

    const viewer = await this.data.df.plot.fromType('Biostructure', {
      biostructureIdColumnName: this.data.idColumnName,
      biostructureDataProvider: this.data.dataProviderNqName,
    });
    this.view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Biostructure', 0.4);
    grok.shell.windows.showContextPanel = true;
    grok.shell.o = viewer;
  }
}
