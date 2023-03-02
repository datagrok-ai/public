import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package';
import {Unsubscribable} from 'rxjs';

export class NglForGridTestApp {
  private readonly appName: string;
  df: DG.DataFrame;

  constructor(appName: string = 'nglForGridTestApp') {
    this.appName = appName;
  }

  async init(data?: { df: DG.DataFrame }): Promise<void> {
    if (data) {
      await this.setData(data.df);
    } else {
      const [df] = await NglForGridTestApp.loadData();
      await this.setData(df);
    }
  }

  static async loadData(): Promise<[DG.DataFrame]> {
    const dfCsv: string = await _package.files.readAsText('pdb_data.csv');
    const df: DG.DataFrame = DG.DataFrame.fromCsv(dfCsv);
    return [df];
  }

  async setData(df: DG.DataFrame): Promise<void> {
    this.df = df;

    await this.buildView();
  }

  // -- View --

  private view: DG.TableView;
  private viewSubs: Unsubscribable[] = [];

  async buildView(): Promise<void> {
    this.view = grok.shell.addTableView(this.df);
    this.view.path = this.view.basePath = `/func/${_package.name}.${this.appName}`;
  }
}
