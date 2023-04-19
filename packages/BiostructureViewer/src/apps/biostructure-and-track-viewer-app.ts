import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PdbHelper} from '../utils/pdb-helper';
import {IPdbHelper} from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import {_package, getPdbHelper} from '../package';
import {BiostructureAndTrackViewer} from '../viewers/twins/molstar-twin-viewer';

/** The app for biostructure viewer and biotrack viewer */
export class BiostructureAndTrackViewerApp {
  private readonly appFuncName: string = '';

  constructor(appFuncName: string = 'biostructureAndTrackViewerApp') {
    this.appFuncName = appFuncName;
  }

  async init(data?: { df: DG.DataFrame }): Promise<void> {
    if (data) {
      await this.setData(data.df);
    } else {
      const [df] = await BiostructureAndTrackViewerApp.loadData();
      await this.setData(df);
    }
  }

  static async loadData(): Promise<[DG.DataFrame]> {
    const ph: IPdbHelper = await getPdbHelper();
    const pdbStr: string = await _package.files.readAsText('samples/1bdq.pdb');
    const pdbDf: DG.DataFrame = await ph.pdbToDf(pdbStr, '1bdq');
    return [pdbDf];
  }

  async setData(df: DG.DataFrame): Promise<void> {
    this.df = df;

    await this.buildView();
  }

  private df: DG.DataFrame;

  // -- View --

  private view: DG.TableView;

  async buildView(): Promise<void> {
    this.view = grok.shell.addTableView(this.df);
    this.view.path = this.view.basePath = `func/${_package.name}.${this.appFuncName}`;

    const viewer = new BiostructureAndTrackViewer(

    );
    this.view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Biotrack', 0.4);
  }
}
