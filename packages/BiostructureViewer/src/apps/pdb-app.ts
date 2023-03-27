import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PdbHelper} from '../utils/pdb-helper';
import {IPdbHelper} from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import {getPdbHelper} from '../package';

/** The app for .pdb file handler */
export class PdbApp {
  private _funcName: string = '';

  constructor() {}

  /** {@link df} created with pdbToDf() */
  async init(df: DG.DataFrame, funcName: string = 'pdbApp'): Promise<void> {
    this._funcName = funcName;
    await this.loadData(df);
  }

  async loadData(df: DG.DataFrame): Promise<void> {
    const ph: IPdbHelper = await getPdbHelper();
    await this.setData(df);
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

    const viewer: DG.JsViewer = (await this.view.dataFrame.plot.fromType('Biostructure', {})) as DG.JsViewer;
    this.view.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'NGL', 0.4);
  }
}
