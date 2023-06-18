import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IWebLogoViewer} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {_package} from '../package';

export class WebLogoApp {
  private _funcName: string = '';

  df: DG.DataFrame;
  view: DG.TableView;

  constructor() {}

  async init(df: DG.DataFrame, funcName: string): Promise<void> {
    this._funcName = funcName;
    this.df = df;

    await this.buildView();
  }

  // -- View --

  async buildView(): Promise<void> {
    this.view = grok.shell.addTableView(this.df);
    this.view.path = this.view.basePath = `func/${_package.name}.${this._funcName}`;

    const viewer: DG.Viewer & IWebLogoViewer = (await this.view.dataFrame.plot.fromType('WebLogo', {
      sequenceColumnName: 'sequence',
    }));
    this.view.dockManager.dock(viewer, DG.DOCK_TYPE.DOWN, null, 'WebLogo', 0.35);
  }
}
