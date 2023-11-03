import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {IWebLogoViewer, WebLogoProps} from '@datagrok-libraries/bio/src/viewers/web-logo';

import {PROPS as wlPROPS} from '../viewers/web-logo-viewer';

import {_package} from '../package';

export class WebLogoApp {
  df: DG.DataFrame;
  view: DG.TableView;

  constructor(
    private readonly urlParams: URLSearchParams,
    private readonly funcName: string,
    private readonly options: Partial<WebLogoProps> = {}
  ) {}

  async init(df: DG.DataFrame): Promise<void> {
    this.df = df;

    await this.buildView();
  }

  // -- View --

  async buildView(): Promise<void> {
    const urlParamsTxt = wu(this.urlParams.entries())
      .map(([key, value]) => `${key}=${encodeURIComponent(value)}`)
      .toArray().join('&');

    this.view = grok.shell.addTableView(this.df);
    this.view.path = this.view.basePath = `func/${_package.name}.${this.funcName}?${urlParamsTxt}`;

    const options: { [p: string]: any } = {...this.options, ...{sequenceColumnName: 'sequence'}};
    for (const [optName, optValue] of this.urlParams.entries()) {
      switch (optName) {
        // boolean
        case wlPROPS.fixWidth:
        case wlPROPS.fitArea:
          options[optName] = ((v) => { return ['1', 'on', 'true'].includes(v.toLowerCase()); })(optValue);
          break;
        default:
          options[optName] = optValue;
      }
    }
    const viewer: DG.Viewer & IWebLogoViewer = (await this.view.dataFrame.plot
      .fromType('WebLogo', options)) as DG.Viewer & IWebLogoViewer;
    this.view.dockManager.dock(viewer, DG.DOCK_TYPE.DOWN, null, 'WebLogo', 0.35);
  }
}
