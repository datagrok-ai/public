/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {_package} from '../../../package';
import {AppUIBase} from './app-ui-base';

export abstract class IsolatedAppUIBase extends AppUIBase {
  constructor(appName: string) {
    super(appName);
    this.view = DG.View.create();
    this.configureView();
  }

  protected view: DG.View;
  async addView(): Promise<void> {
    await this.initView();
    const name = this.parentAppName ? this.parentAppName + '/' + this.appName : this.appName;
    this.view.path = `/apps/${_package.name}/${name.replace(/\s/g, '')}/`;
    grok.shell.addView(this.view);
  }

  protected abstract getContent(): Promise<HTMLDivElement>;
  async initView(): Promise<void> {
    const content = await this.getContent();
    this.view.append(content);
  }

  protected configureView(): void {
    this.view.box = true;
    this.view.name = this.appName;

    const windows = grok.shell.windows;
    windows.showProperties = false;
    windows.showToolbox = false;
    windows.showHelp = false;
  }

  getView(): DG.View {
    return this.view;
  }
}

