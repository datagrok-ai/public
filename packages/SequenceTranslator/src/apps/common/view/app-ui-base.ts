/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

export abstract class AppUIBase {
  constructor(protected appName: string, protected parentAppName?: string) { }
  protected abstract constructView(): Promise<DG.ViewBase>;

  async getAppView(): Promise<DG.ViewBase> {
    const progressIndicator: DG.TaskBarProgressIndicator =
      DG.TaskBarProgressIndicator.create(`Loading ${this.appName}...`);

    const currentView = grok.shell.v?.root;
    if (currentView)
      ui.setUpdateIndicator(currentView, true);

    try {
      const view = await this.constructView();
      return view;
    } finally {
      progressIndicator.close();
      if (currentView)
        ui.setUpdateIndicator(currentView, false);
    }
  }
}

