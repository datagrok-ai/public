/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {tryCatch} from '../model/helpers';

export abstract class AppUIBase {
  constructor(protected appName: string, protected parentAppName?: string) { }
  abstract addView(): Promise<void>;

  async initializeAppLayout(): Promise<void> {
    const progressIndicator: DG.TaskBarProgressIndicator =
      DG.TaskBarProgressIndicator.create(`Loading ${this.appName}...`);

    const currentView = grok.shell.v?.root;
    if (currentView)
      ui.setUpdateIndicator(currentView, true);

    await tryCatch(async () => {
      await this.addView();
    }, () => progressIndicator.close());

    if (currentView)
      ui.setUpdateIndicator(currentView, false);
  }
}

