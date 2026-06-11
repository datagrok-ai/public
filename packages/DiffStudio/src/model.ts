// Items for Model catalog

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// @ts-ignore
import '../css/app-styles.css';
import {DiffStudio, UiOptions} from './app';

export type ModelInfo = {
  equations: string;
  info: string;
  uiOptions: UiOptions;
};

/** Diff Studio model  */
export class Model {
  private model: string;
  private uiOptions: UiOptions;
  private info: string;

  constructor(modelInfo: ModelInfo) {
    this.model = modelInfo.equations;
    this.uiOptions = modelInfo.uiOptions;
    this.info = modelInfo.info;
  }

  private async showHelpPanel(): Promise<void> {
    const helpMD = ui.markdown(this.info);
    helpMD.classList.add('diff-studio-demo-app-div-md');
    const divHelp = ui.div([helpMD], 'diff-studio-demo-app-div-help');

    await new Promise<void>((r) => setTimeout(r, 1000));
    grok.shell.windows.help.showHelp(divHelp);
    grok.shell.windows.context.visible = true;
    grok.shell.windows.showContextPanel = false;
    grok.shell.windows.showProperties = false;
    grok.shell.windows.help.visible = true;
  }

  public run(): Promise<void> {
    const solver = new DiffStudio(true, undefined, undefined, undefined, this.uiOptions);
    return solver.runModel(this.model);
  }

  public async runDemo(): Promise<void> {
    const solver = new DiffStudio(true, undefined, undefined, undefined, this.uiOptions);
    await solver.runModel(this.model);
    await this.showHelpPanel();
  }
}; // Model
