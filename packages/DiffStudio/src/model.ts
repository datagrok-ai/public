// Items for Model catalog

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../css/app-styles.css';
import {DiffStudio, UiOptions} from './app';

/** Diff Studio model  */
export class Model {
  private model: string;
  private uiOptions: UiOptions;
  private info: string;

  constructor(model: string, uiOptions: UiOptions, info: string) {
    this.model = model;
    this.uiOptions = uiOptions;
    this.info = info;
  }

  private showHelpPanel(): void {
    grok.shell.windows.help.visible = true;
    const helpMD = ui.markdown(this.info);
    helpMD.classList.add('diff-studio-demo-app-div-md');
    const divHelp = ui.div([helpMD], 'diff-studio-demo-app-div-help');
    grok.shell.windows.help.showHelp(divHelp);
    grok.shell.windows.context.visible = true;
    grok.shell.windows.showContextPanel = false;
    grok.shell.windows.showProperties = false;
    grok.shell.windows.help.visible = true;
  }

  public async run(): Promise<void> {
    const solver = new DiffStudio(true, undefined, undefined, undefined, this.uiOptions);
    await solver.runModel(this.model);
  }

  public async runDemo(): Promise<void> {
    const solver = new DiffStudio(true, undefined, undefined, undefined, this.uiOptions);
    await solver.runModel(this.model);
    this.showHelpPanel();
  }
}; // Model
