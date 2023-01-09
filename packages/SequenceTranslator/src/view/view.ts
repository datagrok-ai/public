/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MAIN, AXOLABS, SDF, DEFAULT_SEQUENCE} from './view-const';

/** Class responsible for the UI of the application */
export class View {
  constructor() {
    this.layout = grok.shell.newView('Sequence Translator', []);
    this.layout.box = true;
  }

  /** The master view containing app's main interface elements */
  private readonly layout: DG.View;
  private readonly tabControl: DG.TabControl;

  get view() { return this.layout; };
};
