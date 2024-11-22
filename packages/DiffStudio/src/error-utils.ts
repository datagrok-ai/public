// Tools for Diff Studio model errors

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../css/app-styles.css';

export class ModelError extends Error {
  public helpUrl: string;

  constructor(message: string, helpUrl: string) {
    super(message);
    this.helpUrl = helpUrl;
  }
};

export function showModelErrorHint(err: ModelError, root: HTMLElement) {
  const clearBtn = ui.button('Clear', () => popup.remove());

  const msg = ui.divV([
    ui.h1('Model Error'),
    ui.divV([
      ui.markdown(err.message),
      ui.link('See help', () => {
        grok.shell.windows.help.visible = true;
        grok.shell.windows.help.showHelp(err.helpUrl);
      }),
      ui.divH([clearBtn], 'diff-studio-hint-btns-div'),
    ]),
  ]);

  root.hidden = false;
  const popup = ui.hints.addHint(root, msg, ui.hints.POSITION.RIGHT);
  root.hidden = true;
}
