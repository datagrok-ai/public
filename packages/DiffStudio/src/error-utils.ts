// Tools for Diff Studio model errors

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {LINK, MISC} from './ui-constants';

import '../css/app-styles.css';

export class ModelError extends Error {
  public helpUrl: string;
  public toHighlight: string = undefined;

  constructor(message: string, helpUrl: string, toHighlight?: string) {
    super(message);
    this.helpUrl = helpUrl;
    this.toHighlight = toHighlight;
  }
};

export function showModelErrorHint(err: ModelError, tabControl: DG.TabControl) {
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

  const header = tabControl.header;
  header.hidden = false;
  const popup = ui.hints.addHint(header, msg, ui.hints.POSITION.RIGHT);
  header.hidden = true;

  if (err.toHighlight !== undefined)
    highLight(tabControl.root, err.toHighlight);
}

function highLight(root: HTMLElement, text: string) {
  const lines = root.querySelectorAll('div.cm-line');
  const inds: number[] = [];

  lines.forEach((line, idx) => {
    if (line.textContent.includes(text))
      inds.push(idx);
  });

  if (inds.length === 1) {
    const line = lines[inds[0]];
    line.classList.add('diff-studio-highlight-text');
    line.scrollIntoView();
  }
}

export function getIsNotDefined(msg: string): ModelError {
  const idx = msg.indexOf(MISC.IS_NOT_DEF);
  const toHighlight = msg.slice(0, idx - 1);

  return new ModelError(
    `Error in formulas: **${toHighlight}** ${msg.slice(idx)}. Correct the model.`,
    LINK.BASIC_MODEL,
    toHighlight,
  );
}

export function getUnexpected(msg: string): ModelError {
  const start = msg.indexOf(`'`) + 1;
  const finish = msg.lastIndexOf(`'`);
  const toHighlight = msg.slice(start, finish);

  return new ModelError(
    `Error in formulas: undefined identifier **${toHighlight}**. Correct the model.`,
    LINK.BASIC_MODEL,
    toHighlight,
  );
}
