// Tools for Diff Studio model errors

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {LINK, MISC, UI_TIME} from './ui-constants';

import '../css/app-styles.css';

/** Diff Studio model error */
export class ModelError extends Error {
  private helpUrl: string;
  private toHighlight: string = undefined;

  constructor(message: string, helpUrl: string, toHighlight?: string) {
    super(message);
    this.helpUrl = helpUrl;
    this.toHighlight = toHighlight;
  }

  public getHelpUrl() {
    return this.helpUrl;
  }

  public getToHighlight() {
    return this.toHighlight;
  }
}; // ModelError

/** Show hint indicating the model error */
export function showModelErrorHint(err: ModelError, tabControl: DG.TabControl) {
  const clearBtn = ui.button('Clear', () => popup.remove());
  const helpUrl = err.getHelpUrl();
  const toHighlight = err.getToHighlight();

  const msg = ui.divV([
    ui.h1('Model Error'),
    ui.divV([
      ui.markdown(err.message),
      ui.link('See help', () => {
        grok.shell.windows.help.visible = true;
        grok.shell.windows.help.showHelp(helpUrl);
      }),
      ui.divH([clearBtn], 'diff-studio-hint-btns-div'),
    ]),
  ]);

  const header = tabControl.header;
  header.hidden = false;
  const popup = ui.hints.addHint(header, msg, ui.hints.POSITION.RIGHT);
  header.hidden = true;
  tabControl.root.appendChild(popup);

  if (toHighlight !== undefined)
    highLightText(tabControl.root, toHighlight);
} // showModelErrorHint

/** Find the specified text in the root element and highlight it */
function highLightText(root: HTMLElement, text: string) {
  const lines = root.querySelectorAll('div.cm-line');
  const inds: number[] = [];

  lines.forEach((line, idx) => {
    if (line.textContent.includes(text))
      inds.push(idx);
  });

  if (inds.length === 1) {
    const line = lines[inds[0]];
    const numElements = root.querySelectorAll('div.cm-gutterElement');
    numElements[inds[0]].scrollIntoView();
    setTimeout(() => line.classList.add('diff-studio-highlight-text'), UI_TIME.WGT_CLICK);
  }
} // highLight

/** Return ModelError corresponding to ".. is not defined" */
export function getIsNotDefined(msg: string): ModelError {
  const idx = msg.indexOf(MISC.IS_NOT_DEF);
  const toHighlight = msg.slice(0, idx - 1);

  return new ModelError(
    `Error in formulas: **${toHighlight}** ${msg.slice(idx)}. Correct the model.`,
    LINK.MODEL_COMPONENTS,
    toHighlight,
  );
}

/** Return unexpected entity name */
export function getUnexpected(msg: string): ModelError {
  const start = msg.indexOf(`'`) + 1;
  const finish = msg.lastIndexOf(`'`);
  const toHighlight = msg.slice(start, finish);

  return new ModelError(
    `Error in formulas: undefined identifier **${toHighlight}**. Correct the model.`,
    LINK.MODEL_COMPONENTS,
    toHighlight,
  );
}

/** Return ModelError corresponding to "Cannot set properties of null" */
export function getNullOutput(): ModelError {
  return new ModelError(
    `Non-existent item in the **#output** block. Use entities specified in **#equations** or **#expressions**.`,
    LINK.INTERFACE,
    '#output',
  );
}
