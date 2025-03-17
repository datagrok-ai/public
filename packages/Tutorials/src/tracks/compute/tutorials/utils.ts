import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import '../../../../css/tutorial.css';

/** Max waiting time */
const MAX_TIME = 60000;

/** Timeout tick */
const TICK = 100;

/** Wait for current view has element */
export async function getViewWithElement(selector: string) {
  let view: DG.ViewBase | null = null;  
  let totalTime = 0;

  while (true) {
    view = grok.shell.v;

    if (view.root.querySelector(selector))
      break;

    if (totalTime > MAX_TIME)
      return null;

    await new Promise((resolve) => setTimeout(resolve, TICK));
    totalTime += TICK;
  }

  return view;  
}

/** Get the specified child */
export async function getElement(parent: HTMLElement | Document, selectors: string) {
  let element: HTMLElement | null;
  let totalTime = 0;

  while (true) {
    element = parent.querySelector(selectors);

    if (element)
      break;

    if (totalTime > MAX_TIME)
      return null;

    await new Promise((resolve) => setTimeout(resolve, TICK));
    totalTime += TICK;
  }

  return element;
}

/** Get the specified view */
export async function getView(name: string) {
  let view: DG.ViewBase;
  let totalTime = 0;

  while (true) {
    view = grok.shell.v;

    if (view.name.includes(name))
      return view;

    if (totalTime > MAX_TIME)
      return null;

    await new Promise((resolve) => setTimeout(resolve, TICK));
    totalTime += TICK;
  }
}

/** Describe viewers and return the Done button */
export function describeElements(roots: HTMLElement[], description: string[]): HTMLButtonElement {
  if (roots.length !== description.length)
    throw new Error('Non-equal size of viewer roots and descriptions');

  let idx = 0;
  let closeIcn: HTMLElement;
  let msg: HTMLDivElement;
  let popup: HTMLDivElement;

  const nextBtn = ui.button('next', () => {
    popup.remove();
    ++idx;
    step();
  }, 'Go to the next viewer');

  const prevBtn = ui.button('prev', () => {
    idx -= 1;    
    popup.remove();
    step();
  }, 'Go to the previous viewer');

  const doneBtn = ui.button('done', () => popup.remove(), 'Go to the next step');

  const btnsDiv = ui.divH([prevBtn, nextBtn, doneBtn]);
  btnsDiv.classList.add('tutorials-sci-comp-btns-div');

  const step = () => {
    if (idx < roots.length) {
      msg = ui.divV([ui.markdown(description[idx]), btnsDiv]);
      popup = ui.hints.addHint(roots[idx], msg, 'left');
      doneBtn.hidden = (idx < roots.length - 1);
      nextBtn.hidden = (idx === roots.length - 1);
      prevBtn.hidden = (idx < 1);
      
      closeIcn = popup.querySelector('i') as HTMLElement;
      closeIcn.onclick = () => doneBtn.click();
    }
  };

  step();

  return doneBtn;
}

/** Description of a single element */
export function singleDescription(root: HTMLElement, description: string, tooltip: string, position: ui.hints.POSITION = ui.hints.POSITION.LEFT) {
  const clearBtn = ui.button('ok', () => closeIcn.click(), tooltip);
  const btnDiv = ui.divH([clearBtn]);
  const msg = ui.divV([
    ui.markdown(description),
    btnDiv,
  ]);
  const popup = ui.hints.addHint(root, msg, position);
  const closeIcn = popup.querySelector('i') as HTMLElement;
  btnDiv.classList.add('tutorials-sci-comp-btns-div');

  return closeIcn;
}

/** Close windows */
export function closeWindows() {
  grok.shell.windows.showToolbox = false;
  grok.shell.windows.showHelp = false;
  grok.shell.windows.showContextPanel = false;
  grok.shell.windows.showConsole = false;
  grok.shell.windows.showVariables = false;
  grok.shell.windows.showTables = false;
  grok.shell.windows.showColumns = false;
}
