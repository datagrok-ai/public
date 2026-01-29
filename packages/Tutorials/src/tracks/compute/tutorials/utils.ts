import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import '../../../../css/tutorial.css';
import '../../../../css/ui-describer.css';
import {_package} from '../../../package';

/** Max waiting time */
const MAX_TIME = 60000;

/** Timeout tick */
const TICK = 100;

/** Starting pause */
export const PAUSE = 500;

/** Editor delay */
export const DELAY = 2000;

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
  try {
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
  } catch (err) {
    return null;
  }
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

/** Return icon from the file */
export function getIcn(fileName: string, width: number, height: number, tooltip?: string): HTMLDivElement {
  const path = `${_package.webRoot}images/${fileName}`;
  const img = ui.image(path, width, height);
  const div = ui.div(img);

  if (tooltip != null)
    ui.tooltip.bind(div, tooltip);

  return div;
}

export function getTextWithSlider(textBefore: string, textAfter: string): HTMLElement {
  const circle = ui.span([]);
  circle.classList.add('tutorials-ui-describer-diff-studio-slider-in-text');
  circle.style.backgroundColor = 'var(--blue-1)';
  const sliderWithText = ui.divH([
    ui.divText(textBefore),
    circle,
    ui.divText(textAfter),
  ]);
  sliderWithText.classList.add('tutorials-ui-describer-diff-studio-text-with-slider');

  return sliderWithText;
} // getTextWithSlider

export function simulateMouseEventsWithMove(el: HTMLElement, moveX: number, moveY: number) {
  const rect = el.getBoundingClientRect();
  const startX = rect.left + rect.width / 2;
  const startY = rect.top + rect.height / 2;

  const downEvent = new MouseEvent("mousedown", {
    bubbles: true,
    clientX: startX,
    clientY: startY,
  });
  el.dispatchEvent(downEvent);

  const moveEvent = new MouseEvent("mousemove", {
    bubbles: true,
    clientX: startX + moveX,
    clientY: startY + moveY,
  });

  setTimeout(() => {
    el.dispatchEvent(moveEvent);
  }, 50);

  const upEvent = new MouseEvent("mouseup", {
    bubbles: true,
    clientX: startX + 1,
    clientY: startY + 1,
  });

  setTimeout(() => {
    el.dispatchEvent(upEvent);
  }, 100);
}

export function getLegendDiv(title: string, lines: string[]): HTMLElement {    
  return ui.divV([
    ui.markdown(title),
    ...lines.map((line) => {
      const md = ui.markdown(line);
      md.classList.add('tutorials-ui-describer-diff-studio-legend-line');

      return md;
    }),
  ]);
}

export function getBallFlightModelLegend(): HTMLElement {
  return getLegendDiv('# Simulation 🏀\n\nThis model takes the ball and thrown parameters, and computes:', [
    '* the ball flight trajectory',
    '* max height and distance'      
  ]);
}  

export function buildToggleOverlay(toggle: HTMLElement): HTMLElement {
  const overlay = document.createElement('div');
  overlay.style.position = 'absolute';
  overlay.style.left = '1px';
  overlay.style.width = '28px';
  overlay.style.height = '28px';
  overlay.style.background = 'transparent';
  overlay.style.pointerEvents = 'none';
  toggle.appendChild(overlay);

  return overlay;
}
