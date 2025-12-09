// Tutorial features
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../../../css/ui-describer.css';
import { getIcn } from './utils';

type Rect = {
  left: number,
  top: number,
  width: number,
  height: number,
};

interface IRectSettings {
  width?: number,
  height?: number,
  paddingBottom?: number,
};

export type HighlightElement = HTMLElement | Rect;

export type DescriptionPage = {
  root: HTMLElement,
  description: string | HTMLElement,
  position: string,
  nextBtnAction?: () => void,
  prevBtnAction?: () => void,
  elementsToHighlight: HighlightElement[],
};

export type ButtonsText = {
  next: string,
  prev: string,
  done: string,
};

const DEFAULT_BTNS_TEXT: ButtonsText = {
  next: 'next',
  prev: 'prev',
  done: 'done',
};

export type Tour = {
  pages: DescriptionPage[],
  doneBtnAction?: () => void,
  btnsText?: ButtonsText,
};

function getHolesCount(pages: DescriptionPage[]): number {
  let res = 0;

  for (const page of pages) {

  }

  return res;
}

export function getRect(elem: HTMLElement, settings: IRectSettings = {}): Rect {
  const rect = elem.getBoundingClientRect();

  return {
    left: rect.left,
    top: rect.top,
    width: settings.width ?? rect.width,
    height: settings.height ?? (rect.height + (settings.paddingBottom ?? 0)),
  };
}

/** Run ui describer */
export function runDescriber(tour: Tour): HTMLButtonElement {
  const pages = tour.pages;
  const pagesCount = pages.length;
  if (pagesCount < 1)
    throw new Error('Empty description pages');

  const overlay = createOverlay();

  let idx = 0;
  let closeIcn: HTMLElement;
  let msg: HTMLElement;
  let popup: HTMLDivElement;

  const btnsText = tour.btnsText ?? DEFAULT_BTNS_TEXT;

  const nextBtn = ui.button(btnsText.next, () => {
    popup.remove();
    clearSpotlight();

    const action = pages[idx].nextBtnAction;
    if (action != null)
      action();

    ++idx;
    step();
  });

  const prevBtn = ui.button(btnsText.prev, () => {
    const action = pages[idx].prevBtnAction;
    if (action != null)
      action();

    idx -= 1;
    popup.remove();
    clearSpotlight();
    step();
  });

  const doneBtn = ui.button(btnsText.done, () => {
    popup.remove();
    clearSpotlight();
    overlay.remove();

    if (tour.doneBtnAction != null)
      tour.doneBtnAction();
  });

  const btnsDiv = ui.divH([prevBtn, nextBtn, doneBtn]);
  btnsDiv.style.marginLeft = 'auto';
  btnsDiv.style.marginRight = '0px';

  const step = () => {
    if (idx < pagesCount) {
      const descr = pages[idx].description;

      msg = ui.divV([descr instanceof HTMLElement ? descr : ui.markdown(descr), btnsDiv]);

      popup = ui.hints.addHint(pages[idx].root, msg, pages[idx].position as ui.hints.POSITION);

      const elementsToHighlight: HighlightElement[] = [popup];
      elementsToHighlight.push(...(pages[idx].elementsToHighlight ?? []));

      spotlight(elementsToHighlight);

      doneBtn.hidden = (idx < pagesCount - 1);
      nextBtn.hidden = (idx === pagesCount - 1);
      prevBtn.hidden = (idx < 1);

      closeIcn = popup.querySelector('i') as HTMLElement;
      closeIcn.onclick = () => doneBtn.click();
    }
  };

  step();

  return doneBtn;
} // runDescriber

/** Return div with circle and text */
export function getCircleWithText(color: string, text: string): HTMLElement {
  const circle = ui.span([]);
  circle.classList.add('bio-rx-sim-tutorial-circle');
  circle.style.backgroundColor = color;

  return ui.divH([circle, ui.divText(text)]);
} // getCircleWithText

/** Return legend */
export function getLegend(fileName: string, text: string): HTMLElement {
  const img = getIcn(fileName, 25, 20);
  img.style.marginRight = '4px';
  img.style.marginLeft = '20px';
  img.style.marginBottom = '5px';

  return ui.divH([img, ui.divText(text)]);
} // getLegend

function createOverlay(): HTMLElement {
  const existing = document.getElementById('tutorials-ui-describer-tour-overlay');
  if (existing)
    return existing;

  const overlay = document.createElement('div');
  overlay.id = 'tutorials-ui-describer-tour-overlay';
  overlay.style.display = 'none';


  const svgNS = 'http://www.w3.org/2000/svg';

  const svg = document.createElementNS(svgNS, 'svg');

  const defs = document.createElementNS(svgNS, 'defs');
  const mask = document.createElementNS(svgNS, 'mask');
  mask.setAttribute('id', 'mask');

  const full = document.createElementNS(svgNS, 'rect');
  full.setAttribute('x', '0');
  full.setAttribute('y', '0');
  full.setAttribute('width', '100%');
  full.setAttribute('height', '100%');
  full.setAttribute('fill', 'white');

  const hole1 = document.createElementNS(svgNS, 'rect');
  hole1.id = 'tutorials-ui-describer-tour-hole1';
  hole1.setAttribute('fill', 'black');

  const hole2 = document.createElementNS(svgNS, 'rect');
  hole2.id = 'tutorials-ui-describer-tour-hole2';
  hole2.setAttribute('fill', 'black');

  const hole3 = document.createElementNS(svgNS, 'rect');
  hole3.id = 'tutorials-ui-describer-tour-hole3';
  hole3.setAttribute('fill', 'black');

  const hole4 = document.createElementNS(svgNS, 'rect');
  hole4.id = 'tutorials-ui-describer-tour-hole4';
  hole4.setAttribute('fill', 'black');

  const hole5 = document.createElementNS(svgNS, 'rect');
  hole5.id = 'tutorials-ui-describer-tour-hole5';
  hole5.setAttribute('fill', 'black');

  mask.appendChild(full);
  mask.appendChild(hole1);
  mask.appendChild(hole2);
  mask.appendChild(hole3);
  mask.appendChild(hole4);
  mask.appendChild(hole5);
  defs.appendChild(mask);
  svg.appendChild(defs);

  const overlayRect = document.createElementNS(svgNS, 'rect');
  overlayRect.setAttribute('x', '0');
  overlayRect.setAttribute('y', '0');
  overlayRect.setAttribute('width', '100%');
  overlayRect.setAttribute('height', '100%');
  overlayRect.setAttribute('fill', 'rgba(0,0,0,0.6)');
  overlayRect.setAttribute('mask', 'url(#mask)');

  svg.appendChild(overlayRect);
  overlay.appendChild(svg);

  document.body.appendChild(overlay);
  return overlay;
} // createOverlay

function spotlight(elements: HighlightElement[]) {
  const count = elements.length;
  if (count < 1)
    return;

  const overlay = createOverlay();
  overlay.style.display = 'block';

  const holes = [
    document.getElementById('tutorials-ui-describer-tour-hole1')!,
    document.getElementById('tutorials-ui-describer-tour-hole2')!,
    document.getElementById('tutorials-ui-describer-tour-hole3')!,
    document.getElementById('tutorials-ui-describer-tour-hole4')!,
    document.getElementById('tutorials-ui-describer-tour-hole5')!,
  ];

  holes.forEach((hole, idx) => {
    const elem = elements[(idx < count) ? idx : 0];
    const rect = (elem instanceof HTMLElement) ? elem.getBoundingClientRect() : elem;
    const x = rect.left;
    const y = rect.top;
    const width = rect.width;
    const height = rect.height;

    hole.setAttribute('x', x.toString());
    hole.setAttribute('y', y.toString());
    hole.setAttribute('width', width.toString());
    hole.setAttribute('height', height.toString());
  });
} // spotlight

function clearSpotlight() {
  const overlay = document.getElementById('tutorials-ui-describer-tour-overlay');
  if (overlay) overlay.style.display = 'none';
} // clearSpotlight
