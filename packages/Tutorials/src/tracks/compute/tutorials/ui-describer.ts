// Tutorial features
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../../../css/ui-describer.css';

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

export type Spotlight = HTMLElement | Rect;

export type SpotlightElements = {
  major: Spotlight,
  extra?: Spotlight,
};

export type DescriptionPage = {
  root: HTMLElement,
  description: string | HTMLElement,
  position: string,
  nextBtnAction?: () => void,
  prevBtnAction?: () => void,
  elements: SpotlightElements,
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

      spotlight(popup, pages[idx].elements);

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

  const popupHole = document.createElementNS(svgNS, 'rect');
  popupHole.id = 'tutorials-ui-describer-tour-popup-hole';
  popupHole.setAttribute('fill', 'black');

  const majorHole = document.createElementNS(svgNS, 'rect');
  majorHole.id = 'tutorials-ui-describer-tour-major-hole';
  majorHole.setAttribute('fill', 'black');

  const extraHole = document.createElementNS(svgNS, 'rect');
  extraHole.id = 'tutorials-ui-describer-tour-extra-hole';
  extraHole.setAttribute('fill', 'black');

  mask.appendChild(full);
  mask.appendChild(popupHole);
  mask.appendChild(majorHole);
  mask.appendChild(extraHole);
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

function spotlight(popup: HTMLElement, elements: SpotlightElements) {
  const overlay = createOverlay();
  overlay.style.display = 'block';

  const popupHole = document.getElementById('tutorials-ui-describer-tour-popup-hole');
  if (popupHole == null)
    throw new Error('Failed to spotlight elements: the popup hole not found');

  const majorHole = document.getElementById('tutorials-ui-describer-tour-major-hole');
  if (majorHole == null)
    throw new Error('Failed to spotlight elements: the major element hole not found');

  const extraHole = document.getElementById('tutorials-ui-describer-tour-extra-hole');
  if (extraHole == null)
    throw new Error('Failed to spotlight elements: the extra element hole not found');

  const spotlight = (hole: HTMLElement, elem: any) => {
    const rect = (elem instanceof HTMLElement) ? elem.getBoundingClientRect() : elem;
    const x = rect.left;
    const y = rect.top;
    const width = rect.width;
    const height = rect.height;

    hole.setAttribute('x', x.toString());
    hole.setAttribute('y', y.toString());
    hole.setAttribute('width', width.toString());
    hole.setAttribute('height', height.toString());    
  };

  spotlight(popupHole, popup);
  spotlight(majorHole, elements.major);
  spotlight(extraHole, elements.extra ?? popup);
} // spotlight

function clearSpotlight() {
  const overlay = document.getElementById('tutorials-ui-describer-tour-overlay');
  if (overlay) overlay.style.display = 'none';
} // clearSpotlight
