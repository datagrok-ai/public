// Tutorial features for UI elements spotlighting

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../../../css/ui-describer.css';

/**
 * Represents a rectangular area.
 */
type Rect = {
  left: number,
  top: number,
  width: number,
  height: number,
};

/**
 * Optional settings for computing a rectangle.
 */
interface IRectSettings {
  width?: number,
  height?: number,
  paddingBottom?: number,
};

/**
 * Spotlight target: can be an HTMLElement or a Rect.
 */
export type Spotlight = HTMLElement | Rect;

/**
 * Elements to highlight in the spotlight.
 */
export type SpotlightElements = {
  major: Spotlight,
  extra?: Spotlight,
};

/**
 * Represents one page/step in the tutorial.
 */
export type DescriptionPage = {
  root: HTMLElement,                  // Root element for hint positioning
  description: string | HTMLElement,  // Text or HTML content for the hint
  position: string,                   // Hint position relative to root
  nextBtnAction?: () => void,         // Optional callback for Next button
  prevBtnAction?: () => void,         // Optional callback for Previous button
  elements: SpotlightElements,        // Elements to spotlight
};

/**
 * Text labels for the tutorial buttons.
 */
export type ButtonsText = {
  next: string,
  prev: string,
  done: string,
};

/**
 * Default button labels.
 */
const DEFAULT_BTNS_TEXT: ButtonsText = {
  next: 'next',
  prev: 'prev',
  done: 'done',
};

/**
 * Represents the full tutorial.
 */
export type Tour = {
  pages: DescriptionPage[],    // Array of tutorial pages
  doneBtnAction?: () => void,  // Optional callback when tutorial ends
  btnsText?: ButtonsText,      // Optional custom button labels
};

/**
 * Returns the rectangle for an element with optional adjustments.
 * @param elem The HTMLElement to measure.
 * @param settings Optional settings to override width/height or add padding.
 * @returns The computed Rect.
 */
export function getRect(elem: HTMLElement, settings: IRectSettings = {}): Rect {
  const rect = elem.getBoundingClientRect();

  return {
    left: rect.left,
    top: rect.top,
    width: settings.width ?? rect.width,
    height: settings.height ?? (rect.height + (settings.paddingBottom ?? 0)),
  };
}

/**
 * Runs the tutorial describer with hints and spotlight.
 * @param tour The Tour object containing pages and optional callbacks.
 * @returns The "Done" button HTMLElement.
 */
export function runDescriber(tour: Tour): HTMLButtonElement {
  const pages = tour.pages;
  const pagesCount = pages.length;
  if (pagesCount < 1)
    throw new Error('Empty description pages');

  // Create the overlay for dimming background
  const overlay = createOverlay();

  let idx = 0;  // Current page index
  let closeIcn: HTMLElement;
  let msg: HTMLElement;
  let popup: HTMLDivElement;

  const btnsText = tour.btnsText ?? DEFAULT_BTNS_TEXT;

  // Create "Next" button
  const nextBtn = ui.button(btnsText.next, () => {
    popup.remove();
    clearSpotlight();

    const action = pages[idx].nextBtnAction;
    if (action != null)
      action();

    ++idx;
    step();
  });

  // Create "Prev" button
  const prevBtn = ui.button(btnsText.prev, () => {
    const action = pages[idx].prevBtnAction;
    if (action != null)
      action();

    idx -= 1;
    popup.remove();
    clearSpotlight();
    step();
  });

  // Create "Done" button
  const doneBtn = ui.button(btnsText.done, () => {
    popup.remove();
    clearSpotlight();
    overlay.remove();

    if (tour.doneBtnAction != null)
      tour.doneBtnAction();
  });

  // Container for buttons, aligned to the right
  const btnsDiv = ui.divH([prevBtn, nextBtn, doneBtn]);
  btnsDiv.style.marginLeft = 'auto';
  btnsDiv.style.marginRight = '0px';

  /**
   * Displays the current step/page.
   */
  const step = () => {
    if (idx < pagesCount) {
      const descr = pages[idx].description;

      // Create the message container (HTML element or Markdown)
      msg = ui.divV([descr instanceof HTMLElement ? descr : ui.markdown(descr), btnsDiv]);

      // Add hint popup
      popup = ui.hints.addHint(pages[idx].root, msg, pages[idx].position as ui.hints.POSITION);

      // Highlight elements with spotlight
      spotlight(popup, pages[idx].elements);

      // Show/hide buttons based on current step
      doneBtn.hidden = (idx < pagesCount - 1);
      nextBtn.hidden = (idx === pagesCount - 1);
      prevBtn.hidden = (idx < 1);

      // Setup close icon to finish tutorial
      closeIcn = popup.querySelector('i') as HTMLElement;
      closeIcn.onclick = () => doneBtn.click();
    }
  }; // step

  step();

  return doneBtn;
} // runDescriber

/**
 * Creates the dimming overlay with SVG mask for spotlight effect.
 * @returns The overlay HTMLElement.
 */
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

  // Black rectangles represent holes in overlay
  const popupHole = document.createElementNS(svgNS, 'rect');
  popupHole.id = 'tutorials-ui-describer-tour-popup-hole';
  popupHole.setAttribute('fill', 'black');

  const majorHole = document.createElementNS(svgNS, 'rect');
  majorHole.id = 'tutorials-ui-describer-tour-major-hole';
  majorHole.setAttribute('fill', 'black');

  const extraHole = document.createElementNS(svgNS, 'rect');
  extraHole.id = 'tutorials-ui-describer-tour-extra-hole';
  extraHole.setAttribute('fill', 'black');

  // Assemble mask
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

/**
 * Highlights the specified elements by adjusting the SVG mask holes.
 * @param popup The popup element to highlight.
 * @param elements The major and optional extra elements to spotlight.
 */
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
  
  const setHole = (hole: HTMLElement, elem: any) => {
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

  setHole(popupHole, popup);
  setHole(majorHole, elements.major);
  setHole(extraHole, elements.extra ?? popup);
} // spotlight

/**
 * Hides the spotlight overlay.
 */
function clearSpotlight() {
  const overlay = document.getElementById('tutorials-ui-describer-tour-overlay');
  if (overlay) overlay.style.display = 'none';
} // clearSpotlight
