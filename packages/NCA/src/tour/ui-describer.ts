/**
 * Spotlight-tour primitive — vendored from
 * `C:/Users/vmaka/Datagrok/public/packages/Tutorials/src/tracks/compute/tutorials/ui-describer.ts`
 * (Datagrok Tutorials package). Adaptations for nca-studio:
 *  - DOM IDs renamed `tutorials-ui-describer-*` → `nca-tour-*` so the overlay
 *    does not collide when the Tutorials package is also loaded in the same
 *    session.
 *  - CSS path adjusted to `../../css/tour.css`.
 *  - `import type` for type-only imports (`verbatimModuleSyntax: true`).
 *  - No behavioural changes — `runDescriber` is the unchanged spotlight tour
 *    engine: dim background via SVG mask, cut holes around the popup +
 *    spotlighted elements, navigate pages via prev / next / done.
 *
 * Public surface:
 *   - {@link Spotlight}, {@link SpotlightElements}, {@link DescriptionPage},
 *     {@link Tour}, {@link ButtonsText}
 *   - {@link runDescriber} — start a tour
 */

import * as ui from 'datagrok-api/ui';

// @ts-ignore  — CSS module not declared in tsconfig; webpack resolves at build time.
import '../../css/tour.css';

/** A rectangular area in viewport coordinates. */
type Rect = {
  left: number,
  top: number,
  width: number,
  height: number,
};

/** Optional measurement adjustments for {@link getRect}. */
interface IRectSettings {
  width?: number,
  height?: number,
  paddingBottom?: number,
}

/** Spotlight target: either a live HTMLElement or a fixed Rect. */
export type Spotlight = HTMLElement | Rect;

/** Elements to highlight in a spotlight step. `extra` is optional (defaults
 *  to the popup itself, so only the major target stands out). */
export type SpotlightElements = {
  major: Spotlight,
  extra?: Spotlight,
};

/** One step / page of a tour. */
export type DescriptionPage = {
  /** Anchor element the hint popup is attached to. */
  root: HTMLElement,
  /** Body content — Markdown string or pre-built HTMLElement. */
  description: string | HTMLElement,
  /** Hint position relative to `root` (passed through to `ui.hints.POSITION`). */
  position: string,
  /** Optional callback when the user clicks Next on this page. */
  nextBtnAction?: () => void,
  /** Optional callback when the user clicks Prev on this page. */
  prevBtnAction?: () => void,
  /** Optional callback fired just BEFORE the hint is rendered for this page.
   *  Use to bring `root` into view (e.g. `scrollIntoView`) so the spotlight
   *  geometry is computed against an on-screen rect. */
  beforeShow?: () => void,
  /** Elements to spotlight (major + optional extra). */
  elements: SpotlightElements,
};

/** Text labels for the tour buttons. */
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

/** Full tour definition. */
export type Tour = {
  pages: DescriptionPage[],
  doneBtnAction?: () => void,
  btnsText?: ButtonsText,
};

/** Return an element's bounding rect with optional size / padding overrides.
 *  Module-private — re-export from this file when a tour step needs to
 *  spotlight a custom sub-region (e.g. a column inside a grid). */
// eslint-disable-next-line @typescript-eslint/no-unused-vars
function getRect(elem: HTMLElement, settings: IRectSettings = {}): Rect {
  const rect = elem.getBoundingClientRect();
  return {
    left: rect.left,
    top: rect.top,
    width: settings.width ?? rect.width,
    height: settings.height ?? (rect.height + (settings.paddingBottom ?? 0)),
  };
}

/** Start a spotlight tour. Returns the Done button so callers can drive
 *  programmatic close (e.g. on view-close cleanup). */
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

  const step = (): void => {
    if (idx < pagesCount) {
      const beforeShow = pages[idx].beforeShow;
      if (beforeShow != null) beforeShow();
      const descr = pages[idx].description;
      msg = ui.divV([descr instanceof HTMLElement ? descr : ui.markdown(descr), btnsDiv]);
      popup = ui.hints.addHint(pages[idx].root, msg, pages[idx].position as ui.hints.POSITION);
      spotlight(popup, pages[idx].elements);
      doneBtn.hidden = (idx < pagesCount - 1);
      nextBtn.hidden = (idx === pagesCount - 1);
      prevBtn.hidden = (idx < 1);
      closeIcn = popup.querySelector('i') as HTMLElement;
      closeIcn.onclick = (): void => doneBtn.click();
    }
  };

  step();
  return doneBtn;
}

const OVERLAY_ID = 'nca-tour-overlay';
const POPUP_HOLE_ID = 'nca-tour-popup-hole';
const MAJOR_HOLE_ID = 'nca-tour-major-hole';
const EXTRA_HOLE_ID = 'nca-tour-extra-hole';

/** Build or retrieve the dimming SVG overlay used to spotlight elements. */
function createOverlay(): HTMLElement {
  const existing = document.getElementById(OVERLAY_ID);
  if (existing)
    return existing;

  const overlay = document.createElement('div');
  overlay.id = OVERLAY_ID;
  overlay.style.display = 'none';

  const svgNS = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(svgNS, 'svg');

  const defs = document.createElementNS(svgNS, 'defs');
  const mask = document.createElementNS(svgNS, 'mask');
  mask.setAttribute('id', 'nca-tour-mask');

  const full = document.createElementNS(svgNS, 'rect');
  full.setAttribute('x', '0');
  full.setAttribute('y', '0');
  full.setAttribute('width', '100%');
  full.setAttribute('height', '100%');
  full.setAttribute('fill', 'white');

  const popupHole = document.createElementNS(svgNS, 'rect');
  popupHole.id = POPUP_HOLE_ID;
  popupHole.setAttribute('fill', 'black');

  const majorHole = document.createElementNS(svgNS, 'rect');
  majorHole.id = MAJOR_HOLE_ID;
  majorHole.setAttribute('fill', 'black');

  const extraHole = document.createElementNS(svgNS, 'rect');
  extraHole.id = EXTRA_HOLE_ID;
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
  overlayRect.setAttribute('mask', 'url(#nca-tour-mask)');

  svg.appendChild(overlayRect);
  overlay.appendChild(svg);

  document.body.appendChild(overlay);
  return overlay;
}

/** Position the SVG-mask holes over the popup + spotlight targets. */
function spotlight(popup: HTMLElement, elements: SpotlightElements): void {
  const overlay = createOverlay();
  overlay.style.display = 'block';

  const popupHole = document.getElementById(POPUP_HOLE_ID);
  if (popupHole == null)
    throw new Error('Failed to spotlight elements: the popup hole not found');

  const majorHole = document.getElementById(MAJOR_HOLE_ID);
  if (majorHole == null)
    throw new Error('Failed to spotlight elements: the major element hole not found');

  const extraHole = document.getElementById(EXTRA_HOLE_ID);
  if (extraHole == null)
    throw new Error('Failed to spotlight elements: the extra element hole not found');

  const setHole = (hole: HTMLElement, elem: Spotlight): void => {
    const rect = (elem instanceof HTMLElement) ? elem.getBoundingClientRect() : elem;
    hole.setAttribute('x', rect.left.toString());
    hole.setAttribute('y', rect.top.toString());
    hole.setAttribute('width', rect.width.toString());
    hole.setAttribute('height', rect.height.toString());
  };

  setHole(popupHole, popup);
  setHole(majorHole, elements.major);
  setHole(extraHole, elements.extra ?? popup);
}

/** Hide (but do not remove) the spotlight overlay between steps. */
function clearSpotlight(): void {
  const overlay = document.getElementById(OVERLAY_ID);
  if (overlay) overlay.style.display = 'none';
}
