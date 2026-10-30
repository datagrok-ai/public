import {computePosition, flip, offset, shift} from '@floating-ui/dom';
import {Overlay} from './overlay.js';
import {Ownership} from './owner.js';
import type {Owner} from './owner.js';

const WARM_MS = 100;
const TOOLTIP = Symbol.for('u2.tooltip');

/** The single hover tooltip. One reused element in the shared {@link Overlay} layer — it is
 * never disposed, so a bound element only ever costs two listeners on its owner's scope.
 * Positioned directly (not via `Overlay.show`, whose outside-pointerdown dismissal and close
 * event belong to popups, not to hover). */
export class Tooltip {
  private static _root: HTMLElement | undefined;
  private static _anchor: HTMLElement | undefined;
  private static _pending: HTMLElement | undefined;
  private static _timer = 0;
  private static _hiddenAt = 0;

  /** How long a hover waits before the tooltip shows, in ms; 0 shows it at once (a test). */
  static delay = 300;

  static get root(): HTMLElement {
    if (!Tooltip._root || !Tooltip._root.isConnected) {
      const el = document.createElement('div');
      el.className = 'u2-tooltip';
      el.dataset.u2 = 'tooltip';
      el.id = 'u2-tooltip';
      el.setAttribute('role', 'tooltip');
      el.style.cssText = 'position: fixed; left: 0; top: 0; display: none; pointer-events: none;';
      Overlay.host.append(el);
      Tooltip._root = el;
    }
    return Tooltip._root;
  }

  static get isVisible(): boolean {
    return Tooltip._anchor !== undefined;
  }

  /** The elements whose `title` {@link titleIfCut} wrote. */
  private static readonly _cutTitles = new WeakSet<HTMLElement>();

  /** Full text on hover for a line that is cut: `el`'s whole text as its `title` while it runs past its box, gone
   * once it fits. A title the element carries of its own is left alone. */
  static titleIfCut(el: HTMLElement): void {
    const ours = Tooltip._cutTitles.has(el);
    if (el.scrollWidth > el.clientWidth) {
      if (!ours && el.title !== '')
        return;
      el.title = (el.textContent ?? '').trim();
      Tooltip._cutTitles.add(el);
    }
    else if (ours) {
      el.removeAttribute('title');
      Tooltip._cutTitles.delete(el);
    }
  }

  /** Shows the tooltip for `el` after a hover delay; `content` is evaluated at show time.
   * The listeners are released with `owner` — by default whoever is building (an `each()` row). */
  static bind(el: HTMLElement, content: string | (() => string | HTMLElement), options?: {owner?: Owner}): void {
    const owner = Ownership.scopeOf(options?.owner, 'Tooltip.bind()');
    const enter = () => {
      const show = () => Tooltip.show(typeof content === 'function' ? content() : content, el);
      // Warm state: leaving one bound element and entering another switches with no delay.
      if (Tooltip.isVisible || Date.now() - Tooltip._hiddenAt < WARM_MS || Tooltip.delay <= 0)
        show();
      else {
        window.clearTimeout(Tooltip._timer);
        Tooltip._pending = el;
        Tooltip._timer = window.setTimeout(show, Tooltip.delay);
      }
    };
    const leave = () => Tooltip.hide();
    // keyboard focus shows it too, so a keyboard user reaches the same text; a click's focus does not
    const focus = (e: Event) => {
      if (Tooltip._keyboardFocus(e.target as Element))
        enter();
    };
    el.addEventListener('pointerenter', enter);
    el.addEventListener('pointerleave', leave);
    el.addEventListener('pointerdown', leave);
    el.addEventListener('focusin', focus);
    el.addEventListener('focusout', leave);
    // what the engine's operability probe compares an icon button's accessible name with
    const slot = el as unknown as Record<symbol, unknown>;
    slot[TOOLTIP] = content;
    owner.own(() => {
      if (slot[TOOLTIP] === content)
        delete slot[TOOLTIP];
      el.removeEventListener('pointerenter', enter);
      el.removeEventListener('pointerleave', leave);
      el.removeEventListener('pointerdown', leave);
      el.removeEventListener('focusin', focus);
      el.removeEventListener('focusout', leave);
      // the pending show too: an element torn down under the pointer never gets its pointerleave,
      // so a timer left running would show the tooltip against an anchor that is no longer there
      if (Tooltip._anchor === el || Tooltip._pending === el)
        Tooltip.hide();
    });
  }

  static show(content: string | HTMLElement, anchor: HTMLElement): void {
    window.clearTimeout(Tooltip._timer);
    Tooltip._pending = undefined;
    if (!anchor.isConnected || content === '')
      return;
    const root = Tooltip.root;
    root.textContent = '';
    if (typeof content === 'string')
      root.textContent = content;
    else
      root.append(content);
    root.style.display = 'block';
    if (Tooltip._anchor !== undefined && Tooltip._anchor !== anchor)
      Tooltip._describe(Tooltip._anchor, false);
    Tooltip._anchor = anchor;
    Tooltip._describe(anchor, true);
    document.addEventListener('pointerdown', Tooltip._onPointerDown, true);
    document.addEventListener('keydown', Tooltip._onKeyDown);
    document.addEventListener('scroll', Tooltip._onScroll, true);
    computePosition(anchor, root, {
      strategy: 'fixed',
      placement: 'bottom-start',
      middleware: [offset(12), flip({padding: 4}), shift({padding: 4})],
    }).then(({x, y}) => {
      if (Tooltip._anchor !== anchor)
        return;
      root.style.left = `${x}px`;
      root.style.top = `${y}px`;
    });
  }

  static hide(): void {
    window.clearTimeout(Tooltip._timer);
    Tooltip._pending = undefined;
    if (!Tooltip.isVisible)
      return;
    Tooltip._describe(Tooltip._anchor!, false);
    Tooltip._anchor = undefined;
    Tooltip._hiddenAt = Date.now();
    document.removeEventListener('pointerdown', Tooltip._onPointerDown, true);
    document.removeEventListener('keydown', Tooltip._onKeyDown);
    document.removeEventListener('scroll', Tooltip._onScroll, true);
    const root = Tooltip.root;
    root.style.display = 'none';
    root.textContent = '';
  }

  /** The anchor names the open tooltip in its `aria-describedby`, beside whatever it already lists. */
  private static _describe(anchor: HTMLElement, on: boolean): void {
    const ids = (anchor.getAttribute('aria-describedby') ?? '').split(/\s+/)
      .filter((id) => id !== '' && id !== 'u2-tooltip');
    if (on)
      ids.push('u2-tooltip');
    if (ids.length > 0)
      anchor.setAttribute('aria-describedby', ids.join(' '));
    else
      anchor.removeAttribute('aria-describedby');
  }

  private static _keyboardFocus(el: Element | null): boolean {
    try {
      return el?.matches(':focus-visible') ?? false;
    }
    catch (_) {
      return false;
    }
  }

  private static _onPointerDown = () => Tooltip.hide();

  private static _onScroll = () => Tooltip.hide();

  private static _onKeyDown = (e: KeyboardEvent) => {
    if (e.key === 'Escape')
      Tooltip.hide();
  };
}
