import {computePosition, flip, offset, shift} from '@floating-ui/dom';
import {Overlay} from './overlay.js';
import {Scope} from './scope.js';

const SHOW_DELAY_MS = 300;
const WARM_MS = 100;

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

  /** Shows the tooltip for `el` after a hover delay; `content` is evaluated at show time.
   * The listeners are released with `scope` — by default whoever is building (an `each()` row). */
  static bind(el: HTMLElement, content: string | (() => string | HTMLElement), scope?: Scope): void {
    const owner = scope ?? Scope.ambient;
    if (owner === undefined)
      throw new Error('u2: Tooltip.bind needs an owner — pass a scope, or bind while a control builds');
    const enter = () => {
      const show = () => Tooltip.show(typeof content === 'function' ? content() : content, el);
      // Warm state: leaving one bound element and entering another switches with no delay.
      if (Tooltip.isVisible || Date.now() - Tooltip._hiddenAt < WARM_MS)
        show();
      else {
        window.clearTimeout(Tooltip._timer);
        Tooltip._pending = el;
        Tooltip._timer = window.setTimeout(show, SHOW_DELAY_MS);
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
    owner.own(() => {
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
