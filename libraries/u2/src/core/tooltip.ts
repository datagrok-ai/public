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
  private static _timer = 0;
  private static _hiddenAt = 0;

  static get root(): HTMLElement {
    if (!Tooltip._root || !Tooltip._root.isConnected) {
      const el = document.createElement('div');
      el.className = 'u2-tooltip';
      el.dataset.u2 = 'tooltip';
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
   * Both listeners are released with `scope`. */
  static bind(el: HTMLElement, content: string | (() => string | HTMLElement), scope: Scope): void {
    const enter = () => {
      const show = () => Tooltip.show(typeof content === 'function' ? content() : content, el);
      // Warm state: leaving one bound element and entering another switches with no delay.
      if (Tooltip.isVisible || Date.now() - Tooltip._hiddenAt < WARM_MS)
        show();
      else {
        window.clearTimeout(Tooltip._timer);
        Tooltip._timer = window.setTimeout(show, SHOW_DELAY_MS);
      }
    };
    const leave = () => Tooltip.hide();
    el.addEventListener('pointerenter', enter);
    el.addEventListener('pointerleave', leave);
    el.addEventListener('pointerdown', leave);
    scope.own(() => {
      el.removeEventListener('pointerenter', enter);
      el.removeEventListener('pointerleave', leave);
      el.removeEventListener('pointerdown', leave);
      if (Tooltip._anchor === el)
        Tooltip.hide();
    });
  }

  static show(content: string | HTMLElement, anchor: HTMLElement): void {
    window.clearTimeout(Tooltip._timer);
    const root = Tooltip.root;
    root.textContent = '';
    if (typeof content === 'string')
      root.textContent = content;
    else
      root.append(content);
    root.style.display = 'block';
    Tooltip._anchor = anchor;
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
    if (!Tooltip.isVisible)
      return;
    Tooltip._anchor = undefined;
    Tooltip._hiddenAt = Date.now();
    document.removeEventListener('pointerdown', Tooltip._onPointerDown, true);
    document.removeEventListener('keydown', Tooltip._onKeyDown);
    document.removeEventListener('scroll', Tooltip._onScroll, true);
    const root = Tooltip.root;
    root.style.display = 'none';
    root.textContent = '';
  }

  private static _onPointerDown = () => Tooltip.hide();

  private static _onScroll = () => Tooltip.hide();

  private static _onKeyDown = (e: KeyboardEvent) => {
    if (e.key === 'Escape')
      Tooltip.hide();
  };
}
