import {autoUpdate, computePosition, flip, offset, size} from '@floating-ui/dom';
import {Scope} from './scope.js';

/** Fired on the content element when the overlay closes for any reason (outside pointerdown,
 * Esc, anchor removal, scope disposal) — lets the owning control resync its own open state. */
export const OVERLAY_CLOSE_EVENT = 'u2-overlay-close';

/** The single body-level layer host. One container for every u2 popup, so z-index and
 * dismissal are decided in one place instead of per control. */
export class Overlay {
  private static _host: HTMLElement | undefined;
  private static _layer = 0;

  /** Monotonic layer counter shared by everything in the host (popups, dialogs, their
   * backdrops): whatever opened last paints on top, so a combobox inside a dialog is
   * never occluded by it. Scoped to the u2 overlay host — the platform's own surfaces
   * keep their banded z-orders. */
  static nextLayer(): number {
    return ++Overlay._layer;
  }

  static get host(): HTMLElement {
    if (!Overlay._host || !Overlay._host.isConnected) {
      const el = document.createElement('div');
      el.className = 'u2-overlay';
      el.style.cssText = 'position: fixed; left: 0; top: 0; width: 0; height: 0; z-index: var(--dg-z-overlay);';
      document.body.append(el);
      Overlay._host = el;
    }
    return Overlay._host;
  }

  /** Shows `content` under `anchor` in the shared layer. The returned close is idempotent and
   * also registered on `scope`, so disposing the owner tears the overlay down; closing drops that
   * registration again, so an input reopened all day never piles dead closures up on its scope. */
  static show(anchor: HTMLElement, content: HTMLElement, scope: Scope): () => void {
    // Fixed strategy: viewport coordinates, independent of offsetParent chains — the app's
    // dock/view containers are positioned ancestors that break absolute-strategy math.
    content.style.position = 'fixed';
    content.style.left = '0';
    content.style.top = '0';
    content.style.zIndex = String(Overlay.nextLayer());
    Overlay.host.append(content);

    let closed = false;
    let stopAutoUpdate: () => void = () => {};

    const onPointerDown = (e: Event) => {
      const target = e.target as Node;
      if (!content.contains(target) && !anchor.contains(target))
        close();
    };
    const onKeyDown = (e: KeyboardEvent) => {
      if (e.key === 'Escape')
        close();
    };
    // autoUpdate reacts to scroll/resize but not to detach; dock_spawn reparents constantly.
    const detachObserver = new MutationObserver(() => {
      if (!anchor.isConnected)
        close();
    });

    const close = () => {
      if (closed)
        return;
      closed = true;
      stopAutoUpdate();
      detachObserver.disconnect();
      document.removeEventListener('pointerdown', onPointerDown, true);
      document.removeEventListener('keydown', onKeyDown);
      content.remove();
      scope.disown(close);
      content.dispatchEvent(new CustomEvent(OVERLAY_CLOSE_EVENT));
    };

    const position = () => {
      if (!anchor.isConnected) {
        close();
        return;
      }
      computePosition(anchor, content, {
        strategy: 'fixed',
        placement: 'bottom-start',
        middleware: [
          offset(1),
          flip({padding: 4}),
          size({
            padding: 4,
            apply: ({rects, availableHeight, elements}) => {
              elements.floating.style.minWidth = `${rects.reference.width}px`;
              elements.floating.style.maxHeight = `${Math.max(availableHeight, 64)}px`;
            },
          }),
        ],
      }).then(({x, y}) => {
        content.style.left = `${x}px`;
        content.style.top = `${y}px`;
      });
    };

    detachObserver.observe(document.body, {childList: true, subtree: true});
    document.addEventListener('pointerdown', onPointerDown, true);
    document.addEventListener('keydown', onKeyDown);
    // animationFrame: dock_spawn moves/reparents anchors outside scroll/resize (pane
    // activation, docking) — the default observers miss pure element moves.
    stopAutoUpdate = autoUpdate(anchor, content, position, {animationFrame: true});
    scope.own(close);
    return close;
  }
}
