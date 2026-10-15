import {autoUpdate, computePosition, flip, offset, size} from '@floating-ui/dom';
import {Scope, debug} from './scope.js';
import {Control} from './component.js';

/** Fired on the content element when the overlay closes for any reason (outside pointerdown,
 * Esc, anchor removal, scope disposal) — lets the owning control resync its own open state. */
export const OVERLAY_CLOSE_EVENT = 'u2-overlay-close';

export interface OverlayOptions {
  /** A ceiling on the height the overlay is given — a CSS length or pixels. The room under the
   * anchor is written as an inline `max-height`, so a popup that wants to stay shorter than that
   * says so here rather than outranking the inline value from a stylesheet. */
  maxHeight?: number | string;
}

/** A dialog as the overlay stack sees it. */
export interface OverlayDialog {
  readonly root: HTMLElement;
  /** A modal dialog on top also takes an Escape whose focus fell out onto the page. */
  readonly modal: boolean;
  /** The Escape meant for this dialog: cancels it, or aborts its pending OK. */
  cancelByKey(): void;
  close(): void;
}

interface Entry {
  content: HTMLElement;
  close: () => void;
  anchor?: HTMLElement;
  /** The owner's element keyboard focus lands on; it carries `aria-controls`/`aria-expanded`. */
  focus?: HTMLElement;
  dialog?: OverlayDialog;
}

const FOCUSABLE = 'input, textarea, select, button, [role="combobox"], [tabindex]:not([tabindex="-1"])';
const POPUP_ROLES = ['listbox', 'menu', 'tree', 'grid', 'dialog'];
/** Everything in the dialog band, the platform's dialogs included. */
const DIALOGS = '.d4-dialog, .d4-dialog-center-wrapper, .u2-dialog, .u2-dialog-backdrop';
/** `--dg-z-dialog`: where the dialog band starts, for u2 and the platform alike. */
const DIALOG_Z = 3000;

/** The `u2-<n>` id counter, shared with `Control.rootId` and every engine copy on the page. */
const ids: {next: number} = (globalThis as Record<PropertyKey, any>)[Symbol.for('u2.ids')] ??= {next: 0};

/** The one overlay stack: every open popup and dialog, the last one on top. Escape and outside
 * clicks are decided here, by one listener set that exists only while something is open. Popups
 * live in a body-level host above everything; dialogs live in the dialog band they share with
 * the platform's `DG.Dialog`. */
export class Overlay {
  private static _host: HTMLElement | undefined;
  private static _layer = 0;
  private static readonly _stack: Entry[] = [];
  private static _listeners: Scope | undefined;

  /** Monotonic layer counter of the host: whatever opened there last paints on top. */
  static nextLayer(): number {
    return ++Overlay._layer;
  }

  /** The z-index that puts a dialog above every dialog on the page, `DG.Dialog` included: one more
   * than the highest in the band, and never under its floor. The platform raises its own dialogs
   * above u2's the same way, so whichever was opened or clicked last is on top. */
  static nextDialogZ(): number {
    let top = DIALOG_Z - 1;
    for (const el of document.querySelectorAll<HTMLElement>(DIALOGS)) {
      const z = parseInt(el.style.zIndex, 10);
      if (z > top)
        top = z;
    }
    return top + 1;
  }

  static get host(): HTMLElement {
    if (!Overlay._host || !Overlay._host.isConnected) {
      const el = document.createElement('div');
      el.className = 'u2-overlay';
      el.style.cssText = 'position: fixed; left: 0; top: 0; width: 0; height: 0; z-index: var(--dg-z-overlay);';
      document.body.append(el);
      // an empty layer reads as a leftover in a DOM dump, and every surface that puts something
      // here takes it away somewhere else — a popup's close, a balloon's timer. The host watching
      // itself is the one teardown all of them share; the getter builds a fresh one the moment
      // anything needs one again.
      new MutationObserver(() => {
        if (el.isConnected && el.childNodes.length === 0)
          el.remove();
      }).observe(el, {childList: true});
      Overlay._host = el;
    }
    return Overlay._host;
  }

  /** Shows `content` under `anchor` in the shared layer. The returned close is idempotent and
   * also registered on `scope`, so disposing the owner tears the overlay down; closing drops that
   * registration again, so an input reopened all day never piles dead closures up on its scope. */
  static show(anchor: HTMLElement, content: HTMLElement, scope: Scope,
    options: OverlayOptions = {}): () => void {
    // nothing to anchor to — and autoUpdate's first, synchronous update would close through
    // `position()` before the stop it returns is held, leaving its observers on for good
    if (!anchor.isConnected)
      return () => {};
    // The owner link: overlay content is portaled to the body host, so DOM ancestry cannot answer
    // "whose popup is this" — data-u2-owner carries the nearest named ancestor's automation id
    // instead, and hierarchical selectors follow that edge across the portal boundary.
    const owner = anchor.closest('[data-u2-name]') as HTMLElement | null;
    if (owner?.dataset.u2Name !== undefined)
      content.dataset.u2Owner = owner.dataset.u2Name;
    else
      delete content.dataset.u2Owner;
    const ownerControl = Control.forElement(anchor);
    if (ownerControl !== undefined)
      content.dataset.u2OwnerId = ownerControl.rootId;
    else
      delete content.dataset.u2OwnerId;
    content.id ||= `u2-${++ids.next}-popup`;
    content.dataset.u2Part = 'popup';
    const focus = Overlay._focusTarget(anchor);
    const controls = focus.getAttribute('aria-controls');
    focus.setAttribute('aria-controls', content.id);
    focus.setAttribute('aria-expanded', 'true');
    if (!focus.hasAttribute('aria-haspopup'))
      focus.setAttribute('aria-haspopup', POPUP_ROLES.find((r) => content.getAttribute('role') === r) ?? 'true');
    // Fixed strategy: viewport coordinates, independent of offsetParent chains — the app's
    // dock/view containers are positioned ancestors that break absolute-strategy math.
    content.style.position = 'fixed';
    content.style.left = '0';
    content.style.top = '0';
    content.style.zIndex = String(Overlay.nextLayer());
    Overlay.host.append(content);

    let closed = false;
    let stopAutoUpdate: () => void = () => {};
    const close = () => {
      if (closed)
        return;
      closed = true;
      stopAutoUpdate();
      Overlay._remove(entry);
      focus.setAttribute('aria-expanded', 'false');
      if (controls === null)
        focus.removeAttribute('aria-controls');
      else
        focus.setAttribute('aria-controls', controls);
      content.remove();
      scope.disown(close);
      content.dispatchEvent(new CustomEvent(OVERLAY_CLOSE_EVENT));
    };
    const entry: Entry = {content, close, anchor, focus};

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
              const available = `${Math.max(availableHeight, 64)}px`;
              const cap = options.maxHeight;
              elements.floating.style.maxHeight = cap === undefined ? available :
                `min(${available}, ${typeof cap === 'number' ? `${cap}px` : cap})`;
            },
          }),
        ],
      }).then(({x, y}) => {
        content.style.left = `${x}px`;
        content.style.top = `${y}px`;
      });
    };

    Overlay._push(entry);
    // event-driven: scroll, resize, element resize and layout shift of the anchor — no frame loop
    stopAutoUpdate = autoUpdate(anchor, content, position);
    scope.own(close);
    return close;
  }

  /** Puts a shown dialog on the stack; answers its removal. */
  static pushDialog(dialog: OverlayDialog): () => void {
    const entry: Entry = {content: dialog.root, close: () => dialog.close(), dialog};
    Overlay._push(entry);
    return () => Overlay._remove(entry);
  }

  /** Moves a stacked dialog above the other dialogs (the popups open above them stay on top):
   * the dialog the user brought forward is the one Escape reaches. */
  static raiseDialog(dialog: OverlayDialog): void {
    const stack = Overlay._stack;
    const i = stack.findIndex((e) => e.dialog === dialog);
    if (i < 0)
      return;
    const [entry] = stack.splice(i, 1);
    let j = stack.length - 1;
    while (j >= 0 && stack[j].dialog === undefined)
      j--;
    stack.splice(j + 1, 0, entry);
  }

  /** Whether a dialog of the platform (`.d4-dialog`) paints above `el`. */
  static platformDialogAbove(el: HTMLElement): boolean {
    const z = parseInt(el.style.zIndex, 10) || 0;
    return [...document.querySelectorAll<HTMLElement>('.d4-dialog')]
      .some((d) => (parseInt(d.style.zIndex, 10) || 0) > z);
  }

  private static _push(entry: Entry): void {
    Overlay._stack.push(entry);
    if (Overlay._listeners !== undefined)
      return;
    const scope = Overlay._listeners = new Scope(undefined, {label: 'overlay'});
    debug.markSingleton(scope);
    // capture: a popup owns its Escape before the control under it or the dialog around it sees one
    scope.listen(document, 'keydown', (e) => Scope.runWith(undefined, () => Overlay._onKeyDown(e as KeyboardEvent)),
      {capture: true});
    scope.listen(document, 'pointerdown', (e) => Scope.runWith(undefined, () => Overlay._onPointerDown(e)),
      {capture: true});
    // a window resize reflows the page under an anchored popup and strands it where its anchor was
    scope.listen(window, 'resize', () => Scope.runWith(undefined, () => Overlay._closePopups()));
    // autoUpdate reacts to scroll/resize/layout shift, not to detach; dock_spawn reparents constantly
    const observer = new MutationObserver(() => Overlay._sweep());
    observer.observe(document.body, {childList: true, subtree: true});
    scope.own(() => observer.disconnect());
  }

  private static _remove(entry: Entry): void {
    const stack = Overlay._stack;
    const i = stack.indexOf(entry);
    if (i < 0)
      return;
    stack.splice(i, 1);
    const listeners = Overlay._listeners;
    if (stack.length === 0 && listeners !== undefined) {
      Overlay._listeners = undefined;
      listeners.dispose();
    }
  }

  /** Escape goes to the top entry. A popup closes and the event stops there, so the control and
   * the dialog under it never see it. A dialog takes it from its own root (so a control inside can
   * claim it first); here only when the focus fell out onto the page and nothing is above it. */
  private static _onKeyDown(e: KeyboardEvent): void {
    if (e.key !== 'Escape' || e.defaultPrevented)
      return;
    const top = Overlay._stack[Overlay._stack.length - 1];
    if (top === undefined)
      return;
    const target = e.target instanceof Node ? e.target : null;
    const lost = Overlay._focusLost(target);
    if (top.dialog === undefined) {
      const own = target !== null && (top.content.contains(target) || top.anchor!.contains(target));
      // a platform dialog holding the focus keeps its Escape
      if (!own && !lost && target instanceof Element && target.closest('.d4-dialog') !== null)
        return;
      e.preventDefault();
      e.stopPropagation();
      // focus inside the popup would fall to the page with it: it goes back to the owner
      const active = document.activeElement;
      const refocus = active !== null && top.content.contains(active);
      top.close();
      if (refocus && top.focus!.isConnected)
        top.focus!.focus();
      return;
    }
    if (lost && top.dialog.modal && !Overlay.platformDialogAbove(top.dialog.root)) {
      e.preventDefault();
      e.stopPropagation();
      top.dialog.cancelByKey();
    }
  }

  /** Closes the popups above the entry the press landed in; dialogs never close on an outside press. */
  private static _onPointerDown(e: Event): void {
    const target = e.target as Node;
    for (const entry of [...Overlay._stack].reverse()) {
      if (entry.content.contains(target) || entry.anchor?.contains(target))
        return;
      if (entry.dialog === undefined)
        entry.close();
    }
  }

  private static _closePopups(): void {
    for (const entry of [...Overlay._stack].reverse()) {
      if (entry.dialog === undefined)
        entry.close();
    }
  }

  /** Closes what lost its place in the document: a popup whose anchor (or itself) was removed, a
   * dialog whose root was. */
  private static _sweep(): void {
    for (const entry of [...Overlay._stack]) {
      if (!entry.content.isConnected || (entry.anchor !== undefined && !entry.anchor.isConnected))
        entry.close();
    }
  }

  private static _focusLost(node: Node | Element | null): boolean {
    return node === null || node === document.body || node === document.documentElement || node === document;
  }

  /** The anchor, or the element in it that announces the popup or takes the focus. */
  private static _focusTarget(anchor: HTMLElement): HTMLElement {
    if (anchor.hasAttribute('aria-haspopup') || anchor.matches(FOCUSABLE))
      return anchor;
    return anchor.querySelector<HTMLElement>('[aria-haspopup]') ?? anchor.querySelector<HTMLElement>(FOCUSABLE) ??
      anchor;
  }
}
