/* Dialog. Renders into the shared Overlay.host layer, so one place decides layering; each show
   takes the next z offset, which is what makes stacking work. The fluent surface mirrors
   DG.Dialog so platform users read it unaided. */

import {signal, ReadonlySignal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {Scope} from '../../core/scope.js';
import {Overlay} from '../../core/overlay.js';
import {button} from '../../core/elements.js';

export interface DialogShowOptions {
  modal?: boolean;
  x?: number;
  y?: number;
  width?: number;
  height?: number;
}

const FOCUSABLE = 'a[href], button:not([disabled]), input:not([disabled]), select:not([disabled]), ' +
  'textarea:not([disabled]), [tabindex]:not([tabindex="-1"])';

export class Dialog extends Control {
  private static _seq = 0;

  readonly isOpen: ReadonlySignal<boolean>;

  private readonly _open = signal(false);
  private readonly _bar = document.createElement('div');
  private readonly _content = document.createElement('div');
  private readonly _buttons = document.createElement('div');
  private readonly _backdrop = document.createElement('div');
  private _ok: HTMLButtonElement | undefined;
  private _cancel: HTMLButtonElement | undefined;
  private _onOK: (() => unknown) | undefined;
  private _onCancel: (() => void) | undefined;
  private _scope: Scope | undefined;
  private _restore: HTMLElement | undefined;
  private _drag: {dx: number, dy: number} | undefined;

  constructor(title: string, options?: {name?: string}) {
    super();
    this.isOpen = this._open;
    if (options?.name !== undefined)
      this.name = options.name;
    this.root.classList.add('u2-dialog');
    this.root.dataset.u2 = 'dialog';
    this.root.tabIndex = -1;
    this.root.setAttribute('role', 'dialog');

    const id = `u2-dialog-${++Dialog._seq}-title`;
    const text = document.createElement('div');
    text.className = 'u2-dialog-title-text';
    text.id = id;
    text.textContent = title;
    this.root.setAttribute('aria-labelledby', id);

    const close = document.createElement('button');
    close.type = 'button';
    close.className = 'u2-dialog-close';
    close.textContent = '✕';
    close.setAttribute('aria-label', 'Close');

    this._bar.className = 'u2-dialog-title';
    this._bar.append(text, close);
    this._content.className = 'u2-dialog-content';
    this._buttons.className = 'u2-dialog-buttons';
    this._backdrop.className = 'u2-dialog-backdrop';
    this.root.append(this._bar, this._content, this._buttons);

    this._listen(close, 'click', () => this._finish(this._onCancel));
    this._listen(this._bar, 'pointerdown', (e) => this._onPointerDown(e as PointerEvent));
    this._listen(this._bar, 'pointermove', (e) => this._onPointerMove(e as PointerEvent));
    this._listen(this._bar, 'pointerup', () => this._endDrag());
    this._listen(this._bar, 'pointercancel', () => this._endDrag());
    this._listen(this._backdrop, 'pointerdown', (e) => e.preventDefault());
  }

  static create(title: string, options?: {name?: string}): Dialog {
    return new Dialog(title, options);
  }

  add(el: HTMLElement | Control | string): Dialog {
    this._content.append(Control.is(el) ? el.root : el);
    return this;
  }

  addButton(text: string, onClick: () => void, options?: {primary?: boolean}): Dialog {
    this._buttons.append(this._button(text, onClick, options?.primary));
    return this;
  }

  /** `fn` answering `false` keeps the dialog open — for a pick that is not one yet. */
  onOK(fn: () => unknown): Dialog {
    this._onOK = fn;
    this._okButton();
    return this;
  }

  /** OK — the button and Enter alike — follows `on`: a pick the dialog cannot finish without. */
  okEnabled(on: ReadonlySignal<boolean>): Dialog {
    const ok = this._okButton();
    this.effect(() => ok.disabled = !on.value);
    return this;
  }

  private _okButton(): HTMLButtonElement {
    if (!this._ok) {
      this._ok = this._button('OK', () => this._finish(this._onOK), true);
      this._buttons.prepend(this._ok);
    }
    return this._ok;
  }

  onCancel(fn: () => void): Dialog {
    this._onCancel = fn;
    if (!this._cancel) {
      this._cancel = this._button('CANCEL', () => this._finish(this._onCancel));
      this._buttons.insertBefore(this._cancel, this._ok ? this._ok.nextSibling : this._buttons.firstChild);
    }
    return this;
  }

  show(options: DialogShowOptions = {}): Dialog {
    if (this._scope)
      return this;
    const scope = new Scope();
    this._scope = scope;
    this._restore = (document.activeElement as HTMLElement | null) ?? undefined;

    const backdropLayer = Overlay.nextLayer();
    this.root.style.zIndex = String(Overlay.nextLayer());
    if (options.width !== undefined)
      this.root.style.width = `${options.width}px`;
    if (options.height !== undefined)
      this.root.style.height = `${options.height}px`;
    if (options.modal) {
      this._backdrop.style.zIndex = String(backdropLayer);
      this.root.setAttribute('aria-modal', 'true');
      Overlay.host.append(this._backdrop);
    } else
      this.root.removeAttribute('aria-modal');

    Overlay.host.append(this.root);
    this._open.value = true;
    this._place(options.x, options.y);
    this._listen(this.root, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent), scope);
    const first = Dialog._focusable(this._content)[0];
    (first ?? this._ok ?? this.root).focus();
    return this;
  }

  close(): void {
    const scope = this._scope;
    if (!scope)
      return;
    this._scope = undefined;
    this._drag = undefined;
    scope.dispose();
    this._backdrop.remove();
    this.root.remove();
    this._open.value = false;
    const restore = this._restore;
    this._restore = undefined;
    if (restore && restore.isConnected)
      restore.focus();
  }

  dispose(): void {
    this.close();
    super.dispose();
  }

  private _button(text: string, onClick: () => void, primary?: boolean): HTMLButtonElement {
    return this.runInScope(() => button(text, onClick, {primary}));
  }

  private _finish(fn: (() => unknown) | undefined): void {
    if (fn?.() === false)
      return;
    this.close();
  }

  private _place(x?: number, y?: number): void {
    const rect = this.root.getBoundingClientRect();
    const left = x ?? (window.innerWidth - rect.width) / 2;
    const top = y ?? (window.innerHeight - rect.height) / 2;
    this.root.style.left = `${Math.max(0, Math.min(left, window.innerWidth - rect.width))}px`;
    this.root.style.top = `${Math.max(0, Math.min(top, window.innerHeight - rect.height))}px`;
  }

  private _onPointerDown(e: PointerEvent): void {
    if (e.button !== 0 || (e.target as HTMLElement).closest('.u2-dialog-close'))
      return;
    e.preventDefault();
    const rect = this.root.getBoundingClientRect();
    this._drag = {dx: e.clientX - rect.left, dy: e.clientY - rect.top};
    this._bar.setPointerCapture(e.pointerId);
  }

  private _onPointerMove(e: PointerEvent): void {
    if (this._drag)
      this._place(e.clientX - this._drag.dx, e.clientY - this._drag.dy);
  }

  private _endDrag(): void {
    this._drag = undefined;
  }

  private _onKeyDown(e: KeyboardEvent): void {
    const target = e.target as HTMLElement;
    if (e.key === 'Escape') {
      e.preventDefault();
      this._finish(this._onCancel);
    } else if (e.key === 'Enter' && this._ok && target.tagName !== 'TEXTAREA' && target.tagName !== 'BUTTON') {
      e.preventDefault();
      if (this._ok.disabled)
        return;
      this._finish(this._onOK);
    } else if (e.key === 'Tab')
      this._trapTab(e);
  }

  private _trapTab(e: KeyboardEvent): void {
    const els = Dialog._focusable(this.root);
    if (els.length === 0)
      return;
    const active = document.activeElement;
    if (e.shiftKey && active === els[0]) {
      e.preventDefault();
      els[els.length - 1].focus();
    } else if (!e.shiftKey && active === els[els.length - 1]) {
      e.preventDefault();
      els[0].focus();
    }
  }

  private _listen(target: EventTarget, type: string, handler: (e: Event) => void, scope?: Scope): void {
    target.addEventListener(type, handler);
    (scope ?? this.scope).own(() => target.removeEventListener(type, handler));
  }

  private static _focusable(el: HTMLElement): HTMLElement[] {
    const els = Array.from(el.querySelectorAll(FOCUSABLE)) as HTMLElement[];
    return els.filter((x) => x.offsetParent !== null);
  }
}
