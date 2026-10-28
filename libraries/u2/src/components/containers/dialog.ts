/* Dialog. Lives in the dialog band it shares with the platform's DG.Dialog (z-index from
   Overlay.nextDialogZ, raised on press), so whichever was opened or clicked last is on top; popups
   opened from it stay above, in the overlay host. The overlay stack routes Escape. Closing disposes
   the dialog and what it holds unless `keepAlive`; shown, it no longer dies with the scope that built
   it (a row, an effect run). The fluent surface mirrors DG.Dialog so platform users read it unaided. */

import {signal, computed, batch, ReadonlySignal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {Scope} from '../../core/scope.js';
import {Ownership} from '../../core/owner.js';
import type {Owner} from '../../core/owner.js';
import {Overlay, OverlayDialog} from '../../core/overlay.js';
import {Tooltip} from '../../core/tooltip.js';
import {button} from '../../core/elements.js';
import type {ButtonOptions} from '../../core/elements.js';
import {Reason} from '../../core/reason.js';
import {icon} from '../display/icon.js';
import {Form} from '../forms/form.js';
import {statusLine} from '../display/badge.js';

export interface DialogOptions {
  /** The automation id (`data-u2-name`). */
  name?: string;
  /** Keep the dialog and its content alive after close, to show it again; it then stays owned by the
   * scope that built it. Default: closing disposes. */
  keepAlive?: boolean;
  /** One way out: a CLOSE button, never a CANCEL — for a dialog whose work applies as it is done (a tracker),
   * where there is nothing to take back. CLOSE, ✕ and Escape all close it as OK does (`onOK` runs if given). */
  closeOnly?: boolean;
}

export interface DialogShowOptions {
  modal?: boolean;
  x?: number;
  y?: number;
  width?: number;
  /** A fixed height: the body fills it, the footer stays at the bottom. */
  height?: number;
  /** The user resizes it by its bottom-right corner; the body fills whatever size it has. */
  resizable?: boolean;
  /** Hidden — not closed — while false, a modal's backdrop too: a dialog pinned to a view hides while another
   * view is current (u2/dg `ViewHost.show`). */
  visible?: ReadonlySignal<boolean>;
  /** Closes the dialog when it is disposed: a dialog that must not outlive a view or a session. */
  owner?: Owner;
}

/** Why a dialog closed: OK (the button or Enter), cancel (CANCEL, ✕ or Escape), or `close()`. */
export type DialogCloseReason = 'ok' | 'cancel' | 'code';

/** What an OK handler gets: `signal` aborts when the user presses Escape or ✕ while it runs. */
export interface DialogOkContext {
  signal: AbortSignal;
}

const FOCUSABLE = 'a[href], button:not([disabled]), input:not([disabled]), select:not([disabled]), ' +
  'textarea:not([disabled]), [tabindex]:not([tabindex="-1"])';

/** Where an Enter is the control's own — a button, a link, a box, a select, a textarea, a chip, or
 * anything playing one of those — never the dialog's default action. */
const ENTER_OWNERS = 'button, a, input[type="checkbox"], input[type="radio"], select, textarea, .u2-chip, ' +
  '.u2-chips-chip, [role="button"], [role="checkbox"], [role="option"]';

export class Dialog extends Control implements OverlayDialog {
  static isRoot = true;
  /** Why the footer waits while an OK runs. */
  static readonly WAITING = 'Wait for OK to finish';
  private static _seq = 0;

  readonly isOpen: ReadonlySignal<boolean>;
  /** An async OK is running: the footer is disabled and OK shows a spinner. */
  readonly busy: ReadonlySignal<boolean>;
  /** The error row: the message of the last OK that threw or rejected; null when there is none. */
  readonly error: ReadonlySignal<string | null>;

  private readonly _open = signal(false);
  private readonly _busy = signal(false);
  private readonly _error = signal<string | null>(null);
  private readonly _title = document.createElement('div');
  private readonly _bar = document.createElement('div');
  private readonly _closeBtn = Dialog._titleButton('u2-dialog-close', 'close', 'times', 'Close');
  private _helpBtn: HTMLButtonElement | undefined;
  private _onHelp: (() => unknown) | undefined;
  private readonly _content = document.createElement('div');
  private readonly _errorRow = document.createElement('div');
  private readonly _buttons = document.createElement('div');
  private readonly _backdrop = document.createElement('div');
  private readonly _spinner = document.createElement('span');
  private readonly _keepAlive: boolean;
  private readonly _onClose: ((reason: DialogCloseReason) => void)[] = [];
  private _okBtn: HTMLButtonElement | undefined;
  private _cancelBtn: HTMLButtonElement | undefined;
  private _reasonLine: HTMLElement | undefined;
  /** The footer line's quiet text, while OK can run: what it will do. */
  private readonly _quiet = signal<ReadonlySignal<string | null> | string | null>(null);
  private readonly _okOff = signal<Reason | undefined>(undefined);
  /** The form in the body, found on show: OK waits for its verdict, as the platform's does. */
  private readonly _form = signal<Form | undefined>(undefined);
  /** Why OK is off; null while it is on. */
  private readonly _okReason = computed<string | null>(() =>
    // the form's fields first: the reasons read top to bottom, as the user fills the dialog in
    this._form.value?.problem.value ?? Reason.read(this._okOff.value));
  private readonly _closeOnly: boolean;
  private _onOK: ((ctx: DialogOkContext) => unknown) | undefined;
  private _onCancel: (() => unknown) | undefined;
  private _closeGuard: (() => boolean) | undefined;
  /** Set while shown: takes the dialog off the overlay stack. */
  private _unstack: (() => void) | undefined;
  /** Set while an async OK runs. */
  private _pending: AbortController | undefined;
  /** The footer buttons a busy OK disabled, enabled again when it settles. */
  private _held: HTMLButtonElement[] = [];
  /** The reason a handler in flight closes with, when it closes or disposes the dialog itself. */
  private _closing: DialogCloseReason | undefined;
  private _restore: HTMLElement | undefined;
  private _drag: {dx: number, dy: number} | undefined;
  /** Where the dialog was put — an axis left undefined stays centred, through window resizes too. */
  private _at: {x?: number, y?: number} = {};
  private _modal = false;
  /** What one showing owns (the `visible` follow, the link to an `owner`): disposed on close. */
  private _showing: Scope | undefined;

  constructor(title: string, options: DialogOptions = {}) {
    super();
    this.isOpen = this._open;
    this.busy = this._busy;
    this.error = this._error;
    this._keepAlive = options.keepAlive === true;
    this._closeOnly = options.closeOnly === true;
    if (options.name !== undefined)
      this.name = options.name;
    this.root.classList.add('u2-dialog');
    this.root.dataset.u2 = 'dialog';
    this.root.tabIndex = -1;
    this.root.setAttribute('role', 'dialog');

    const id = `u2-dialog-${++Dialog._seq}-title`;
    this._title.className = 'u2-dialog-title-text';
    this._title.id = id;
    this._title.dataset.u2Part = 'title';
    this._title.textContent = title;
    this.root.setAttribute('aria-labelledby', id);

    const close = this._closeBtn;
    this._bar.className = 'u2-dialog-title';
    this._bar.append(this._title, close);
    this._content.className = 'u2-dialog-content';
    this._content.dataset.u2Part = 'body';
    this._errorRow.className = 'u2-dialog-error';
    this._errorRow.dataset.u2Part = 'error-row';
    this._errorRow.setAttribute('role', 'alert');
    this._errorRow.hidden = true;
    this._buttons.className = 'u2-dialog-buttons';
    this._buttons.dataset.u2Part = 'footer';
    this._backdrop.className = 'u2-dialog-backdrop';
    this._spinner.className = 'u2-dialog-busy';
    this._spinner.dataset.u2Part = 'busy';
    this._spinner.setAttribute('aria-hidden', 'true');
    this.root.append(this._bar, this._content, this._errorRow, this._buttons);

    const scope = this.scope;
    scope.listen(close, 'click', () => this.cancelByKey());
    scope.listen(this._bar, 'pointerdown', (e) => this._onPointerDown(e as PointerEvent));
    scope.listen(this._bar, 'pointermove', (e) => this._onPointerMove(e as PointerEvent));
    scope.listen(this._bar, 'pointerup', () => this._endDrag());
    scope.listen(this._bar, 'pointercancel', () => this._endDrag());
    scope.listen(this._backdrop, 'pointerdown', (e) => e.preventDefault());
    scope.listen(this.root, 'pointerdown', () => this._raise());
    scope.listen(this.root, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    scope.listen(this._content, 'input', () => this._error.value = null);
    scope.listen(this._content, 'change', () => this._error.value = null);
    const refit = () => {
      if (this._unstack !== undefined)
        this._place(this._at.x, this._at.y);
    };
    scope.listen(window, 'resize', refit);
    const resize = new ResizeObserver(refit);
    resize.observe(this.root);
    this.own(() => resize.disconnect());
    this.effect(() => {
      const message = this._error.value;
      this._errorRow.textContent = message ?? '';
      this._errorRow.hidden = message === null;
    });
    if (this._closeOnly)
      this._okButton('CLOSE');
  }

  static create(title: string, options?: DialogOptions): Dialog {
    return new Dialog(title, options);
  }

  /** Shown modal and not hidden by `visible`: a hidden dialog blocks nothing and takes no Escape. */
  get modal(): boolean {
    return this._modal && !this.root.hidden;
  }

  /** Content handed to the dialog is the dialog's: a control built where the dialog is being filled
   * (adopted by the scope building it, or by none), and every such control inside an added element,
   * moves to the dialog, so closing frees it. Plain elements' handlers and live text stay with the
   * scope that built them: content built in a short-lived scope (a row, an effect run) is better
   * passed as `add(() => …)`, which builds it owned by the dialog. */
  add(el: HTMLElement | Control | string | (() => HTMLElement | Control | string)): this {
    const built = typeof el === 'function' ? this.build(el) : el;
    const node = Control.is(built) ? built.root : built;
    if (typeof node !== 'string')
      this._adopt(node);
    this._content.append(node);
    return this;
  }

  /** A footer button. `onClick` gets the event and the dialog and runs owned by whoever added the
   * button; it does not close the dialog. `disabled` is `button()`'s: off while it says why. `start`
   * puts it at the footer's start, away from the outcomes: a tool of the dialog ("Modify Stages"). */
  addButton(text: string, onClick: (e: MouseEvent, dialog: Dialog) => unknown,
    options?: {primary?: boolean, disabled?: ButtonOptions['disabled'], start?: boolean}): this {
    const run = Dialog._owned(onClick);
    const el = this._button(text, (e) => run(e, this), options?.primary, options?.disabled);
    if (options?.start === true) {
      el.classList.add('u2-dialog-start');
      this._buttons.insertBefore(el, Dialog._firstAction(this._buttons, true));
    } else
      this._buttons.insertBefore(el, this._okBtn ?? null);
    this._cancelButton();
    return this;
  }

  /** The ? in the title bar, beside ✕: the platform dialog's help. */
  onHelp(fn: () => unknown): this {
    this._onHelp = Dialog._owned(fn);
    if (!this._helpBtn) {
      this._helpBtn = Dialog._titleButton('u2-dialog-help', 'help', 'question', 'Help');
      this.scope.listen(this._helpBtn, 'click', () => this._onHelp?.());
      this._bar.insertBefore(this._helpBtn, this._closeBtn);
    }
    return this;
  }

  /** A button by its caption, case-insensitive: the footer's, OK and CANCEL included, and the ✕
   * (`Close`) and ? (`Help`). */
  getButton(text: string): HTMLButtonElement | undefined {
    const wanted = text.trim().toLowerCase();
    return [...this._bar.querySelectorAll<HTMLButtonElement>('button'),
      ...this._buttons.querySelectorAll<HTMLButtonElement>('button')].find((b) =>
      (b.textContent ?? '').trim().toLowerCase() === wanted || b.getAttribute('aria-label')?.toLowerCase() === wanted);
  }

  /** The OK action — the button and Enter alike. Answering `false` keeps the dialog open (a pick
   * that is not one yet). A promise keeps it open and busy until it settles: it closes unless the
   * promise answers `false`; a rejection, like a throw, stays open with the message in the error
   * row. Escape or ✕ while busy aborts `signal`. `text` is the button's caption; `primary` raises
   * it, for an OK that commits something. */
  onOK(fn: (ctx: DialogOkContext) => unknown, options?: {text?: string, primary?: boolean}): this {
    this._onOK = Dialog._owned(fn);
    const ok = this._okButton(options?.text);
    if (options?.primary !== undefined)
      ok.classList.toggle('u2-btn-primary', options.primary);
    return this;
  }

  /** OK — the button and Enter alike — is off while `reason` says why: it becomes the button's tooltip, so
   * the user learns what is missing (`'Select at least one monomer'`); `null` leaves it on. A form in the body
   * adds its own: "Name: Value can't be empty". `footer` also shows the reason as a warning beside the buttons. */
  okDisabled(reason: Reason, options?: {footer?: boolean}): this {
    this._okButton();
    this._okOff.value = reason;
    if (options?.footer === true)
      this._footerLine();
    return this;
  }

  /** A quiet line at the footer's start saying what the dialog's action will do ("Submits 260 molecules to Demo
   * File Submit"); while OK is off, the line says why instead, as `okDisabled(reason, {footer: true})` does. */
  status(text: ReadonlySignal<string | null> | string): this {
    this._quiet.value = text;
    this._footerLine().classList.add('u2-dialog-status');
    return this;
  }

  private _footerLine(): HTMLElement {
    if (this._reasonLine === undefined) {
      const quiet = computed(() => {
        const text = this._quiet.value;
        return text === null || typeof text === 'string' ? text : text.value;
      });
      const line = this._reasonLine = this.build(() =>
        statusLine(computed(() => this._okReason.value || null), quiet));
      line.classList.add('u2-dialog-reason', 'u2-status-line-warning');
      line.dataset.u2Part = 'reason';
      this._buttons.insertBefore(line, Dialog._firstAction(this._buttons, false));
    }
    return this._reasonLine;
  }

  /** Where the footer's actions begin: after the start buttons (and, unless `beforeLine`, the status line). */
  private static _firstAction(footer: HTMLElement, beforeLine: boolean): Element | null {
    return Array.from(footer.children).find((el) => !el.classList.contains('u2-dialog-start') &&
      (beforeLine || !el.classList.contains('u2-dialog-reason'))) ?? null;
  }

  /** Vetoes every close — ✕, Esc and the buttons — while `fn` answers false: work in flight the
   * dialog must outlive. */
  closeGuard(fn: () => boolean): this {
    this._closeGuard = fn;
    return this;
  }

  /** CANCEL, ✕ and Escape; answering `false` keeps the dialog open. */
  onCancel(fn: () => void): this {
    this._onCancel = Dialog._owned(fn);
    this._cancelButton();
    return this;
  }

  /** CANCEL comes with the first action, as in every platform dialog (`dialog.dart` commandBar.onCancel); a
   * close-only dialog has none. */
  private _cancelButton(): void {
    if (this._cancelBtn || this._closeOnly)
      return;
    this._cancelBtn = this._button('CANCEL', () => this._dismiss());
    this._buttons.insertBefore(this._cancelBtn, Dialog._firstAction(this._buttons, false));
  }

  /** Runs on every close, with the reason. */
  onClose(fn: (reason: DialogCloseReason) => void): this {
    this._onClose.push(Dialog._owned(fn));
    return this;
  }

  show(options: DialogShowOptions = {}): this {
    if (this._unstack !== undefined)
      return this;
    if (this.scope.isDisposed) {
      console.warn(`u2: dialog "${this._title.textContent}" was disposed when it closed — ` +
        'create it with {keepAlive: true} to show it again');
      return this;
    }
    const adopter = this._u2.adopter;
    if (adopter !== undefined && !this._keepAlive) {
      this._u2.adopter = undefined;
      adopter.scope.disown(adopter.release);
    }
    this._restore = (document.activeElement as HTMLElement | null) ?? undefined;
    this._error.value = null;
    if (options.width !== undefined)
      this.root.style.width = `${options.width}px`;
    if (options.height !== undefined)
      this.root.style.height = `${options.height}px`;
    // under a fixed or a user-set height the body fills the dialog instead of stopping at its own cap
    this.root.classList.toggle('u2-dialog-sized', options.height !== undefined || options.resizable === true);
    this.root.classList.toggle('u2-dialog-resizable', options.resizable === true);
    this._modal = options.modal === true;
    const z = Overlay.nextDialogZ();
    if (this._modal) {
      this._backdrop.style.zIndex = String(z);
      this.root.style.zIndex = String(z + 1);
      this.root.setAttribute('aria-modal', 'true');
      document.body.append(this._backdrop);
    } else {
      this.root.style.zIndex = String(z);
      this.root.removeAttribute('aria-modal');
    }
    document.body.append(this.root);
    this._unstack = Overlay.pushDialog(this);
    this._open.value = true;
    this._form.value = Dialog._formOf(this._content);
    this._alignLabels();
    this._at = {x: options.x, y: options.y};
    this._place(options.x, options.y);
    if (options.resizable === true) {
      // the corner changes the size, not the place: a dialog kept centred would slide under the pointer
      const rect = this.root.getBoundingClientRect();
      this._at = {x: rect.left, y: rect.top};
    }
    const first = Dialog._focusable(this._content).find((el) => !Dialog._readOnly(el));
    (first ?? this._okBtn ?? this.root).focus();
    if (first instanceof HTMLInputElement || first instanceof HTMLTextAreaElement)
      first.select();
    this._pin(options.visible, options.owner);
    return this;
  }

  /** What this showing follows: `visible` hides it, a disposed `owner` closes it; both end with it. */
  private _pin(visible?: ReadonlySignal<boolean>, owner?: Owner): void {
    if (visible === undefined && owner === undefined)
      return;
    const showing = this._showing = this.scope.child('dialog-showing');
    if (visible !== undefined)
      showing.effect(() => this._hide(!visible.value));
    if (owner !== undefined)
      showing.own(Ownership.adopt(owner, () => this.close()));
  }

  private _hide(hidden: boolean): void {
    if (this.root.hidden === hidden)
      return;
    this.root.hidden = hidden;
    this._backdrop.hidden = hidden;
    if (!hidden)
      this._place(this._at.x, this._at.y);
  }

  /** A field shown for reading — a native `readOnly`, `aria-readonly`, or inside a read-only u2 input: the first
   * focus passes it by, so a dialog opening on a summary does not select the summary's text. */
  private static _readOnly(el: HTMLElement): boolean {
    return (el as HTMLInputElement).readOnly === true || el.getAttribute('aria-readonly') === 'true' ||
      el.closest('.u2-input-readonly') !== null;
  }

  /** Hides the dialog and, unless `keepAlive`, disposes it with its content. A pending OK is aborted. */
  close(reason: DialogCloseReason = 'code'): void {
    const unstack = this._unstack;
    if (unstack === undefined)
      return;
    this._unstack = undefined;
    unstack();
    this._pending?.abort();
    this._setBusy(undefined);
    this._drag = undefined;
    this._showing?.dispose();
    this._showing = undefined;
    this._hide(false);
    this._backdrop.remove();
    this.root.remove();
    this._open.value = false;
    const restore = this._restore;
    this._restore = undefined;
    // a dialog still open below takes focus unless what held it is inside that dialog: the page under
    // it (a ribbon icon cannot take focus, so the opener reads as the body) would leave Escape nowhere
    const below = Dialog._topmost();
    if (below !== undefined && !(restore?.isConnected && below.contains(restore)))
      below.focus();
    else if (restore?.isConnected)
      restore.focus();
    const why = this._closing ?? reason;
    for (const fn of this._onClose) {
      try {
        fn(why);
      } catch (e) {
        this.scope.reportError(e);
      }
    }
    this.fireEvent('close', why);
    if (!this._keepAlive)
      this.dispose();
  }

  dispose(): void {
    this.close();
    super.dispose();
  }

  /** Escape and ✕: aborts a pending OK, else cancels — in a close-only dialog, closes as OK does. */
  cancelByKey(): void {
    if (this._pending !== undefined)
      this._pending.abort();
    else if (this._closeOnly)
      this._accept();
    else
      this._dismiss();
  }

  /** The base state plus whether it is open, busy, the error row and the validity of a form in it. */
  getState(): ReturnType<Control['getState']> {
    return {...super.getState(), open: this._open.peek(), busy: this._busy.peek(), error: this._error.peek(),
      valid: Dialog._formOf(this._content)?.validity?.peek() == null};
  }

  private static _formOf(content: HTMLElement): Form | undefined {
    const formEl = content.querySelector<HTMLElement>('.u2-form');
    const form = formEl === null ? undefined : Control.forElement(formEl);
    return form instanceof Form ? form : undefined;
  }

  /** The title — read off the element, so a dialog not shown (not in the document) still has it. */
  protected accessibleLabel(): string | undefined {
    return Control.text(this._title) || undefined;
  }

  /** The platform's footer: CANCEL first, OK right-most, both flat. */
  private _okButton(text?: string): HTMLButtonElement {
    if (!this._okBtn) {
      const ok = this._okBtn = this._button(text ?? 'OK', () => this._accept());
      Tooltip.bind(ok, () => Reason.of(ok) ?? '', {owner: this.scope});
      this._buttons.append(ok);
      this._cancelButton();
      this.effect(() => {
        const reason = this._okReason.value;
        if (!this._busy.peek())
          Reason.apply(ok, reason);
      });
    } else if (text !== undefined)
      this._okBtn.textContent = text;
    return this._okBtn;
  }

  private _button(text: string, onClick: (e: MouseEvent) => void, primary?: boolean,
    disabled?: ButtonOptions['disabled']): HTMLButtonElement {
    return this.build(() => button(text, onClick, {primary, disabled}));
  }

  private _accept(): void {
    if (this._pending !== undefined || this._closeGuard?.() === false)
      return;
    this._error.value = null;
    const abort = new AbortController();
    let answer: unknown;
    this._closing = 'ok';
    try {
      answer = this._onOK?.({signal: abort.signal});
    } catch (e) {
      this._fail(e);
      return;
    } finally {
      this._closing = undefined;
    }
    if (this._unstack === undefined)
      return;
    const work = Scope.thenable(answer);
    if (work !== undefined)
      void this._await(work, abort);
    else if (answer !== false)
      this.close('ok');
  }

  private async _await(work: PromiseLike<unknown>, abort: AbortController): Promise<void> {
    this._setBusy(abort);
    let answer: unknown;
    try {
      answer = await work;
    } catch (e) {
      if (this._pending === abort) {
        this._setBusy(undefined);
        this._fail(e);
      }
      return;
    }
    if (this._pending !== abort)
      return;
    this._setBusy(undefined);
    if (answer !== false)
      this.close('ok');
  }

  private _dismiss(): void {
    if (this._closeGuard?.() === false)
      return;
    let answer: unknown;
    this._closing = 'cancel';
    try {
      answer = this._onCancel?.();
    } catch (e) {
      this.scope.reportError(e);
    } finally {
      this._closing = undefined;
    }
    if (answer !== false)
      this.close('cancel');
  }

  /** An aborted OK shows nothing: the user asked for it. */
  private _fail(e: unknown): void {
    if ((e as {name?: unknown} | null)?.name !== 'AbortError')
      this._error.value = e instanceof Error ? e.message : String(e);
  }

  private _setBusy(abort: AbortController | undefined): void {
    this._pending = abort;
    const busy = abort !== undefined;
    if (this._busy.peek() === busy)
      return;
    this._busy.value = busy;
    const ok = this._okBtn;
    if (busy) {
      // the dialog holds the focus while busy, not a footer button a key could press again
      if (this._buttons.contains(document.activeElement))
        this.root.focus();
      this._held = [...this._buttons.querySelectorAll<HTMLButtonElement>('button')].filter((b) => !b.disabled);
      for (const b of this._held)
        Reason.apply(b, Dialog.WAITING);
      this.root.setAttribute('aria-busy', 'true');
      ok?.setAttribute('aria-busy', 'true');
      ok?.append(this._spinner);
      return;
    }
    for (const b of this._held)
      Reason.apply(b, null);
    this._held = [];
    if (ok !== undefined)
      Reason.apply(ok, this._okReason.peek());
    this.root.removeAttribute('aria-busy');
    ok?.removeAttribute('aria-busy');
    this._spinner.remove();
  }

  private _adopt(node: HTMLElement): void {
    const builder = Scope.ambient;
    if (builder === this.scope)
      return;
    for (const el of [node, ...node.querySelectorAll<HTMLElement>('.u2-component')]) {
      const control = Control.forElement(el);
      if (control === undefined || control.root !== el || control === this)
        continue;
      const u2 = control._u2;
      const previous = u2.adopter;
      // unowned and built with no scope around — unless it never joins one (a platform DG.Widget)
      const orphan = previous === undefined && builder === undefined &&
        (control as unknown as {ambientAdoption: boolean}).ambientAdoption;
      if (!orphan && (previous === undefined || previous.scope !== builder))
        continue;
      if (previous !== undefined)
        previous.scope.disown(previous.release);
      const release = () => control.dispose();
      u2.adopter = {scope: this.scope, release};
      this.scope.own(release);
    }
  }

  /** A press anywhere on the dialog brings it above the other dialogs, the platform's included. */
  private _raise(): void {
    if (this._unstack === undefined)
      return;
    const next = Overlay.nextDialogZ();
    if (parseInt(this.root.style.zIndex, 10) < next - 1)
      this.root.style.zIndex = String(next);
    Overlay.raiseDialog(this);
  }

  private _place(x?: number, y?: number): void {
    // a fixed box shrinks to the room right of its left edge: measured there, a dialog dragged
    // right would keep the width the window squeezed it to
    this.root.style.left = '0px';
    const rect = this.root.getBoundingClientRect();
    const left = x ?? (window.innerWidth - rect.width) / 2;
    const top = y ?? (window.innerHeight - rect.height) / 2;
    this.root.style.left = `${Math.max(0, Math.min(left, window.innerWidth - rect.width))}px`;
    this.root.style.top = `${Math.max(0, Math.min(top, window.innerHeight - rect.height))}px`;
  }

  private _onPointerDown(e: PointerEvent): void {
    if (e.button !== 0 || (e.target as HTMLElement).closest('button'))
      return;
    e.preventDefault();
    const rect = this.root.getBoundingClientRect();
    this._drag = {dx: e.clientX - rect.left, dy: e.clientY - rect.top};
    this._bar.setPointerCapture(e.pointerId);
  }

  private _onPointerMove(e: PointerEvent): void {
    if (!this._drag)
      return;
    this._at = {x: e.clientX - this._drag.dx, y: e.clientY - this._drag.dy};
    this._place(this._at.x, this._at.y);
  }

  /** Captions straight in the body share one column, as wide as the longest of them up to a cap;
   * a caption cut there keeps its full text as the tooltip, and a choice gets the width a form
   * would give it. A Form keeps its own column. */
  private _alignLabels(): void {
    const rows = Array.from(this._content.children)
      .filter((row): row is HTMLElement => row.classList.contains('u2-input-root'));
    Form.fitChoices(this._content, rows);
    const labels = rows
      .map((row) => row.querySelector<HTMLElement>('.u2-input-label'))
      .filter((label): label is HTMLElement => label !== null);
    if (labels.length === 0)
      return;
    this._content.style.setProperty('--u2-dialog-label-width', 'auto');
    let width = 0;
    for (const label of labels)
      width = Math.max(width, Math.ceil(label.getBoundingClientRect().width));
    this._content.style.setProperty('--u2-dialog-label-width', `${width}px`);
    for (const label of labels)
      Tooltip.titleIfCut(label);
  }

  private _endDrag(): void {
    this._drag = undefined;
  }

  /** Escape here, from inside the dialog, after the control holding the focus had its say; a popup
   * open above took its Escape before (the overlay stack). */
  private _onKeyDown(e: KeyboardEvent): void {
    if (e.key === 'Escape') {
      if (e.defaultPrevented)
        return;
      e.preventDefault();
      e.stopPropagation();
      this.cancelByKey();
    } else if (this._okBtn && Dialog.entersDefault(e)) {
      e.preventDefault();
      // an OK held off never runs: Enter shows every field's verdict instead, or — the fields all
      // valid — focuses OK, whose tooltip says why
      if (!this._okBtn.disabled)
        this._accept();
      else if (!this._busy.peek() && this._form.peek()?.validate() !== false)
        this._okBtn.focus();
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

  /** The open dialog on top of the others — where focus goes back when the one that held it is gone. */
  private static _topmost(): HTMLElement | undefined {
    let top: HTMLElement | undefined;
    for (const el of Array.from(document.querySelectorAll<HTMLElement>('.u2-dialog'))) {
      if (el.isConnected && !el.hidden &&
        (top === undefined || Number(el.style.zIndex || 0) >= Number(top.style.zIndex || 0)))
        top = el;
    }
    return top;
  }

  /** Whether an Enter means the dialog's default action: an unhandled one from a plain text input
   * or the dialog body, not one on or inside a control that acts on Enter itself. */
  static entersDefault(e: KeyboardEvent): boolean {
    if (e.key !== 'Enter' || e.defaultPrevented)
      return false;
    const target = e.target;
    return !(target instanceof Element) || target.closest(ENTER_OWNERS) === null;
  }

  /** A handler run the way its registrant would run it: with the scope ambient where it was
   * registered (none once that scope is gone, so a dialog that outlived its opener still answers),
   * in one batch. A throw reaches the caller. */
  private static _owned<A extends unknown[], R>(fn: (...args: A) => R): (...args: A) => R {
    const owner = Scope.ambient;
    return (...args: A) => batch(() => Scope.runWith(owner?.isDisposed ? undefined : owner, () => fn(...args)));
  }

  private static _titleButton(cls: string, part: string, glyph: string, label: string): HTMLButtonElement {
    const b = document.createElement('button');
    b.type = 'button';
    b.className = cls;
    b.dataset.u2Part = part;
    b.setAttribute('aria-label', label);
    b.append(icon(glyph));
    return b;
  }

  private static _focusable(el: HTMLElement): HTMLElement[] {
    const els = Array.from(el.querySelectorAll(FOCUSABLE)) as HTMLElement[];
    return els.filter((x) => x.offsetParent !== null);
  }
}
