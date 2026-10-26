/* Multi-step wizard. Panels are hidden, never detached, so a step's state survives navigation;
   lazy content is built once under the step's own scope. In dialog mode the wizard's own footer
   carries the buttons — the Dialog's button row stays empty — and any close that is not a finish
   reports a cancel, which covers CANCEL, the ✕ and Esc in one place; while a step's commit, an
   action or the finish runs, none of them closes and the rail does not move: Esc and ✕ abort the
   run's signal instead. A run that throws or rejects keeps the step and shows the message. */
import {signal, computed, untracked, batch, ReadonlySignal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {Scope} from '../../core/scope.js';
import {button} from '../../core/elements.js';
import {Dialog} from './dialog.js';
import {statusLine} from '../display/badge.js';
import type {IItemsState} from '../../core/widget-like.js';

/** What a step's commit, an action and the finish get: `signal` aborts when the user presses Esc
 * or ✕ while the work runs. */
export interface WizardRunContext {
  signal: AbortSignal;
}

export interface WizardStep {
  id: string;
  title: string;
  /** A function is lazy: built on first visit. */
  content: HTMLElement | (() => HTMLElement);
  /** Gates NEXT/FINISH: `false` blocks silently, a string blocks and shows the reason next to the
   * buttons. A `ReadonlySignal<string | null>` (null = ok) composes directly with `Form.validity`. */
  canProceed?: ReadonlySignal<boolean> | ReadonlySignal<string | null> | (() => boolean | string | null);
  /** A step reached only once the work is done — a report over what was written. There is
   * nothing to go back to and nothing to cancel, so only the closing button stands, reading
   * CLOSE; a plain Enter does not press it. */
  done?: boolean;
  /** What the closing button reads on the last step (FINISH by default). */
  finishText?: string;
  /** What NEXT reads on this step (NEXT by default) — CREATE on a step whose commit registers; a
   * signal where the wording follows the state (SAVE, then REOPEN once a save's outcome is unknown).
   * A named NEXT is raised, as is one with a {@link commit} and FINISH; a plain NEXT is flat. */
  nextText?: string | ReadonlySignal<string>;
  /** Work NEXT runs before leaving this step: the wizard waits on it and nothing closes meanwhile;
   * `false` keeps the step (a commit the server refused), and so does a throw or a rejection, whose
   * message the footer shows. */
  commit?: (ctx: WizardRunContext) => void | boolean | Promise<void | boolean>;
  /** Buttons of this step's own, before NEXT in the footer — a check the user runs before
   * finishing. While one runs, the footer waits. */
  actions?: WizardAction[];
  /** Runs on every visit. What it builds lives until the step's next visit or the wizard's end. */
  onActivate?: (step: WizardStep) => void;
}

export interface WizardAction {
  text: string;
  run: (ctx: WizardRunContext) => void | Promise<void>;
  /** Follows the signal; enabled otherwise. */
  enabled?: ReadonlySignal<boolean>;
  /** Follows the signal; shown on its step otherwise. */
  visible?: ReadonlySignal<boolean>;
}

export interface WizardOptions {
  steps: WizardStep[];
  /** The step to open on (the first by default); the ones before it count as visited — the host
   * pre-filled them — so BACK and the rail reach them. */
  start?: string;
  /** Answering `false` keeps the wizard open — a finish the server refused; the footer waits
   * on a promise. Also the `finish` event (`wizard.on('finish', …)`), fired once it went through. */
  onFinish?: (ctx: WizardRunContext) => void | boolean | Promise<void | boolean>;
  /** Also the `cancel` event. */
  onCancel?: () => void;
}

interface Step {
  options: WizardStep;
  marker: HTMLElement;
  circle: HTMLElement;
  panel: HTMLElement;
  scope: Scope;
  /** What the last {@link WizardStep.onActivate} built. */
  visit?: Scope;
  actions: HTMLButtonElement[];
  builder?: () => HTMLElement;
}

let wizardCount = 0;

/** Events: `step` (the step id, on every visit), `finish`, `cancel`. */
export class Wizard extends Control {
  readonly currentStep: ReadonlySignal<string>;
  readonly completed: ReadonlySignal<boolean>;
  /** A commit, an action or the finish is running. */
  readonly busy: ReadonlySignal<boolean>;
  /** The message of the last run that threw or rejected; null when there is none. */
  readonly error: ReadonlySignal<string | null>;
  /** The footer's status line, left of the gate's reason — the host writes what it is doing. */
  readonly status = document.createElement('span');

  private readonly _idPrefix = `u2-wizard-${++wizardCount}`;
  private readonly _steps: Step[] = [];
  private readonly _visited = new Set<string>();
  private readonly _options: WizardOptions;
  /** Where the options' callbacks run: the scope that built the wizard. */
  private readonly _owner: Scope | undefined;
  private readonly _rail = document.createElement('ol');
  private readonly _content = document.createElement('div');
  private readonly _footer = document.createElement('div');
  private readonly _reason: HTMLElement;
  private readonly _errorRow = document.createElement('span');
  private readonly _back: HTMLButtonElement;
  private readonly _next: HTMLButtonElement;
  private _cancel: HTMLButtonElement | undefined;
  private readonly _index = signal(0);
  private readonly _completed = signal(false);
  private readonly _recheck = signal(0);
  private readonly _busy = signal(false);
  private readonly _error = signal<string | null>(null);
  private readonly _failed = signal<ReadonlySet<string>>(new Set());
  /** The run in flight: Esc and ✕ abort it. */
  private _pending: AbortController | undefined;
  /** The last step whose commit went through: what it registered stands, so nothing goes back to it. */
  private _committed = -1;
  private readonly _blockReason: ReadonlySignal<string | null>;
  private _dialog: Dialog | undefined;

  constructor(options: WizardOptions) {
    super();
    this._options = options;
    this._owner = Scope.ambient;
    this.currentStep = computed(() => this._steps[this._index.value].options.id);
    this.completed = this._completed;
    this.busy = this._busy;
    this.error = this._error;
    this._blockReason = computed(() => {
      this._recheck.value;
      return this._busy.value ? '' : Wizard._gate(this._steps[this._index.value].options.canProceed);
    });

    this.root.classList.add('u2-wizard');
    this.root.dataset.u2 = 'wizard';
    this._rail.className = 'u2-wizard-steps';
    this._rail.dataset.u2Part = 'steps';
    this._content.className = 'u2-wizard-content';
    this._content.dataset.u2Part = 'body';
    this._footer.className = 'u2-wizard-buttons';
    this._footer.dataset.u2Part = 'footer';
    this._errorRow.className = 'u2-wizard-error';
    this._errorRow.dataset.u2Part = 'error-row';
    this._errorRow.setAttribute('role', 'alert');
    this._errorRow.hidden = true;
    this.status.className = 'u2-wizard-status';
    this.status.dataset.u2Part = 'status';
    const left = document.createElement('span');
    left.className = 'u2-wizard-left';
    left.append(this.status, this._errorRow);

    this._back = this._button('BACK', () => this.back());
    this._next = this._button('NEXT', () => this.next());
    this._footer.append(left, this._back, this._next);
    this.root.append(this._rail, this._content, this._footer);

    const scope = this.scope;
    scope.listen(this._rail, 'click', (e) => this._onRailClick(e));
    scope.listen(this._rail, 'keydown', (e) => this._onRailKeyDown(e as KeyboardEvent));
    scope.listen(this._content, 'keydown', (e) => this._onContentKeyDown(e as KeyboardEvent));

    for (const step of options.steps)
      this._addStep(step);
    const start = options.start === undefined ? 0 : options.steps.findIndex((s) => s.id === options.start);
    if (start < 0)
      throw new Error(`Wizard: no step "${options.start}" to start on`);
    for (const step of options.steps.slice(0, start))
      this._visited.add(step.id);
    this._index.value = start;
    this._reason = this.runInScope(() => statusLine(computed(() => this._blockReason.value || null)));
    this._reason.classList.add('u2-wizard-reason', 'u2-status-line-warning');
    this._reason.dataset.u2Part = 'reason';
    this._errorRow.parentElement!.insertBefore(this._reason, this._errorRow);

    this.effect(() => this._applyStep());
    this.effect(() => this._applyNextText());
    this.effect(() => this._applyGate());
    this.effect(() => {
      const message = this._error.value;
      this._errorRow.textContent = message ?? '';
      this._errorRow.hidden = message === null;
    });
  }

  next(): void {
    this._recheck.value = this._recheck.peek() + 1;
    if (this._blockReason.peek() !== null)
      return;
    const index = this._index.peek();
    if (index === this._steps.length - 1) {
      this._finish();
      return;
    }
    const commit = this._steps[index].options.commit;
    const advance = (ok: unknown): void => {
      if (ok === false || this.scope.isDisposed || this._index.peek() !== index)
        return;
      if (commit)
        this._committed = index;
      this._index.value = index + 1;
    };
    if (commit)
      this._run(commit, advance);
    else
      advance(undefined);
  }

  back(): void {
    const index = this._index.peek();
    if (index - 1 > this._committed && !this._busy.peek())
      this._index.value = index - 1;
  }

  /** Back to a visited step after the last committed one; forward to a visited one only through the
   * gate of every step on the way and past no committing step, or to the next one through {@link next}.
   * Nothing moves while the wizard is busy. */
  goTo(id: string): void {
    const target = this._steps.findIndex((s) => s.options.id === id);
    const index = this._index.peek();
    if (target < 0 || target === index || target <= this._committed || this._busy.peek())
      return;
    if (this._visited.has(id) && (target < index ||
        this._steps.slice(index, target).every((s) => !s.options.commit && Wizard._gate(s.options.canProceed) === null)))
      this._index.value = target;
    else if (target === index + 1)
      this.next();
  }

  /** Renames the step at `index` in place — a bound step title following its source. */
  setStepTitle(index: number, title: string): void {
    const step = this._steps[index];
    if (step === undefined)
      return;
    step.options = {...step.options, title};
    step.marker.querySelector('.u2-wizard-title')!.textContent = title;
  }

  /** Picks up where `from` was — its visited steps and its open one, matched by id: what a rebuild of
   * the same steps (a designer patch) carries over. */
  resumeFrom(from: Wizard): void {
    for (const id of from._visited) {
      if (this._steps.some((s) => s.options.id === id))
        this._visited.add(id);
    }
    const at = this._steps.findIndex((s) => s.options.id === from.currentStep.peek());
    if (at >= 0)
      this._index.value = at;
  }

  /** Marks a step as failed: its marker reads "!" in the error colour until cleared. */
  setStepFailed(id: string, failed: boolean): void {
    const next = new Set(this._failed.peek());
    if (failed)
      next.add(id);
    else
      next.delete(id);
    this._failed.value = next;
  }

  /** Shows the wizard as a modal dialog; repeated calls reopen the same one. */
  openInDialog(title: string, options: {width?: number, height?: number} = {}): Dialog {
    if (!this._dialog) {
      const dialog = this.runInScope(() => Dialog.create(title, {keepAlive: true}).add(this)
        .closeGuard(() => this._guard())
        // ✕ and Esc on a done step are the CLOSE button; anywhere else they are a cancel
        .onClose(() => {
          if (this._completed.peek())
            return;
          if (this._steps[this._index.peek()].options.done === true)
            this._finish();
          else
            this._cancelled();
        }));
      this._dialog = dialog;
      this._cancel = this._button('CANCEL', () => {
        if (!this._busy.peek())
          dialog.close();
      });
      this._cancel.classList.add('u2-wizard-cancel');
      this._footer.insertBefore(this._cancel, this._back);
    }
    return this._dialog.show({modal: true, width: options.width, height: options.height});
  }

  /** The base state plus the current step, whether a run is in flight, the gate and the error row,
   * and the steps as items. */
  getState(): ReturnType<Control['getState']> {
    const reason = this._blockReason.peek();
    const error = this._error.peek();
    return {...super.getState(), current: this.currentStep.peek(), busy: this._busy.peek(),
      valid: reason === null && error === null, error: error ?? (reason || null)};
  }

  getItems(): IItemsState {
    const index = this._index.peek();
    const list = this._steps.map((step, i) => ({key: step.options.id, label: step.options.title, index: i,
      current: i === index, completed: i < index && this._visited.has(step.options.id),
      disabled: step.marker.getAttribute('aria-disabled') === 'true', rendered: true}));
    return {count: list.length, list, truncated: false};
  }

  private _addStep(step: WizardStep): void {
    const index = this._steps.length;
    const marker = document.createElement('li');
    marker.className = 'u2-wizard-step';
    marker.id = `${this._idPrefix}-step-${index}`;
    marker.dataset.id = step.id;
    marker.dataset.u2Item = step.id;
    marker.tabIndex = index === 0 ? 0 : -1;

    const circle = document.createElement('span');
    circle.className = 'u2-wizard-marker';
    circle.textContent = String(index + 1);
    const title = document.createElement('span');
    title.className = 'u2-wizard-title';
    title.dataset.u2Part = 'item-label';
    title.textContent = step.title;
    marker.append(circle, title);
    this._rail.append(marker);

    const panel = document.createElement('div');
    panel.className = 'u2-wizard-panel';
    panel.style.display = 'none';
    panel.setAttribute('role', 'group');
    panel.setAttribute('aria-labelledby', marker.id);
    this._content.append(panel);

    const scope = this.scope.child('wizard-step');
    const lazy = typeof step.content === 'function' ? step.content : undefined;
    if (!lazy)
      panel.append(step.content as HTMLElement);
    const actions = (step.actions ?? []).map((action) => {
      const b = this._button(action.text, () => this._run(action.run));
      this.effect(() => b.style.display =
        this._index.value === index && action.visible?.value !== false ? '' : 'none');
      this._footer.insertBefore(b, this._next);
      return b;
    });
    this._steps.push({options: step, marker, circle, panel, scope, actions, builder: lazy});
  }

  /** Runs a commit, an action or the finish as its registrant would, and hands a plain answer to
   * `done` at once — a second NEXT in the same task finds the step already gone. A thenable makes
   * the footer wait: NEXT and the actions are disabled and the gate says nothing until it settles.
   * A throw or a rejection keeps everything where it is and shows the message (an abort, nothing). */
  private _run(work: (ctx: WizardRunContext) => unknown, done?: (ok: unknown) => void): void {
    if (this._busy.peek())
      return;
    this._error.value = null;
    const abort = new AbortController();
    let answer: unknown;
    try {
      answer = this._call(() => work({signal: abort.signal}));
    } catch (e) {
      this._fail(e);
      return;
    }
    if (typeof (answer as PromiseLike<unknown> | null)?.then !== 'function') {
      done?.(answer);
      return;
    }
    this._pending = abort;
    this._busy.value = true;
    const settle = (fn: () => void) => {
      if (this._pending !== abort || this.scope.isDisposed)
        return;
      this._pending = undefined;
      batch(() => {
        this._busy.value = false;
        fn();
      });
    };
    (answer as PromiseLike<unknown>).then((ok) => settle(() => done?.(ok)), (e) => settle(() => this._fail(e)));
  }

  /** The callbacks of the options run with the scope that built the wizard ambient (none once it
   * is gone), in one batch. A throw reaches the caller. */
  private _call<T>(fn: () => T): T {
    const owner = this._owner;
    return batch(() => Scope.runWith(owner?.isDisposed ? undefined : owner, fn));
  }

  private _fail(e: unknown): void {
    if ((e as {name?: unknown} | null)?.name !== 'AbortError')
      this._error.value = e instanceof Error ? e.message : String(e);
  }

  /** The dialog's veto: nothing closes while a run is in flight — the Esc or ✕ aborts its signal. */
  private _guard(): boolean {
    this._pending?.abort();
    return !this._busy.peek();
  }

  private _applyStep(): void {
    const index = this._index.value;
    const current = this._steps[index];
    const failed = this._failed.value;
    this._visited.add(current.options.id);
    for (let i = 0; i < this._steps.length; i++) {
      const step = this._steps[i];
      const on = i === index;
      const visited = !on && this._visited.has(step.options.id);
      // a step the user went BACK past stays reachable, but it is ahead of them: only the steps
      // behind the current one read as completed
      const done = visited && i < index;
      step.marker.classList.toggle('u2-wizard-step-current', on);
      step.marker.classList.toggle('u2-wizard-step-done', done);
      step.marker.classList.toggle('u2-wizard-step-failed', failed.has(step.options.id));
      step.marker.setAttribute('aria-disabled', String(!on && (!visited || i <= this._committed)));
      step.marker.tabIndex = on ? 0 : -1;
      step.circle.textContent = failed.has(step.options.id) ? '!' : done ? '✓' : String(i + 1);
      step.panel.style.display = on ? '' : 'none';
      if (on)
        step.marker.setAttribute('aria-current', 'step');
      else
        step.marker.removeAttribute('aria-current');
    }
    const done = current.options.done === true;
    // NEXT reads CLOSE now, so a second Enter must not land on it; focus stays inside the dialog (Esc,
    // the Tab trap) — Chrome drops a button disabled while its commit ran onto the body
    const active = document.activeElement;
    if (done && (active === this._next || active === null || active === document.body))
      current.marker.focus();
    this._back.style.display = index === 0 || done ? 'none' : '';
    if (this._cancel !== undefined)
      this._cancel.style.display = done ? 'none' : '';

    untracked(() => {
      this._error.value = null;
      const builder = current.builder;
      if (builder) {
        current.builder = undefined;
        current.panel.append(Scope.runWith(current.scope, builder));
      }
      const onActivate = current.options.onActivate;
      if (onActivate) {
        current.visit?.dispose();
        const visit = current.visit = current.scope.child('wizard-visit');
        visit.run(() => onActivate(current.options));
      }
      this.fireEvent('step', current.options.id);
    });
  }

  private _applyNextText(): void {
    const index = this._index.value;
    const current = this._steps[index].options;
    const next = typeof current.nextText === 'string' ? current.nextText : current.nextText?.value;
    const last = index === this._steps.length - 1;
    this._next.textContent = !last ? next ?? 'NEXT' : current.done === true ? 'CLOSE' : current.finishText ?? 'FINISH';
    // raised only where pressing it commits
    this._next.classList.toggle('u2-btn-primary', current.done !== true &&
      (last || current.commit !== undefined || current.nextText !== undefined));
  }

  private _applyGate(): void {
    const reason = this._blockReason.value;
    const busy = this._busy.value;
    const current = this._steps[this._index.value];
    this._next.disabled = reason !== null;
    this._next.title = reason ?? '';
    this._back.disabled = busy;
    if (this._cancel !== undefined)
      this._cancel.disabled = busy;
    this.root.classList.toggle('u2-wizard-busy', busy);
    this.root.setAttribute('aria-busy', String(busy));
    current.actions.forEach((b, i) => b.disabled = busy || current.options.actions![i].enabled?.value === false);
  }

  private _finish(): void {
    const onFinish = this._options.onFinish;
    if (onFinish === undefined)
      this._complete();
    else
      this._run(onFinish, (ok) => ok !== false && this._complete());
  }

  private _complete(): void {
    this._completed.value = true;
    this.fireEvent('finish');
    if (this._dialog)
      this._dialog.close();
  }

  private _cancelled(): void {
    try {
      this._call(() => this._options.onCancel?.());
    } catch (e) {
      this.scope.reportError(e);
    }
    this.fireEvent('cancel');
  }

  private _button(text: string, onClick: () => void, primary?: boolean): HTMLButtonElement {
    return this.runInScope(() => button(text, () => onClick(), {primary}));
  }

  private _focusStep(index: number): void {
    for (let i = 0; i < this._steps.length; i++)
      this._steps[i].marker.tabIndex = i === index ? 0 : -1;
    this._steps[index].marker.focus();
  }

  private _onRailClick(e: Event): void {
    const marker = (e.target as HTMLElement).closest('.u2-wizard-step') as HTMLElement | null;
    if (marker)
      this.goTo(marker.dataset.id!);
  }

  private _onRailKeyDown(e: KeyboardEvent): void {
    const current = this._steps.findIndex((s) => s.marker === document.activeElement);
    if (current < 0)
      return;
    const last = this._steps.length - 1;
    let next: number;
    switch (e.key) {
      case 'ArrowLeft': next = current === 0 ? last : current - 1; break;
      case 'ArrowRight': next = current === last ? 0 : current + 1; break;
      case 'Enter':
      case ' ':
        e.preventDefault();
        this.goTo(this._steps[current].options.id);
        return;
      default: return;
    }
    e.preventDefault();
    this._focusStep(next);
  }

  private _onContentKeyDown(e: KeyboardEvent): void {
    if (e.key !== 'Enter' || e.defaultPrevented)
      return;
    const current = this._steps[this._index.peek()];
    const last = current === this._steps[this._steps.length - 1];
    if (e.ctrlKey || e.metaKey) {
      if (!last)
        return;
    } else if (!Dialog.entersDefault(e) || current.options.done === true)
      return;
    e.preventDefault();
    this.next();
  }

  private static _gate(gate: WizardStep['canProceed']): string | null {
    if (!gate)
      return null;
    const value = typeof gate === 'function' ? gate() : gate.value;
    if (value === true || value === null)
      return null;
    return value === false ? '' : value;
  }
}
