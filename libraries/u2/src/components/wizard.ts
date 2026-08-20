/* Multi-step wizard. Panels are hidden, never detached, so a step's state survives navigation;
   lazy content is built once under the step's own scope. In dialog mode the wizard's own footer
   carries the buttons — the Dialog's button row stays empty — and any close that is not a finish
   reports a cancel, which covers CANCEL, the ✕ and Esc in one place. */
import {signal, computed, untracked, ReadonlySignal} from '../core/signals.js';
import {Control} from '../core/component.js';
import {Scope} from '../core/scope.js';
import {button} from '../core/elements.js';
import {Dialog} from './dialog.js';

export interface WizardStep {
  id: string;
  title: string;
  /** A function is lazy: built on first visit. */
  content: HTMLElement | (() => HTMLElement);
  /** Gates NEXT/FINISH: `false` blocks silently, a string blocks and shows the reason next to the
   * buttons. A `ReadonlySignal<string | null>` (null = ok) composes directly with `Form.validity`. */
  canProceed?: ReadonlySignal<boolean> | ReadonlySignal<string | null> | (() => boolean | string | null);
  onActivate?: (step: WizardStep) => void;
}

export interface WizardOptions {
  steps: WizardStep[];
  onFinish?: () => void;
  onCancel?: () => void;
}

interface Step {
  options: WizardStep;
  marker: HTMLElement;
  circle: HTMLElement;
  panel: HTMLElement;
  scope: Scope;
  builder?: () => HTMLElement;
}

let wizardCount = 0;

export class Wizard extends Control {
  readonly currentStep: ReadonlySignal<string>;
  readonly completed: ReadonlySignal<boolean>;

  private readonly _idPrefix = `u2-wizard-${++wizardCount}`;
  private readonly _steps: Step[] = [];
  private readonly _visited = new Set<string>();
  private readonly _options: WizardOptions;
  private readonly _rail = document.createElement('ol');
  private readonly _content = document.createElement('div');
  private readonly _footer = document.createElement('div');
  private readonly _reason = document.createElement('span');
  private readonly _back: HTMLButtonElement;
  private readonly _next: HTMLButtonElement;
  private readonly _index = signal(0);
  private readonly _completed = signal(false);
  private readonly _recheck = signal(0);
  private readonly _blockReason: ReadonlySignal<string | null>;
  private _dialog: Dialog | undefined;

  constructor(options: WizardOptions) {
    super();
    this._options = options;
    this.currentStep = computed(() => this._steps[this._index.value].options.id);
    this.completed = this._completed;
    this._blockReason = computed(() => {
      this._recheck.value;
      return Wizard._gate(this._steps[this._index.value].options.canProceed);
    });

    this.root.classList.add('u2-wizard');
    this.root.dataset.u2 = 'wizard';
    this._rail.className = 'u2-wizard-steps';
    this._content.className = 'u2-wizard-content';
    this._footer.className = 'u2-wizard-buttons';
    this._reason.className = 'u2-wizard-reason';

    this._back = this._button('BACK', () => this.back());
    this._next = this._button('NEXT', () => this.next(), true);
    this._footer.append(this._reason, this._back, this._next);
    this.root.append(this._rail, this._content, this._footer);

    this._listen(this._rail, 'click', (e) => this._onRailClick(e));
    this._listen(this._rail, 'keydown', (e) => this._onRailKeyDown(e as KeyboardEvent));
    this._listen(this._content, 'keydown', (e) => this._onContentKeyDown(e as KeyboardEvent));

    for (const step of options.steps)
      this._addStep(step);

    this.effect(() => this._applyStep());
    this.effect(() => this._applyGate());
  }

  next(): void {
    this._recheck.value = this._recheck.peek() + 1;
    if (this._blockReason.peek() !== null)
      return;
    const index = this._index.peek();
    if (index < this._steps.length - 1)
      this._index.value = index + 1;
    else
      this._finish();
  }

  back(): void {
    const index = this._index.peek();
    if (index > 0)
      this._index.value = index - 1;
  }

  goTo(id: string): void {
    const target = this._steps.findIndex((s) => s.options.id === id);
    const index = this._index.peek();
    if (target < 0 || target === index)
      return;
    if (this._visited.has(id))
      this._index.value = target;
    else if (target === index + 1)
      this.next();
  }

  /** Shows the wizard as a modal dialog; repeated calls reopen the same one. */
  openInDialog(title: string, options: {width?: number, height?: number} = {}): Dialog {
    if (!this._dialog) {
      const dialog = this.run(() => Dialog.create(title).add(this));
      this._dialog = dialog;
      this._footer.insertBefore(this._button('CANCEL', () => dialog.close()), this._back);
      let open = false;
      this.effect(() => {
        const closed = open && !dialog.isOpen.value;
        open = dialog.isOpen.value;
        const onCancel = this._options.onCancel;
        if (closed && onCancel && !this._completed.peek())
          untracked(() => onCancel());
      });
    }
    return this._dialog.show({modal: true, width: options.width, height: options.height});
  }

  private _addStep(step: WizardStep): void {
    const index = this._steps.length;
    const marker = document.createElement('li');
    marker.className = 'u2-wizard-step';
    marker.id = `${this._idPrefix}-step-${index}`;
    marker.dataset.id = step.id;
    marker.tabIndex = index === 0 ? 0 : -1;

    const circle = document.createElement('span');
    circle.className = 'u2-wizard-marker';
    circle.textContent = String(index + 1);
    const title = document.createElement('span');
    title.className = 'u2-wizard-title';
    title.textContent = step.title;
    marker.append(circle, title);
    this._rail.append(marker);

    const panel = document.createElement('div');
    panel.className = 'u2-wizard-panel';
    panel.style.display = 'none';
    panel.setAttribute('role', 'group');
    panel.setAttribute('aria-labelledby', marker.id);
    this._content.append(panel);

    const scope = new Scope();
    this.own(() => scope.dispose());
    const lazy = typeof step.content === 'function' ? step.content : undefined;
    if (!lazy)
      panel.append(step.content as HTMLElement);
    this._steps.push({options: step, marker, circle, panel, scope, builder: lazy});
  }

  private _applyStep(): void {
    const index = this._index.value;
    const current = this._steps[index];
    this._visited.add(current.options.id);
    for (let i = 0; i < this._steps.length; i++) {
      const step = this._steps[i];
      const on = i === index;
      const done = !on && this._visited.has(step.options.id);
      step.marker.classList.toggle('u2-wizard-step-current', on);
      step.marker.classList.toggle('u2-wizard-step-done', done);
      step.marker.setAttribute('aria-disabled', String(!on && !done));
      step.marker.tabIndex = on ? 0 : -1;
      step.circle.textContent = done ? '✓' : String(i + 1);
      step.panel.style.display = on ? '' : 'none';
      if (on)
        step.marker.setAttribute('aria-current', 'step');
      else
        step.marker.removeAttribute('aria-current');
    }
    this._back.style.display = index === 0 ? 'none' : '';
    this._next.textContent = index === this._steps.length - 1 ? 'FINISH' : 'NEXT';

    const builder = current.builder;
    if (builder) {
      current.builder = undefined;
      untracked(() => current.panel.append(Scope.runWith(current.scope, builder)));
    }
    const onActivate = current.options.onActivate;
    if (onActivate)
      untracked(() => onActivate(current.options));
  }

  private _applyGate(): void {
    const reason = this._blockReason.value;
    this._next.disabled = reason !== null;
    this._reason.textContent = reason ?? '';
  }

  private _finish(): void {
    this._completed.value = true;
    if (this._options.onFinish)
      this._options.onFinish();
    if (this._dialog)
      this._dialog.close();
  }

  private _button(text: string, onClick: () => void, primary?: boolean): HTMLButtonElement {
    return this.run(() => button(text, onClick, {primary}));
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
    const last = this._index.peek() === this._steps.length - 1;
    if (e.ctrlKey || e.metaKey) {
      if (!last)
        return;
    } else if (e.target instanceof HTMLTextAreaElement)
      return;
    e.preventDefault();
    this.next();
  }

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
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
