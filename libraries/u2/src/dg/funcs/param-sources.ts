/* The per-param machinery behind FuncCallForm's dynamic routes (W2): an AsyncValue over
   `call.evalParamChoices` re-triggered through a version signal the dependency subscriptions
   bump, plus the on-field state element the loading/error surface renders into. Internal to
   `src/dg/funcs` — func-form.ts is the only importer; not exported through `src/dg/index.ts`. */
import {signal, computed, ReadonlySignal} from '../../core/signals.js';
import {Scope} from '../../core/scope.js';
import {AsyncValue, AsyncValueState} from '../../core/async-value.js';
import {button, span} from '../../core/elements.js';
import type {FuncCallLike, FuncCallParamLike} from './func-form.js';

/** What `evalParamChoices` answers (js-api functions.ts:377-392). */
export interface ChoicesResult {
  items: string[];
  values: Record<string, any>;
  lookup: Record<string, Record<string, any>> | null;
  dependsOn: string[];
}

export type ParamStateSlot = {kind: 'idle'} | {kind: 'loading', refresh?: boolean} |
  {kind: 'error', message: string, retry: () => void};

/** The state element riding a field's input box (`span.u2-param-source`): slot 0 is the choices
 * source, slot 1 the computed-default eval; an error wins over a loading, loading over idle.
 * A `refresh` loading — landed content re-evaluating — never disables the input, and the spinner
 * waits out {@link graceMs} so a fast eval never blinks; an error shows immediately.
 * Owned by the ambient (form) scope. */
export class ParamState {
  static graceMs = 150;
  readonly root = span('', 'u2-param-source');
  /** True while a slot loads a field that has nothing yet — what disables the input alongside
   * the orphan flag. */
  readonly busy: ReadonlySignal<boolean>;
  private readonly _slots = [signal<ParamStateSlot>({kind: 'idle'}),
    signal<ParamStateSlot>({kind: 'idle'})];
  private readonly _notice = signal<string | null>(null);
  private readonly _grace = signal(false);
  private _timer: ReturnType<typeof setTimeout> | undefined;

  constructor() {
    this.busy = computed(() =>
      this._slots.some((s) => s.value.kind === 'loading' && s.value.refresh !== true));
    const loading = computed(() => this._slots.some((s) => s.value.kind === 'loading'));
    Scope.ambient!.effect(() => {
      if (!loading.value) {
        clearTimeout(this._timer);
        this._timer = undefined;
        if (this._grace.peek())
          this._grace.value = false;
      }
      else if (this._timer === undefined && !this._grace.peek())
        this._timer = setTimeout(() => this._grace.value = true, ParamState.graceMs);
    });
    Scope.ambient!.effect(() => this._render());
    Scope.ambient!.own(() => clearTimeout(this._timer));
  }

  set(slot: number, state: ParamStateSlot): void {
    if (state.kind === 'idle' && this._slots[slot].peek().kind === 'idle')
      return;
    this._slots[slot].value = state;
  }

  /** A transient informational message — shown while nothing louder (an error, a grace-past
   * loading) occupies the element; null clears it. */
  notice(message: string | null): void {
    this._notice.value = message;
  }

  private _render(): void {
    const states = this._slots.map((s) => s.value);
    const notice = this._notice.value;
    const shown = states.find((s) => s.kind === 'error') ??
      (this._grace.value ? states.find((s) => s.kind === 'loading') : undefined);
    const root = this.root;
    root.textContent = '';
    root.className = 'u2-param-source';
    root.hidden = shown === undefined && notice === null;
    if (root.hidden)
      return;
    if (shown === undefined) {
      root.classList.add('u2-param-source-notice');
      const message = span(notice!, 'u2-param-source-message');
      message.title = notice!;
      root.append(message);
      return;
    }
    if (shown.kind === 'loading') {
      root.classList.add('u2-param-source-loading');
      root.append(span('', 'u2-loader-spinner'));
      return;
    }
    root.classList.add('u2-param-source-error');
    const message = span(shown.message, 'u2-param-source-message');
    message.title = shown.message;
    root.append(message);
    const retry = button('Retry', shown.retry);
    retry.classList.add('u2-param-source-retry');
    root.append(retry);
  }
}

/** One dynamic-choices source: `evalParamChoices` behind an {@link AsyncValue} whose `deps` read
 * a version signal every dependency subscription bumps, so a dep edit re-runs the eval through
 * the 200 ms debounce, while {@link start} runs the initial one immediately. Deps ride
 * `param.onChanged` — the param stream, never field signals or captions — so a dep with no field
 * (an unsupported or fieldless param) still triggers. Constructed under the form's scope;
 * {@link dispose} on a `source` rebind. */
export class ParamSource {
  private readonly _version = signal(0);
  private readonly _scope = new Scope();
  private readonly _call: FuncCallLike;
  private readonly _name: string;
  private readonly _apply: (r: ChoicesResult) => void;
  private readonly _state: ParamState;
  private readonly _disown: () => void;
  private _av!: AsyncValue<ChoicesResult>;
  private _depNames: string[] = [];
  private _depSubs: (() => void)[] = [];
  private _landed = false;

  constructor(call: FuncCallLike, name: string, apply: (r: ChoicesResult) => void,
    state: ParamState) {
    this._call = call;
    this._name = name;
    this._apply = apply;
    this._state = state;
    const owner = Scope.ambient!;
    const cleanup = () => this._dispose();
    owner.own(cleanup);
    this._disown = () => owner.disown(cleanup);
    Scope.runWith(this._scope, () => {
      this._av = new AsyncValue(() => this._call.evalParamChoices!(this._name),
        {deps: () => { this._version.value; }, debounceMs: 200});
      this._scope.effect(() => this._onState(this._av.state.value));
    });
  }

  /** The initial, immediate run; the promise lands with the first ready-or-error state — what
   * {@link FuncCallForm.settled} aggregates. */
  start(): Promise<void> {
    const landed = settle(this._av, this._scope);
    void this._av.refresh();
    return landed;
  }

  retry(): void {
    void this._av.refresh();
  }

  dispose(): void {
    this._disown();
    this._dispose();
  }

  private _dispose(): void {
    this._unsubscribe();
    this._scope.dispose();
  }

  private _onState(state: AsyncValueState<ChoicesResult>): void {
    // once anything has landed, a refresh keeps the input enabled — only the initial load,
    // over an empty select, disables (the Dart form never disables at all)
    if (state.kind === 'loading')
      this._state.set(0, {kind: 'loading', refresh: this._landed});
    else if (state.kind === 'ready') {
      this._landed = true;
      this._state.set(0, {kind: 'idle'});
      this._reconcile(state.value.dependsOn);
      this._apply(state.value);
    }
    // on error the stale items stay usable and the previous dependency list stays armed
    // (the Dart `finally`, fpe:628-637); a failed INITIAL load leaves nothing usable, so
    // `_landed` moves on ready only and later runs still count as initial
    else if (state.kind === 'error')
      this._state.set(0, {kind: 'error', message: 'Couldn\'t load choices: ' + state.message,
        retry: () => this.retry()});
  }

  /** Same dep names → do nothing: reconciliation must never bump {@link _version} or
   * resubscribe, or it becomes the echo-storm engine. */
  private _reconcile(names: string[]): void {
    if (names.length === this._depNames.length && names.every((n, i) => n === this._depNames[i]))
      return;
    this._unsubscribe();
    this._depNames = names.slice();
    const byName = new Map<string, FuncCallParamLike>();
    for (const param of this._call.inputParams.values())
      byName.set(param.name, param);
    for (const name of names) {
      if (name === this._name)
        continue;
      const param = byName.get(name);
      if (param === undefined)
        continue;
      const sub = param.onChanged.subscribe(() => this._version.value++);
      this._depSubs.push(() => sub.unsubscribe());
    }
  }

  private _unsubscribe(): void {
    for (const off of this._depSubs)
      off();
    this._depSubs = [];
  }
}

/** The first landing — ready or error — read off the state signal: `refresh()`'s own promise
 * also resolves for a run a dependency change aborted, which has not landed. A dispose before
 * the first landing resolves too — a form torn down or rebound mid-flight never hangs `settled`. */
function settle<T>(av: AsyncValue<T>, owner: Scope): Promise<void> {
  return new Promise((resolve) => {
    const scope = new Scope();
    const done = () => {
      resolve();
      scope.dispose();
    };
    owner.own(done);
    scope.effect(() => {
      const kind = av.state.value.kind;
      if (kind !== 'ready' && kind !== 'error')
        return;
      owner.disown(done);
      // deferred: disposing the scope from inside its own effect is what the microtask avoids
      queueMicrotask(done);
    });
  });
}
