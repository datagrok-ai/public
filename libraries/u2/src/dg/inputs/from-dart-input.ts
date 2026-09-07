/* Platform inputs inside u2 — the inbound half of the input bridge (core/docs/features/ui2/
   INTEROP.md). Structural on purpose: nothing here imports datagrok-api, so a real `DG.InputBase`
   and a headless fake are equally acceptable. */
import {Input, InputOptions} from '../../core/input-base.js';
import {signal, Signal} from '../../core/signals.js';

/** Everything the bridge reads off a platform input; `DG.InputBase` satisfies it structurally. */
export interface DartInputLike<T = unknown> {
  root: HTMLElement;
  value: T;
  caption: string;
  onChanged: {subscribe(fn: (value: T) => void): {unsubscribe(): void}};
  dart?: unknown;
}

interface PlatformInputOptions<T> extends InputOptions<T> {
  dgInput: DartInputLike<T>;
}

/** A `DG.InputBase` seen by u2 as an `Input<T>`: the platform root is the whole input, so caption
 * and validation message stay the platform's and nothing renders twice. The platform value is
 * canonical — signal writes reach it through `dgInput.value`, echo-suppressed by comparison — and
 * Dart-side validators (regex, func-param, script) reach `Form.validity` through {@link validity}.
 * The bridge never owns the wrapped input: `dispose()` severs the subscriptions and nothing else. */
export class PlatformInput<T> extends Input<T, PlatformInputOptions<T>> {
  private _messages!: Signal<string | null>;

  constructor(dgInput: DartInputLike<T>, name?: string) {
    super({dgInput, name: name ?? dgInput.caption}, dgInput.value);
    this.root.dataset.u2 = 'dart-input';
    // u2's editor skin (and the form's focus target) belongs on the platform's editor element, not
    // on its root, which is a whole input
    const editor = dgInput.root.querySelector<HTMLElement>('.ui-input-editor');
    if (editor) {
      dgInput.root.classList.remove('u2-input-editor');
      editor.classList.add('u2-input-editor');
    }

    const changed = dgInput.onChanged.subscribe((value) => this.value.value = value);
    this.own(() => changed.unsubscribe());
    this.effect(() => {
      const value = this.value.value;
      if (!Object.is(value, dgInput.value))
        dgInput.value = value;
    });

    const globals = globalThis as any;
    const messagesOf = globals.grok_InputBase_Get_ValidationMessages;
    const onValidated = globals.grok_InputBase_OnValidated;
    if (dgInput.dart == null || typeof messagesOf !== 'function' || typeof onValidated !== 'function')
      return;
    this._messages.value = PlatformInput._join(messagesOf(dgInput.dart));
    const validated = onValidated(dgInput.dart)
      .subscribe((messages: string[]) => this._messages.value = PlatformInput._join(messages));
    this.own(() => validated.unsubscribe());
    this.root.querySelector('.u2-input-error')?.remove();
  }

  get dgInput(): DartInputLike<T> {
    return this.options.dgInput;
  }

  // registered before the base builds its validation effect, so that effect tracks _messages and
  // re-runs the whole validation on every platform validation run
  protected createEditor(): HTMLElement {
    this._messages = signal<string | null>(null);
    this.addValidator(() => this._messages.value);
    return this.options.dgInput.root;
  }

  private static _join(messages: string[] | null): string | null {
    return messages != null && messages.length > 0 ? messages.join('\n') : null;
  }
}

/** Wraps a platform input as a u2 `Input`, keyed in forms by `name` (the platform caption when
 * omitted). An `asDartInput` product is unwrapped rather than wrapped twice (INTEROP.md): the u2
 * input it was built from comes back as is. */
export function fromDartInput<T>(dgInput: DartInputLike<T>, name?: string): Input<T> {
  // structural brand: importing DartInput would close a cycle, since input-bridge imports this file
  const wrapped = (dgInput as {u2?: unknown}).u2;
  return wrapped instanceof Input ? wrapped as Input<T> : new PlatformInput<T>(dgInput, name);
}
