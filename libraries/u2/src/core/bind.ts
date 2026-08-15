import {Signal, ReadonlySignal} from './signals.js';
import {Scope} from './scope.js';

export function bindText(scope: Scope, el: Element, source: ReadonlySignal<unknown>): void {
  scope.effect(() => {
    el.textContent = String(source.value);
  });
}

/** Two-way binding with centralized echo suppression: a user edit updates the signal
 * without the effect writing the (identical) value back into the input mid-keystroke. */
export function bindValue(scope: Scope, input: HTMLInputElement | HTMLTextAreaElement, source: Signal<string>): void {
  let echo = false;
  scope.effect(() => {
    const v = source.value;
    if (!echo)
      input.value = v;
  });
  const handler = () => {
    echo = true;
    try {
      source.value = input.value;
    } finally {
      echo = false;
    }
  };
  input.addEventListener('input', handler);
  scope.own(() => input.removeEventListener('input', handler));
}
