import {Signal, ReadonlySignal} from './signals.js';
import {Scope} from './scope.js';
import {text} from './text.js';

/** `null` and `undefined` show as empty. */
export function bindText(scope: Scope, el: Element, source: ReadonlySignal<unknown>): void {
  scope.effect(() => {
    el.textContent = text(source.value);
  });
}

/** Two-way binding with centralized echo suppression: a user edit updates the signal
 * without the effect writing the (identical) value back into the input mid-keystroke.
 *
 * `commitOn` picks which of the browser's two edit events writes the signal: every keystroke
 * (`input`, the default) or the commit — blur or Enter — the way Dart's `fireChanged` reports one
 * (`commitOn: 'change'`). Either way the field shows what is being typed; only the signal waits. */
export function bindValue(scope: Scope, input: HTMLInputElement | HTMLTextAreaElement,
  source: Signal<string>, commitOn: 'input' | 'change' = 'input'): void {
  // set by a commit that changed the signal, consumed by the effect run it causes — which comes
  // after the handler, when its batch ends
  let echo = false;
  scope.effect(() => {
    // a bound source may hold nothing yet — a data source before its first run — and the DOM
    // would render that as the literal "undefined"
    const v = source.value ?? '';
    if (echo)
      echo = false;
    else
      input.value = v;
  });
  scope.listen(input, commitOn, () => {
    const before = source.peek();
    source.value = input.value;
    echo = !Object.is(source.peek(), before);
  });
}
