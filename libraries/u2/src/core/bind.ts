import {Signal, ReadonlySignal} from './signals.js';
import {Ownership} from './owner.js';
import type {Owner} from './owner.js';
import {text} from './text.js';

/** `null` and `undefined` show as empty. Followed while `owner` lives — whoever is building by default. */
export function bindText(el: Element, source: ReadonlySignal<unknown>, options?: {owner?: Owner}): void {
  Ownership.scopeOf(options?.owner, 'bindText()').effect(() => {
    el.textContent = text(source.value);
  });
}

/** Two-way binding with centralized echo suppression: a user edit updates the signal
 * without the effect writing the (identical) value back into the input mid-keystroke.
 *
 * `commitOn` picks which of the browser's two edit events writes the signal: every keystroke
 * (`input`, the default) or the commit — blur or Enter — the way Dart's `fireChanged` reports one
 * (`commitOn: 'change'`). Either way the field shows what is being typed; only the signal waits. Owned by
 * `owner`, whoever is building by default. */
export function bindValue(input: HTMLInputElement | HTMLTextAreaElement, source: Signal<string>,
  options?: {owner?: Owner, commitOn?: 'input' | 'change'}): void {
  const scope = Ownership.scopeOf(options?.owner, 'bindValue()');
  const commitOn = options?.commitOn ?? 'input';
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
