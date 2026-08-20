/* Wrapper over the vendored signals engine — the only file that imports it, so the
   engine stays swappable (TC39 Signals later). Raw `effect` is deliberately not part of
   the public surface: every effect must have an owner (Scope), or it leaks. */
import {Signal} from '../../vendor/signals-core/signals-core.js';

export {signal, computed, batch, untracked, Signal} from '../../vendor/signals-core/signals-core.js';
export type {ReadonlySignal} from '../../vendor/signals-core/signals-core.js';
export {effect as rawEffect} from '../../vendor/signals-core/signals-core.js';

/** Whether `s` is a signal whose `value` can be written — the two-way gate. A prototype-descriptor
 * check: the vendored `Computed` defines a getter-only `value`, so computeds classify read-only
 * without any engine knowledge. */
export function isWritableSignal(s: unknown): s is Signal<unknown> {
  if (!(s instanceof Signal))
    return false;
  for (let proto = Object.getPrototypeOf(s); proto !== null; proto = Object.getPrototypeOf(proto)) {
    const descriptor = Object.getOwnPropertyDescriptor(proto, 'value');
    if (descriptor)
      return descriptor.set !== undefined;
  }
  return false;
}
