/* Wrapper over the vendored signals engine — the only file that imports it, so the
   engine stays swappable (TC39 Signals later). Raw `effect` is deliberately not part of
   the public surface: every effect must have an owner (Scope), or it leaks. */
export {signal, computed, batch, untracked, Signal} from '../../vendor/signals-core/signals-core.js';
export type {ReadonlySignal} from '../../vendor/signals-core/signals-core.js';
export {effect as rawEffect} from '../../vendor/signals-core/signals-core.js';
