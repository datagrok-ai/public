/* The DD5 binding protocol: a source resolves one path step at a time (bindStep) and
   enumerates what a picker may offer (bindProps) without allocating a single signal —
   laziness is the contract, not an optimization. The shapes are the core's canonical
   `BindProp`/`BindSource`; this module keeps the spec-layer names and the guard. */
import type {BindSource} from '../core/widget-like.js';

export type {BindProp, BindSource} from '../core/widget-like.js';

export function isBindSource(x: unknown): x is BindSource {
  const source = x as BindSource | null;
  return typeof source?.bindStep === 'function' && typeof source?.bindProps === 'function';
}
