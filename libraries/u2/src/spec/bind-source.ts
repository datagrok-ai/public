/* The DD5 binding protocol: a source resolves one path step at a time (bindStep) and
   enumerates what a picker may offer (bindProps) without allocating a single signal —
   laziness is the contract, not an optimization. */
import type {Signal} from '../core/signals.js';
import type {PropertyLike} from '../core/property-like.js';

/** One enumerated binding step: PropertyLike (typed, semType included) plus what the picker
 * needs — no probing, no signal creation during enumeration. */
export interface BindProp extends PropertyLike {
  /** bindStep(name) answers another BindSource — the picker renders an expandable node. */
  walkable?: boolean;
  /** The leaf is a writable signal — two-way capable. */
  writable?: boolean;
  /** What `bindStep('')` answers: a path that stops at the source binds here, so the picker
   * offers the source itself. */
  default?: boolean;
}

export interface BindSource {
  /** Resolve one step; a Signal ends the walk, a BindSource continues it, null is
   * unresolvable. '' asks for the source's default binding. */
  bindStep(name: string): Signal<unknown> | BindSource | null;
  /** What the binding picker enumerates. Must not allocate signals. */
  bindProps(): BindProp[];
}

export function isBindSource(x: unknown): x is BindSource {
  const source = x as BindSource | null;
  return typeof source?.bindStep === 'function' && typeof source?.bindProps === 'function';
}
