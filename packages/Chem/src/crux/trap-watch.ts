/**
 * The trap watch of a Crux worker (crux.worker.ts). The RDKit worker's code the Crux workers run (rdkit-service/:
 * conversions, CoordGen, flattening, beautifying, MMPA, R-groups) catches each molecule's error (`catch {}`), a wasm
 * trap's too, and goes on with the next molecule on an instance the trap may have broken: a stack overflow leaves
 * Rust's stack pointer where the overflow left it, so later deep calls overflow sooner (a 3,000-atom chain's SMILES,
 * written before a 6,000-atom chain's overflow, overflows after it). So every call into the module goes through the
 * watch, which notes a trap even when its caller swallows it, and from then on the worker fails every job, as a trap:
 * CruxService replaces the workers (a trap restarts the pool) and RDKit does the jobs' work.
 */

import {isWasmTrap} from './crux-api';

/** The first trap of the instance: it stays, the instance being broken (CruxService replaces the worker). */
let noted: string | null = null;

/** The first trap noted (its message), or null. */
export function trapNoted(): string | null {
  return noted;
}

type AnyFunction = (...args: any[]) => any;

function watched<F extends AnyFunction>(f: F): F {
  return function watchedCall(this: unknown, ...args: unknown[]) {
    try {
      return f.apply(this, args); // eslint-disable-line no-invalid-this
    } catch (e) {
      if (isWasmTrap(e))
        noted ??= `${(e as Error).name}: ${(e as Error).message}`;
      throw e;
    }
  } as F;
}

/** Whether `f` is a class of the glue (a JSMol, a Collection, ...): its prototype has methods. */
const isClass = (f: Function): boolean =>
  !!f.prototype && Object.getOwnPropertyNames(f.prototype).some((name) => name !== 'constructor');

/**
 * The glue's exports with every call watched: its functions wrapped in a copy of the namespace; its classes' methods
 * wrapped in place, on their prototypes, once (so the objects the module returns are watched wherever they are used).
 */
export function watchTraps<T extends object>(glue: T): T {
  const copy: {[name: string]: unknown} = {};
  for (const [name, value] of Object.entries(glue)) {
    if (typeof value !== 'function')
      copy[name] = value;
    else if (isClass(value)) {
      for (const method of Object.getOwnPropertyNames(value.prototype)) {
        const descriptor = Object.getOwnPropertyDescriptor(value.prototype, method)!;
        if (method !== 'constructor' && typeof descriptor.value === 'function' && descriptor.writable)
          value.prototype[method] = watched(descriptor.value);
      }
      copy[name] = value;
    } else
      copy[name] = watched(value as AnyFunction);
  }
  return copy as T;
}
