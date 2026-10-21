/**
 * Every Crux path has an RDKit path behind it: a Crux failure, a "not supported by Crux" refusal or a molecule Crux
 * cannot read sends that call (that molecule, when drawing) to RDKit. Each call name is logged once per session.
 */
import {noteCruxError} from './crux-module';

const warned = new Set<string>();
const counts = new Map<string, number>();

/** Logs, once per `call`, that RDKit did what Crux could not (and counts every time). */
export function warnCruxFallback(call: string, reason: unknown): void {
  counts.set(call, (counts.get(call) ?? 0) + 1);
  if (warned.has(call))
    return;
  warned.add(call);
  const message = reason instanceof Error ? reason.message : String(reason);
  console.warn(`Chem | Crux could not do ${call}, RDKit did it (logged once per call): ${message}`);
}

/** `crux()`, or `rdkit()` when it throws. */
export function cruxOrRdKit<T>(call: string, crux: () => T, rdkit: () => T): T {
  try {
    return crux();
  } catch (e) {
    noteCruxError(e);
    warnCruxFallback(call, e);
    return rdkit();
  }
}

/** `crux()`, or `rdkit()` when it rejects. */
export async function cruxOrRdKitAsync<T>(call: string, crux: () => Promise<T>, rdkit: () => Promise<T>): Promise<T> {
  try {
    return await crux();
  } catch (e) {
    noteCruxError(e);
    warnCruxFallback(call, e);
    return rdkit();
  }
}

/** How many times RDKit did `call` for Crux so far (tests). */
export function cruxFallbackCount(call: string): number {
  return counts.get(call) ?? 0;
}
