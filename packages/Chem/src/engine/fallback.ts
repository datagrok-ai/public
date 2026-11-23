/**
 * Every Crux path has an RDKit path behind it: a Crux failure, a "not supported by Crux" refusal or a molecule Crux
 * cannot read sends that call (that molecule, when drawing) to RDKit. Each call name is logged once per session: on
 * the console, and to the platform's log (the package logger), so the switch to Crux can be judged on real usage.
 */
import {_package} from '../package';
import {noteCruxError} from './crux-module';

const warned = new Set<string>();
const counts = new Map<string, number>();
const reasons = new Map<string, string>();
/** How many {@link quietFallbacksForTests} runs are under way. */
let quiet = 0;

/** The first fallback of each call name per session, to the platform's log (a warning with the call, the reason and
 * the setting's engine); never to the UI. Tests replace it ({@link setFallbackReporterForTests}). */
let report = (call: string, message: string): void => {
  try {
    _package.logger.warning(`Crux fallback: ${call}`, {call, reason: message, engine: 'Crux', fallback: 'RDKit'});
  } catch {
    // no package logger (a test bundle without the published package): the console line stands
  }
};

/** Logs, once per `call`, that RDKit did what Crux could not (and counts every time). */
export function warnCruxFallback(call: string, reason: unknown): void {
  counts.set(call, (counts.get(call) ?? 0) + 1);
  const message = reason instanceof Error ? reason.message : String(reason);
  reasons.set(call, message);
  if (warned.has(call) || quiet > 0)
    return;
  warned.add(call);
  console.warn(`Chem | Crux could not do ${call}, RDKit did it (logged once per call): ${message}`);
  report(call, message);
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

/**
 * A column job's RDKit side for the rows its Crux side left (`rows`: a molecule Crux does not read, a call it
 * refuses): `rdkit(rows)` computes them, `put(patch, k, row)` stores the k-th
 * result at `row` and tells whether RDKit did it (a row neither engine does is no fallback); the rows RDKit did are
 * logged as one fallback of `call` (`example`: what the log names).
 */
export async function rdkitForRows<P>(call: string, rows: number[], rdkit: (rows: number[]) => Promise<P>,
  put: (patch: P, k: number, row: number) => boolean, example: (row: number) => string): Promise<void> {
  if (!rows.length)
    return;
  const patch = await rdkit(rows);
  const done = rows.filter((row, k) => put(patch, k, row));
  if (done.length)
    warnCruxFallback(call, `RDKit did ${done.length} of the ${rows.length} rows Crux left, e.g. ${example(done[0])}`);
}

/** How many times RDKit did `call` for Crux so far (tests). */
export function cruxFallbackCount(call: string): number {
  return counts.get(call) ?? 0;
}

/** Why Crux could not do `call` the last time RDKit did it, in its own words (tests). */
export function cruxFallbackReason(call: string): string | undefined {
  return reasons.get(call);
}

/** Tests: `run()` with its fallbacks counted (and their reasons kept) but not logged: a fallback a test makes on
 * purpose neither reaches the platform's log, where it reads as one of real use, nor keeps the first real one of its
 * call in the session from reaching it (each call name is logged once). */
export async function quietFallbacksForTests<T>(run: () => Promise<T>): Promise<T> {
  quiet++;
  try {
    return await run();
  } finally {
    quiet--;
  }
}

/** Tests: where the first fallback of each call name goes (the package logger's warning); returns the one replaced. */
export function setFallbackReporterForTests(reporter: (call: string, message: string) => void):
  (call: string, message: string) => void {
  const previous = report;
  report = reporter;
  return previous;
}
