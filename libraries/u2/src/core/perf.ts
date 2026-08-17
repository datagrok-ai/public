import {Component} from './component.js';

export interface PerfEntry {
  label: string;
  /** Synchronous construction time, ms. */
  construct: number;
  /** Time from insertion into the document to the second animation frame, ms. */
  paint: number;
  total: number;
}

/** Construct/paint budgets for a component tree, accumulated so a harness can report a run. */
export class Perf {
  private static readonly _entries: PerfEntry[] = [];

  /** Builds, mounts, measures, then detaches and disposes. Paint is a double-rAF wait, which
   * means nothing without a compositor — outside a browser it stays 0. */
  static async measure(label: string, build: () => Component): Promise<PerfEntry> {
    const started = performance.now();
    const component = build();
    const construct = performance.now() - started;
    let paint = 0;
    if (!Perf._headless) {
      const mounted = performance.now();
      document.body.append(component.root);
      await new Promise<void>((resolve) =>
        requestAnimationFrame(() => requestAnimationFrame(() => resolve())));
      paint = performance.now() - mounted;
      component.root.remove();
    }
    component.dispose();
    const entry: PerfEntry = {label, construct, paint, total: construct + paint};
    Perf._entries.push(entry);
    return entry;
  }

  static report(): PerfEntry[] {
    return Perf._entries.slice();
  }

  private static get _headless(): boolean {
    const runtime = globalThis as {process?: {versions?: {node?: string}}};
    return runtime.process?.versions?.node !== undefined;
  }
}
