// Temporary instrumentation wrapper around @datagrok-libraries/test.
// Re-exports the library and overrides test/category/before/after to log
// timestamped START/END lines for each test body and before/after hook,
// so we can correlate test timelines with CI wall-clock and find gaps.

import {
  test as _test,
  category as _category,
  before as _before,
  after as _after,
  TestOptions,
  CategoryOptions,
} from '@datagrok-libraries/test/src/test';

export * from '@datagrok-libraries/test/src/test';
import {_rdKitModule} from '../utils/chem-common-rdkit';
import * as grok from 'datagrok-api/grok';

let _currentCategory = '';

// Temporary: reads the RDKit WASM heap size and the JS heap so we can plot
// cumulative memory growth across the suite and locate where it balloons —
// the suite's silent freeze looks like a WASM-heap OOM that floats run-to-run.
function memHint(): string {
  let wasmMb = '?';
  let jsMb = '?';
  try {
    const buf = (_rdKitModule as any)?.HEAPU8?.buffer ?? (_rdKitModule as any)?.HEAP8?.buffer;
    if (buf?.byteLength)
      wasmMb = (buf.byteLength / 1048576).toFixed(1);
  } catch {/* ignore */}
  try {
    const used = (performance as any)?.memory?.usedJSHeapSize;
    if (used)
      jsMb = (used / 1048576).toFixed(1);
  } catch {/* ignore */}
  return `wasm=${wasmMb}MB js=${jsMb}MB`;
}

function log(kind: 'TEST' | 'BEFORE' | 'AFTER', phase: 'START' | 'END',
  cat: string, name?: string, dtMs?: number): void {
  const ts = new Date().toISOString();
  const where = name ? `"${cat}" / "${name}"` : `"${cat}"`;
  const tail = dtMs !== undefined ? ` (${dtMs.toFixed(1)} ms)` : '';
  console.log(`[TIMING ${ts}] ${kind.padEnd(6)} ${phase.padEnd(5)} ${where}${tail} [MEM ${memHint()}]`);
}

export function category(name: string, tests_: () => void, options?: CategoryOptions): void {
  _currentCategory = name;
  _category(name, tests_, options);
}

export function test(name: string, fn: () => Promise<any>, options?: TestOptions): void {
  const cat = _currentCategory;
  _test(name, async () => {
    const t0 = performance.now();
    log('TEST', 'START', cat, name);
    try {
      return await fn();
    } finally {
      log('TEST', 'END', cat, name, performance.now() - t0);
      // Chem-wide afterEach: many tests open table views / viewers / filters and never close
      // them, so views, dataframes and their RDKit-backed caches accumulate across the suite
      // until the CI runner runs out of memory and the whole job is killed (the "hanging" tests).
      // closeAll() only closes VIEWS, not tables — the dataframes (with their attached
      // fingerprint/scaffold-highlight caches) stay registered in the shell and keep growing.
      // So also close every remaining table to release them. Tests that need shared state across
      // the category must set it up in before() (only the disabled 3d-hover category does).
      try {
        const viewsBefore = [...grok.shell.tableViews].length;
        const tablesBefore = grok.shell.tables.length;
        grok.shell.closeAll();
        for (const t of [...grok.shell.tables])
          grok.shell.closeTable(t);
        const viewsAfter = [...grok.shell.tableViews].length;
        const tablesAfter = grok.shell.tables.length;
        console.warn(`[LEAK-DIAG] afterEach cleanup: views ${viewsBefore}->${viewsAfter} tables ${tablesBefore}->${tablesAfter}`);
      } catch (e) {
        console.warn(`[LEAK-DIAG] afterEach cleanup threw: ${e}`);
      }
    }
  }, options);
}

export function before(fn: () => Promise<void>): void {
  const cat = _currentCategory;
  _before(async () => {
    const t0 = performance.now();
    log('BEFORE', 'START', cat);
    try {
      await fn();
    } finally {
      log('BEFORE', 'END', cat, undefined, performance.now() - t0);
    }
  });
}

export function after(fn: () => Promise<void>): void {
  const cat = _currentCategory;
  _after(async () => {
    const t0 = performance.now();
    log('AFTER', 'START', cat);
    try {
      await fn();
    } finally {
      log('AFTER', 'END', cat, undefined, performance.now() - t0);
    }
  });
}
