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

let _currentCategory = '';

function log(kind: 'TEST' | 'BEFORE' | 'AFTER', phase: 'START' | 'END',
  cat: string, name?: string, dtMs?: number): void {
  const ts = new Date().toISOString();
  const where = name ? `"${cat}" / "${name}"` : `"${cat}"`;
  const tail = dtMs !== undefined ? ` (${dtMs.toFixed(1)} ms)` : '';
  console.log(`[TIMING ${ts}] ${kind.padEnd(6)} ${phase.padEnd(5)} ${where}${tail}`);
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
