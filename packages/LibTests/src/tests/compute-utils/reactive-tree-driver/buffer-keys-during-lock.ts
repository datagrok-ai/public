import {category, test, before} from '@datagrok-libraries/test/src/test';
import {TestScheduler} from 'rxjs/testing';
import {BehaviorSubject, merge} from 'rxjs';
import {map, switchMap} from 'rxjs/operators';
import {bufferKeysDuringLock} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import {createTestScheduler} from '../../../test-utils';


category('ComputeUtils: Driver bufferKeysDuringLock', async () => {
  let testScheduler: TestScheduler;

  before(async () => {
    testScheduler = createTestScheduler();
  });

  test('Pass through while unlocked', async () => {
    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const lock$ = cold('f-------', {f: false}) as any;
      const source$ = cold('-a-b-c-', {
        a: ['A', 1] as const,
        b: ['B', 2] as const,
        c: ['C', 3] as const,
      });
      const result$ = source$.pipe(bufferKeysDuringLock<string, number>(lock$));
      expectObservable(result$).toBe('-a-b-c-', {
        a: ['A', 1],
        b: ['B', 2],
        c: ['C', 3],
      });
    });
  });

  test('Buffer during lock and flush on unlock', async () => {
    testScheduler.run((helpers) => {
      const {cold, hot, expectObservable} = helpers;
      const lock$ = hot(' f--t--------f--', {t: true, f: false});
      const source$ = hot('---a--b--c-----', {
        a: ['A', 1] as const,
        b: ['B', 2] as const,
        c: ['C', 3] as const,
      });
      const result$ = source$.pipe(bufferKeysDuringLock<string, number>(lock$));
      expectObservable(result$).toBe('------------(abc)', {
        a: ['A', 1],
        b: ['B', 2],
        c: ['C', 3],
      });
    });
  });

  test('Deduplicate by key keeping latest value', async () => {
    testScheduler.run((helpers) => {
      const {cold, hot, expectObservable} = helpers;
      const lock$ = hot(' f-t-----------f--', {t: true, f: false});
      const source$ = hot('--a--b--c--d-----', {
        a: ['A', 1] as const,
        b: ['B', 2] as const,
        c: ['A', 3] as const,
        d: ['B', 4] as const,
      });
      const result$ = source$.pipe(bufferKeysDuringLock<string, number>(lock$));
      expectObservable(result$).toBe('--------------(cd)', {
        c: ['A', 3],
        d: ['B', 4],
      });
    });
  });

  test('No emissions during lock produces nothing on unlock', async () => {
    testScheduler.run((helpers) => {
      const {cold, hot, expectObservable} = helpers;
      const lock$ = hot(' f-t---f----', {t: true, f: false});
      const source$ = hot('a---------b', {
        a: ['A', 1] as const,
        b: ['B', 2] as const,
      });
      const result$ = source$.pipe(bufferKeysDuringLock<string, number>(lock$));
      expectObservable(result$).toBe('a---------b', {
        a: ['A', 1],
        b: ['B', 2],
      });
    });
  });

  test('Multiple lock/unlock cycles', async () => {
    testScheduler.run((helpers) => {
      const {hot, expectObservable} = helpers;
      const lock$ = hot(' f-t-----f--t-----f--', {t: true, f: false});
      const source$ = hot('--a--b-----c--d-----', {
        a: ['A', 1] as const,
        b: ['B', 2] as const,
        c: ['C', 3] as const,
        d: ['D', 4] as const,
      });
      const result$ = source$.pipe(bufferKeysDuringLock<string, number>(lock$));
      expectObservable(result$).toBe('--------(ab)-----(cd)', {
        a: ['A', 1],
        b: ['B', 2],
        c: ['C', 3],
        d: ['D', 4],
      });
    });
  });

  test('Works with BehaviorSubject lock', async () => {
    testScheduler.run((helpers) => {
      const {cold, expectObservable} = helpers;
      const lock$ = new BehaviorSubject(false);
      const source$ = cold('-a-b-c-', {
        a: ['X', 10] as const,
        b: ['Y', 20] as const,
        c: ['Z', 30] as const,
      });
      const result$ = source$.pipe(bufferKeysDuringLock<string, number>(lock$));
      expectObservable(result$).toBe('-a-b-c-', {
        a: ['X', 10],
        b: ['Y', 20],
        c: ['Z', 30],
      });
    });
  });

  test('Passthrough resumes after unlock', async () => {
    testScheduler.run((helpers) => {
      const {hot, expectObservable} = helpers;
      const lock$ = hot(' f-t---f------', {t: true, f: false});
      const source$ = hot('--a-b--c--d--', {
        a: ['A', 1] as const,
        b: ['A', 2] as const,
        c: ['C', 3] as const,
        d: ['D', 4] as const,
      });
      const result$ = source$.pipe(bufferKeysDuringLock<string, number>(lock$));
      // a and b have same key A, deduplicates to b; flushed at t=60
      // c at t=70 and d at t=100 pass through
      expectObservable(result$).toBe('------bc--d--', {
        b: ['A', 2],
        c: ['C', 3],
        d: ['D', 4],
      });
    });
  });

  test('Lock starts as true (initial lock)', async () => {
    testScheduler.run((helpers) => {
      const {hot, expectObservable} = helpers;
      // Lock starts true — simulates tree init where globalROLocked$ starts locked
      const lock$ = hot(' t--------f----', {t: true, f: false});
      const source$ = hot('-a--b--c-------', {
        a: ['A', 1] as const,
        b: ['B', 2] as const,
        c: ['C', 3] as const,
      });
      const result$ = source$.pipe(bufferKeysDuringLock<string, number>(lock$));
      expectObservable(result$).toBe('---------(abc)', {
        a: ['A', 1],
        b: ['B', 2],
        c: ['C', 3],
      });
    });
  });

  test('Lock starts as true with dedup', async () => {
    testScheduler.run((helpers) => {
      const {hot, expectObservable} = helpers;
      // Same key updated multiple times during initial lock
      const lock$ = hot(' t-----------f----', {t: true, f: false});
      const source$ = hot('-a--b--c--d-------', {
        a: ['X', 1] as const,
        b: ['Y', 2] as const,
        c: ['X', 3] as const,
        d: ['Y', 4] as const,
      });
      const result$ = source$.pipe(bufferKeysDuringLock<string, number>(lock$));
      expectObservable(result$).toBe('------------(cd)', {
        c: ['X', 3],
        d: ['Y', 4],
      });
    });
  });

  test('Source completes while locked flushes buffer', async () => {
    testScheduler.run((helpers) => {
      const {hot, expectObservable} = helpers;
      const lock$ = hot(' t---------f---', {t: true, f: false});
      const source$ = hot('-a--b-|', {
        a: ['A', 1] as const,
        b: ['B', 2] as const,
      });
      const result$ = source$.pipe(bufferKeysDuringLock<string, number>(lock$));
      // Source completes at frame 6 while locked.
      // Window closes on source completion → reduce emits buffered items at frame 6.
      expectObservable(result$).toBe('------(ab|)', {
        a: ['A', 1],
        b: ['B', 2],
      });
    });
  });

  test('Source error propagates while unlocked', async () => {
    testScheduler.run((helpers) => {
      const {hot, expectObservable} = helpers;
      const lock$ = new BehaviorSubject(false);
      const source$ = hot('-a-#', {
        a: ['A', 1] as const,
      });
      const result$ = source$.pipe(bufferKeysDuringLock<string, number>(lock$));
      expectObservable(result$).toBe('-a-#', {
        a: ['A', 1],
      });
    });
  });

  test('Rapid lock toggle', async () => {
    testScheduler.run((helpers) => {
      const {hot, expectObservable} = helpers;
      // lock:   f--t---f--t---f---
      // source: ---a------b-------
      // a buffered in first lock window, flushed at frame 7
      // b buffered in second lock window, flushed at frame 14
      const lock$ = hot(' f--t---f--t---f---', {t: true, f: false});
      const source$ = hot('---a------b-------', {
        a: ['A', 1] as const,
        b: ['B', 2] as const,
      });
      const result$ = source$.pipe(bufferKeysDuringLock<string, number>(lock$));
      expectObservable(result$).toBe('-------a------b---', {
        a: ['A', 1],
        b: ['B', 2],
      });
    });
  });

  test('Simulates Compute2 switchMap pattern', async () => {
    testScheduler.run((helpers) => {
      const {hot, expectObservable} = helpers;
      const lock$ = hot('t------------f-----', {t: true, f: false});
      // Simulates Driver emitting Record<string, Observable<T>>,
      // switchMap resetting state, then merging inner observables
      // Individual [k,v] pairs arrive during lock
      const inner1$ = hot('--a---b-----------', {
        a: ['step1', 'running'] as const,
        b: ['step2', 'idle'] as const,
      });
      const inner2$ = hot('----c-----d-------', {
        c: ['step1', 'done'] as const,
        d: ['step2', 'running'] as const,
      });
      const source$ = merge(inner1$, inner2$);
      const result$ = source$.pipe(bufferKeysDuringLock<string, string>(lock$));
      // step1: running→done (dedup keeps done)
      // step2: idle→running (dedup keeps running)
      expectObservable(result$).toBe('-------------(cd)', {
        c: ['step1', 'done'],
        d: ['step2', 'running'],
      });
    });
  });

  test('Simulates tree replacement with switchMap', async () => {
    testScheduler.run((helpers) => {
      const {hot, expectObservable} = helpers;
      // First tree active and unlocked, then new tree init locks
      const lock$ = hot(' f----t---------f---', {t: true, f: false});
      // Old tree keys pass through, new tree keys arrive during lock
      const source$ = hot('-a-b---c--d--e-----', {
        a: ['old1', 10] as const,
        b: ['old2', 20] as const,
        c: ['new1', 30] as const,
        d: ['new2', 40] as const,
        e: ['new1', 50] as const,
      });
      const result$ = source$.pipe(bufferKeysDuringLock<string, number>(lock$));
      // old1(1), old2(3) pass through; new1 deduped 30→50, new2=40 flushed at 15
      expectObservable(result$).toBe('-a-b-----------(ef)', {
        a: ['old1', 10],
        b: ['old2', 20],
        e: ['new1', 50],
        f: ['new2', 40],
      });
    });
  });
});
