// Pool multi-session lifecycle — sequential and parallel fits on a single
// long-lived WorkerPool.
//
// The WorkerPool is multi-session by design: every fit calls
// `setupAll(sessionId)` to prime its session in every slot, dispatches one
// `RunSeed` per starting point, and `dropSession(sessionId)` in
// WorkerExecutor.run's finally. The contract these tests pin:
//
//   - Sequential reuse: after fit A's dropSession, fit B sees no leaked
//     session state, threshold, or pi.canceled bookkeeping.
//   - Parallel reuse: dropSession(A) while B is active doesn't affect B;
//     setupAll(B) doesn't stomp on activeSetups[A]; runs for A and B
//     interleave via `pump` without cross-session contamination.
//
// Robustness by design (avoiding flakiness):
//   - Explicit `new WorkerPool(2)` per test — no singleton state, no
//     hardware-concurrency variation across CI machines.
//   - Byte-identity vs a serial baseline computed on a fresh ephemeral
//     pool. Output is invariant under slot dispatch order (seed-index
//     reorder finalize), so completion-order timing is invisible.
//   - Trivial bodies, samplesCount=4, default 20s/60s setup/run timeouts.
//   - No wall-clock or completion-order assertions.
//   - All pools disposed in finally.

import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {ExecutorArgs, LOSS, OptimizerInputsConfig, OptimizerOutputsConfig,
  WorkerExecutor, WorkerPool, runWithEphemeralPool} from './imports';
import {makeExpDecayFunc, makeMultiOutputFunc} from './script-fixtures';
import {assertResultParity} from './parity-assertions';
import {defaultNmSettings, noEarlyStopping, rangeBound} from './utils';

// Builds a worker-eligible ExecutorArgs against an exp-decay-style fixture.
// Each call evaluates `func.prepare(...)` to materialize the target DF, so
// per-test seeds and configs stay independent.
async function buildExpDecayArgs(seed: number): Promise<ExecutorArgs> {
  const func = makeExpDecayFunc();
  const targetCall = func.prepare({a: 2.5, b: 0.7, N: 20});
  await targetCall.call();
  const targetDf = targetCall.getParamValue('simulation') as DG.DataFrame;
  const inputBounds: OptimizerInputsConfig = {
    a: rangeBound(0.5, 5, 'a'),
    b: rangeBound(0.1, 2, 'b'),
    N: {type: 'const', value: 20},
  };
  const outputTargets: OptimizerOutputsConfig = [{
    propName: 'simulation',
    type: DG.TYPE.DATA_FRAME,
    target: targetDf,
    argName: 't',
    cols: [targetDf.col('y')!],
  }];
  return {
    objectiveFunc: async () => 0, // worker arm ignores this; main-arm fallback only
    inputsBounds: inputBounds,
    samplesCount: 4,
    settings: defaultNmSettings(),
    reproSettings: {reproducible: true, seed},
    earlyStoppingSettings: noEarlyStopping(),
    func,
    outputTargets,
    lossType: LOSS.RMSE,
  };
}

async function buildMultiOutputArgs(seed: number, lossType: LOSS = LOSS.RMSE): Promise<ExecutorArgs> {
  const func = makeMultiOutputFunc();
  const targetCall = func.prepare({a: 1.5, b: 0.75, N: 20});
  await targetCall.call();
  const targetDf = targetCall.getParamValue('simulation') as DG.DataFrame;
  const inputBounds: OptimizerInputsConfig = {
    a: rangeBound(0, 5, 'a'),
    b: rangeBound(0, 5, 'b'),
    N: {type: 'const', value: 20},
  };
  const outputTargets: OptimizerOutputsConfig = [{
    propName: 'simulation',
    type: DG.TYPE.DATA_FRAME,
    target: targetDf,
    argName: 't',
    cols: [targetDf.col('y1')!, targetDf.col('y2')!],
  }];
  return {
    objectiveFunc: async () => 0,
    inputsBounds: inputBounds,
    samplesCount: 4,
    settings: defaultNmSettings(),
    reproSettings: {reproducible: true, seed},
    earlyStoppingSettings: noEarlyStopping(),
    func,
    outputTargets,
    lossType,
  };
}

category('ComputeUtils: Fitting / Pool multi-session', () => {
  test('sequential_two_fits_same_pool_match_serial', async () => {
    // Fit A then fit B on the same pool, awaited in series. Each must
    // match its solo-on-fresh-pool baseline byte-identically. A leaked
    // session, target, or threshold from A would shift B's result.
    const argsA = await buildExpDecayArgs(42);
    const argsB = await buildMultiOutputArgs(99);

    const soloA = await runWithEphemeralPool(argsA);
    const soloB = await runWithEphemeralPool(argsB);

    const pool = new WorkerPool(2);
    try {
      const wkr = new WorkerExecutor(pool);
      const seqA = await wkr.run(argsA);
      const seqB = await wkr.run(argsB);
      assertResultParity(seqA, soloA, 'sequential_A');
      assertResultParity(seqB, soloB, 'sequential_B');
    } finally {
      pool.dispose();
    }
  });

  test('sequential_same_script_different_seeds', async () => {
    // Same body twice, different seeds. If session lifecycle keyed on
    // body identity instead of sessionId, fit B would inherit fit A's
    // bound parameters / target and produce A's result.
    const argsA = await buildExpDecayArgs(42);
    const argsB = await buildExpDecayArgs(99);

    const soloA = await runWithEphemeralPool(argsA);
    const soloB = await runWithEphemeralPool(argsB);

    const pool = new WorkerPool(2);
    try {
      const wkr = new WorkerExecutor(pool);
      const seqA = await wkr.run(argsA);
      const seqB = await wkr.run(argsB);
      assertResultParity(seqA, soloA, 'same_body_seed42');
      assertResultParity(seqB, soloB, 'same_body_seed99');
    } finally {
      pool.dispose();
    }
  });

  test('sequential_three_fits_no_drift', async () => {
    // Three fits in a row. A sessionId leak that grows by 1 per fit (e.g.
    // dropSession failing to remove the session from `slot.sessions`)
    // wouldn't be visible after one round-trip but might surface as
    // dispatch slot exhaustion or a stale-session match by the third.
    const argsA = await buildExpDecayArgs(42);
    const argsB = await buildMultiOutputArgs(99);
    const argsC = await buildExpDecayArgs(7);

    const soloA = await runWithEphemeralPool(argsA);
    const soloB = await runWithEphemeralPool(argsB);
    const soloC = await runWithEphemeralPool(argsC);

    const pool = new WorkerPool(2);
    try {
      const wkr = new WorkerExecutor(pool);
      const seqA = await wkr.run(argsA);
      const seqB = await wkr.run(argsB);
      const seqC = await wkr.run(argsC);
      assertResultParity(seqA, soloA, 'three_A');
      assertResultParity(seqB, soloB, 'three_B');
      assertResultParity(seqC, soloC, 'three_C');
    } finally {
      pool.dispose();
    }
  });

  test('parallel_two_fits_same_pool_match_serial', async () => {
    // Two fits via Promise.all on the same pool. Catches queue-ordering
    // bugs that misroute a run from A to B's session-set check, or
    // setupAll(B) overwriting activeSetups[A], or dropSession(A)
    // leaking into B's slot.sessions set.
    const argsA = await buildExpDecayArgs(42);
    const argsB = await buildMultiOutputArgs(99);

    const soloA = await runWithEphemeralPool(argsA);
    const soloB = await runWithEphemeralPool(argsB);

    const pool = new WorkerPool(2);
    try {
      const wkrA = new WorkerExecutor(pool);
      const wkrB = new WorkerExecutor(pool);
      const [parA, parB] = await Promise.all([wkrA.run(argsA), wkrB.run(argsB)]);
      assertResultParity(parA, soloA, 'parallel_A');
      assertResultParity(parB, soloB, 'parallel_B');
    } finally {
      pool.dispose();
    }
  });

  test('parallel_overlapping_setup_and_run', async () => {
    // Fire B without awaiting A; B's setupAll overlaps A's runs (or
    // potentially A's setupAll, depending on event-loop ordering). Both
    // must complete and match their solo baselines. Probes pump()
    // races where B's runs would dispatch to a slot not yet primed.
    const argsA = await buildExpDecayArgs(42);
    const argsB = await buildMultiOutputArgs(99);

    const soloA = await runWithEphemeralPool(argsA);
    const soloB = await runWithEphemeralPool(argsB);

    const pool = new WorkerPool(2);
    try {
      const wkr = new WorkerExecutor(pool);
      // No await between starts — both promises in flight before either
      // completes. Promise.all is the join.
      const pA = wkr.run(argsA);
      const pB = wkr.run(argsB);
      const [parA, parB] = await Promise.all([pA, pB]);
      assertResultParity(parA, soloA, 'overlap_A');
      assertResultParity(parB, soloB, 'overlap_B');
    } finally {
      pool.dispose();
    }
  });

  test('parallel_two_fits_different_loss_types', async () => {
    // RMSE in parallel with MAD. Catches any per-pool latching of
    // `useRmse` (it's per-session inside the worker, and per-call in
    // `accumulateLoss`, but a regression that hoisted it to module scope
    // would corrupt one of the two fits).
    const argsA = await buildMultiOutputArgs(42, LOSS.RMSE);
    const argsB = await buildMultiOutputArgs(99, LOSS.MAD);

    const soloA = await runWithEphemeralPool(argsA);
    const soloB = await runWithEphemeralPool(argsB);

    const pool = new WorkerPool(2);
    try {
      const wkrA = new WorkerExecutor(pool);
      const wkrB = new WorkerExecutor(pool);
      const [parA, parB] = await Promise.all([wkrA.run(argsA), wkrB.run(argsB)]);
      assertResultParity(parA, soloA, 'mixed_loss_RMSE');
      assertResultParity(parB, soloB, 'mixed_loss_MAD');
    } finally {
      pool.dispose();
    }
  });

  test('sequential_pool_survives_after_fits', async () => {
    // After several fits and dropSessions, the pool must remain
    // functional for a fresh fit. Asserts the pool isn't accidentally
    // disposing itself or losing slots after dropSession.
    const argsA = await buildExpDecayArgs(42);
    const argsB = await buildExpDecayArgs(7);

    const soloB = await runWithEphemeralPool(argsB);

    const pool = new WorkerPool(2);
    try {
      const wkr = new WorkerExecutor(pool);
      // Run A and discard — dropSession fired in finally.
      await wkr.run(argsA);
      // Slot count should still be 2 — no slots were torn down.
      expect(pool._slotsForTest(), 2, 'pool must keep all slots after a clean fit');
      // Now B on the same pool must produce its solo baseline.
      const seqB = await wkr.run(argsB);
      assertResultParity(seqB, soloB, 'survives_B');
    } finally {
      pool.dispose();
    }
  });
});
