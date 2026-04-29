// Byte-identity assertions for fitting results.
//
// Reused by `cross-executor-parity.ts` (main vs worker) and
// `pool-multi-session.ts` (sequential / parallel sessions on the same pool
// vs solo baseline). The math invariant they pin is the same: given the
// same seed list, same NM settings, same cost function, the OptimizationResult
// is bit-identical regardless of dispatch arrangement.

import {expect} from '@datagrok-libraries/test/src/test';
import type {Extremum, OptimizationResult} from './imports';

export function assertExtremumByteIdentical(a: Extremum, b: Extremum, label: string): void {
  expect(Object.is(a.cost, b.cost), true, `${label}: cost not bit-identical (${a.cost} vs ${b.cost})`);
  expect(a.iterCount, b.iterCount, `${label}: iterCount differs`);
  expect(a.point.length, b.point.length, `${label}: point.length differs`);
  for (let j = 0; j < a.point.length; ++j) {
    expect(Object.is(a.point[j], b.point[j]), true,
      `${label}: point[${j}] not bit-identical (${a.point[j]} vs ${b.point[j]})`);
  }
  const len = Math.min(a.iterCosts.length, b.iterCosts.length, a.iterCount);
  for (let j = 0; j < len; ++j) {
    expect(Object.is(a.iterCosts[j], b.iterCosts[j]), true,
      `${label}: iterCosts[${j}] not bit-identical`);
  }
}

export function assertResultParity(a: OptimizationResult, b: OptimizationResult, label: string): void {
  expect(a.extremums.length, b.extremums.length, `${label}: extremum count differs`);
  // No reordering tolerance: WorkerExecutor records replies by seed index
  // and finalizes in index order, matching MainExecutor exactly. The
  // surrounding `runOptimizer` cost-sort is a stable sort on bit-identical
  // costs, so the post-sort orders also match.
  for (let i = 0; i < a.extremums.length; ++i)
    assertExtremumByteIdentical(a.extremums[i], b.extremums[i], `${label}[${i}]`);
}
