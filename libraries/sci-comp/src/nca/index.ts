/**
 * Public entry point for the `nca` namespace of @datagrok-libraries/sci-comp.
 *
 * Pure-math, stateless computation core for Non-Compartmental Analysis.
 * No dependency on `datagrok-api` — DataFrame integration lives in the
 * `packages/NCA/` package.
 *
 * Usage:
 *   import {nca} from '@datagrok-libraries/sci-comp';
 *   // nca.computeNca(inputs, rules)  ← added by Phase 1 Task 1.9
 */

export * from './core';
