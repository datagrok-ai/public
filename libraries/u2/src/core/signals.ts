import {Signal} from 'datagrok-api/u2core';

export {signal, computed, batch, untracked, Signal, rawEffect, isWritableSignal} from 'datagrok-api/u2core';
export type {ReadonlySignal} from 'datagrok-api/u2core';

// Two engines on a page split every `instanceof Signal`, ambient adoption and the leak counts.
const dg = (globalThis as {DG?: {U2?: {Signal?: unknown}}}).DG;
if (dg?.U2?.Signal !== undefined && dg.U2.Signal !== Signal) {
  throw new Error('u2: a second copy of the signal engine is bundled — build with @datagrok/build-config ' +
    '(it externalizes datagrok-api/u2core as DG.U2)');
}
