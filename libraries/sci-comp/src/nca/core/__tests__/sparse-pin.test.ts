/**
 * Byte-identity pin of `sparseAuc` (sparse-iv-bolus-anchor AC-2.1).
 *
 * The snapshot was written on sci-comp 0.11.0, BEFORE the IV-bolus anchor landed. It pins
 * every number of the PO mouse composites and the batch oracle input, and of the IV-bolus
 * composites run WITHOUT a route (the default keeps today's `(0, 0)` anchor). After the
 * change the same calls, and the same calls with an extravascular route, must reproduce
 * it exactly. Warning codes are pinned, not their messages: the terminal note's message
 * gains its estimated excess (F5), its presence does not depend on the anchor.
 */

import {sparseAuc} from '../sparse';
import type {SparseAucOptions, SparseAucResult, SparseInput} from '../sparse';
import {ROUTE_PO} from '../types';
import {mouseInputs, ivInputs} from './sparse-fixture-io';

function pinned(r: SparseAucResult): unknown {
  return {
    auc: r.auc, se: r.se, df: r.df, ci: r.ci, topology: r.topology,
    composite: r.composite, warningCodes: r.warnings.map((w) => w.code),
  };
}

function cases(): Array<{name: string; input: SparseInput}> {
  return [
    ...mouseInputs().map(({dose, input}) => ({name: `PO mouse ${dose} mg/kg`, input})),
    ...ivInputs().map(({block, input}) => ({name: `IV-bolus composite ${block.dose_mg_kg} mg/kg`, input})),
  ];
}

describe('sparseAuc — 0.11.0 pin (AC-2.1)', () => {
  for (const {name, input} of cases()) {
    it(`${name}: no route`, () => {
      expect(pinned(sparseAuc(input, {blqRule: 'set-zero'}))).toMatchSnapshot();
    });
    it(`${name}: no route, t½ supplied`, () => {
      expect(pinned(sparseAuc(input, {blqRule: 'set-zero', estimatedHalfLife: 2}))).toMatchSnapshot();
    });
  }

  // An extravascular route must give the snapshot of the no-route call.
  for (const {name, input} of mouseInputs().map(({dose, input}) => ({name: `PO mouse ${dose} mg/kg`, input}))) {
    it(`${name}: extravascular route equals no route`, () => {
      const opts: SparseAucOptions = {blqRule: 'set-zero', route: ROUTE_PO};
      expect(pinned(sparseAuc(input, opts))).toEqual(pinned(sparseAuc(input, {blqRule: 'set-zero'})));
    });
  }
});
