// Bounds checker shared by the main-arm sampler and the fitting worker.
//
// DG-free on purpose: imports `runFormula` (zero-dependency) and uses
// `import type` for the bound shapes so the worker bundle never reaches
// `optimizer-misc.ts` (which has a value-level `import * as DG`).

import {runFormula} from './formulas-resolver';
import type {BoundFormula, ChangingValue, ConstValue, ValueBoundsData} from './optimizer-misc';

type AccData = {
  constValues: [string, ConstValue][];
  nonFormulaBounds: [string, number, ChangingValue][];
  formulaBounds: [string, number, ChangingValue][];
  boundsIdx: number;
};

export function getAccData(inputs: Record<string, ValueBoundsData>): AccData {
  // Partition inputs by type (fixed / numeric / formula) while keeping the
  // original ordering index of varied inputs.
  return Object.entries(inputs).reduce((acc, [name, val]) => {
    if (val.type === 'const')
      acc.constValues.push([name, val]);
    else if (val.top.type === 'value' && val.bottom.type === 'value') {
      acc.nonFormulaBounds.push([name, acc.boundsIdx, val]);
      acc.boundsIdx++;
    } else {
      acc.formulaBounds.push([name, acc.boundsIdx, val]);
      acc.boundsIdx++;
    }
    return acc;
  }, {constValues: [], nonFormulaBounds: [], formulaBounds: [], boundsIdx: 0} as AccData);
}

export function getFixedContext(constValues: [string, ConstValue][]): Record<string, any> {
  const contextFixed: Record<string, any> = {};
  for (const [name, data] of constValues)
    contextFixed[name] = data.value;
  return contextFixed;
}

export function evalBoundFormula(bound: BoundFormula, context: Record<string, any>): number | null {
  return runFormula(bound.formula, context);
}

export function makeBoundsChecker(
  inputs: Record<string, ValueBoundsData>,
  variedInputNames: string[],
  fixedContextOverride?: Record<string, any>,
): (point: Float64Array) => boolean {
  const variedNameToPosition = new Map(variedInputNames.map((name, pos) => [name, pos]));
  const {constValues, nonFormulaBounds, formulaBounds} = getAccData(inputs);
  // Worker arm passes a pre-reified fixed context built from
  // setup.fixedInputs + setup.fixedDataFrames so const Dayjs/Date/DataFrame
  // values seen by formulas match the main arm's. Main arm omits the
  // override and falls back to the values stored in `inputs` const entries.
  const contextFixed = fixedContextOverride ?? getFixedContext(constValues);

  return function isInsideBounds(point: Float64Array): boolean {
    const context = {...contextFixed};
    for (const [name, , bound] of [...nonFormulaBounds, ...formulaBounds]) {
      const pos = variedNameToPosition.get(name)!;
      const val = point[pos];
      const min = bound.bottom.type === 'value' ? bound.bottom.value : evalBoundFormula(bound.bottom, context);
      const max = bound.top.type === 'value' ? bound.top.value : evalBoundFormula(bound.top, context);
      if (max == null || min == null || val == null || min > max || val < min || val > max)
        return false;
      context[name] = val;
    }
    return true;
  };
}
