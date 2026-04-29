// Bounds checker shared by the main-arm sampler and the fitting worker.
//
// DG-free on purpose: imports `compileFormula` (zero-dependency) and uses
// `import type` for the bound shapes so the worker bundle never reaches
// `optimizer-misc.ts` (which has a value-level `import * as DG`).

import {compileFormula, runFormula, type CompiledFormula} from './formulas-resolver';
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

// One precompiled bound. Both sides resolved as closures — value sides
// bind the constant; formula sides go through compileFormula. `pos` is
// pre-resolved so the hot loop avoids the per-call Map lookup.
type CompiledBound = {
  name: string;
  pos: number;
  minOf: CompiledFormula;
  maxOf: CompiledFormula;
};

function compileSide(side: ChangingValue['bottom'], knownNames: readonly string[]): CompiledFormula {
  if (side.type === 'value') {
    const v = side.value;
    return () => v;
  }
  return compileFormula(side.formula, knownNames);
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

  // Names a formula may reference: every fixed input + every varied input.
  // (A formula at iteration k may legally reference any varied input from
  // iteration < k. Compiling against the union is fine — referenced names
  // get destructured; ones not actually referenced drop out at compile
  // time via the regex intersect inside compileFormula.)
  const knownNames = [...Object.keys(contextFixed), ...variedInputNames];

  // Precompile every bound at fit-setup time. Hot loop becomes a flat
  // array walk with one closure call per side, no parse / `with` cost.
  const compiled: CompiledBound[] = [];
  for (const [name, , bound] of [...nonFormulaBounds, ...formulaBounds]) {
    compiled.push({
      name,
      pos: variedNameToPosition.get(name)!,
      minOf: compileSide(bound.bottom, knownNames),
      maxOf: compileSide(bound.top, knownNames),
    });
  }

  // Reuse one context object across calls instead of `{...contextFixed}`
  // per call. Varied slots are reset each call so a formula referencing a
  // not-yet-iterated varied input gets `undefined` (matching prior
  // semantics — it would have been missing from a fresh clone).
  const context: Record<string, any> = {...contextFixed};

  return function isInsideBounds(point: Float64Array): boolean {
    for (const name of variedInputNames) context[name] = undefined;
    for (let i = 0; i < compiled.length; ++i) {
      const c = compiled[i];
      const val = point[c.pos];
      const min = c.minOf(context);
      const max = c.maxOf(context);
      if (max == null || min == null || val == null || min > max || val < min || val > max)
        return false;
      context[c.name] = val;
    }
    return true;
  };
}
