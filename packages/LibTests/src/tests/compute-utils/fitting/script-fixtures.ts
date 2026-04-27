// Synthetic JS-language `DG.Func` builders for Phase 3a.
//
// Each helper returns a DG.Script created via DG.Script.create(...) with a
// known closed-form body. Reused across cost-functions and end-to-end tests.

import * as DG from 'datagrok-api/dg';

function buildScript(lines: string[]): DG.Func {
  return DG.Script.create(lines.join('\n'));
}

/** Two scalar outputs: y1 = a + b, y2 = a * b. */
export function makeScalarPairFunc(): DG.Func {
  return buildScript([
    '//name: ScalarPair',
    '//language: javascript',
    '//input: double a',
    '//input: double b',
    '//output: double y1',
    '//output: double y2',
    '',
    'y1 = a + b;',
    'y2 = a * b;',
  ]);
}

/** Single scalar output: y = a + b. */
export function makeSingleScalarFunc(): DG.Func {
  return buildScript([
    '//name: SingleScalar',
    '//language: javascript',
    '//input: double a',
    '//input: double b',
    '//output: double y',
    '',
    'y = a + b;',
  ]);
}

/**
 * Exponential decay: y = a * exp(-b * t), with N samples uniformly in [0, 5].
 * One DataFrame output `simulation` with columns [t, y].
 */
export function makeExpDecayFunc(): DG.Func {
  return buildScript([
    '//name: ExpDecay',
    '//language: javascript',
    '//input: double a',
    '//input: double b',
    '//input: int N = 20',
    '//output: dataframe simulation',
    '',
    'const tArr = new Float32Array(N);',
    'const yArr = new Float32Array(N);',
    'for (let i = 0; i < N; i++) {',
    '  tArr[i] = i / (N - 1) * 5;',
    '  yArr[i] = a * Math.exp(-b * tArr[i]);',
    '}',
    'simulation = DG.DataFrame.fromColumns([',
    "  DG.Column.fromFloat32Array('t', tArr),",
    "  DG.Column.fromFloat32Array('y', yArr),",
    ']);',
  ]);
}

/**
 * Two-output linear mixture: y = a*t + b*sin(c*t).
 * One DataFrame output `simulation` with columns [t, y]; one scalar `mean` (mean of y).
 */
export function makeLinearMixtureFunc(): DG.Func {
  return buildScript([
    '//name: LinearMixture',
    '//language: javascript',
    '//input: double a',
    '//input: double b',
    '//input: double c',
    '//input: int N = 30',
    '//output: dataframe simulation',
    '//output: double mean',
    '',
    'const tArr = new Float32Array(N);',
    'const yArr = new Float32Array(N);',
    'let sum = 0;',
    'for (let i = 0; i < N; i++) {',
    '  tArr[i] = i / (N - 1) * 4;',
    '  yArr[i] = a * tArr[i] + b * Math.sin(c * tArr[i]);',
    '  sum += yArr[i];',
    '}',
    'mean = sum / N;',
    'simulation = DG.DataFrame.fromColumns([',
    "  DG.Column.fromFloat32Array('t', tArr),",
    "  DG.Column.fromFloat32Array('y', yArr),",
    ']);',
  ]);
}

/**
 * Multi-column DataFrame output: simulation has [t, y1, y2] where
 * y1 = a*t and y2 = b*t.
 */
export function makeMultiOutputFunc(): DG.Func {
  return buildScript([
    '//name: MultiOutput',
    '//language: javascript',
    '//input: double a',
    '//input: double b',
    '//input: int N = 20',
    '//output: dataframe simulation',
    '',
    'const tArr = new Float32Array(N);',
    'const y1Arr = new Float32Array(N);',
    'const y2Arr = new Float32Array(N);',
    'for (let i = 0; i < N; i++) {',
    '  tArr[i] = i / (N - 1);',
    '  y1Arr[i] = a * tArr[i];',
    '  y2Arr[i] = b * tArr[i];',
    '}',
    'simulation = DG.DataFrame.fromColumns([',
    "  DG.Column.fromFloat32Array('t', tArr),",
    "  DG.Column.fromFloat32Array('y1', y1Arr),",
    "  DG.Column.fromFloat32Array('y2', y2Arr),",
    ']);',
  ]);
}

/**
 * Throws on inputs satisfying a predicate (e.g. negative `a`). Used to test
 * fault tolerance in NM (failed seeds → fails DF, others continue).
 */
export function makeThrowingFunc(): DG.Func {
  return buildScript([
    '//name: ThrowsOnNegA',
    '//language: javascript',
    '//input: double a',
    '//input: double b',
    '//output: double y',
    '',
    'if (a < 0) throw new Error("negative a not allowed");',
    'y = a + b;',
  ]);
}
