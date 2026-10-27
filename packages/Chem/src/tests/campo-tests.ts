import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';
import {
  DESIRABILITY_PROFILE_TYPE, DesirabilityProfile, MPO_WEIGHTING_CAPTIONS, MpoCalculator, NumericalDesirability,
  desirabilityKey, mpo,
} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';

import {campoWeights, cholesky, finalWeights, pearsonMatrix, prepareDesirabilities, robustInverse}
  from '../mpo/campo-math';
import {stageAKey} from '../mpo/campo-controller';
import {MpoProfileCreateView} from '../mpo/mpo-create-profile';
import {fromRow, toRow} from '../mpo/mpo-profile-store';
import {applyDesirabilityTags} from '../mpo/utils';
import {ProfileRow} from '../generated/db';

/** m×m correlation matrix from its upper triangle, row by row. */
function corr(m: number, upper: number[]): Float64Array {
  const R = new Float64Array(m * m);
  let k = 0;
  for (let i = 0; i < m; i++) {
    R[i * m + i] = 1;
    for (let j = i + 1; j < m; j++)
      R[i * m + j] = R[j * m + i] = upper[k++];
  }
  return R;
}

function weightsOf(R: Float64Array, m: number): Float64Array {
  const inv = robustInverse(R, m);
  expect(inv != null, true, 'R is invertible');
  return campoWeights(inv!.Rinv, m, Int32Array.from({length: m}, (_, i) => i), m).base;
}

/** Deterministic pseudo-random numbers in [0, 1). */
function rng(seed: number): () => number {
  let s = seed >>> 0;
  return () => (s = (s * 1664525 + 1013904223) >>> 0) / 4294967296;
}

/** n rows of d columns in [0, 1]; column k leans on column k−1, so neighbours correlate. */
function correlatedColumns(n: number, d: number, seed = 1): Float32Array[] {
  const next = rng(seed);
  const cols: Float32Array[] = [];
  for (let k = 0; k < d; k++) {
    const c = new Float32Array(n);
    for (let i = 0; i < n; i++)
      c[i] = k === 0 ? next() : 0.6 * cols[k - 1][i] + 0.4 * next();
    cols.push(c);
  }
  return cols;
}

const numerical = (line: number[][], p: Partial<NumericalDesirability> = {}): NumericalDesirability =>
  ({functionType: 'numerical', weight: 1, line, ...p});

/* Correlation-aware MPO weights (caMPO): the computation of campo-math.ts, the desirability cache it reads, the
   profile format it adds and the switch between weighting modes. The panel and the messages are covered by the
   manual scenarios in docs/campo_spec_en.md. */
category('MPO: correlation-aware weights', () => {
  test('Four properties, one independent and three nearly identical (paper)', async () => {
    const w = weightsOf(corr(4, [0, 0, 0, 0.98, 0.98, 0.98]), 4);
    expectFloat(w[0], 0.497, 1e-3, 'A');
    for (let i = 1; i < 4; i++)
      expectFloat(w[i], 0.168, 1e-3, `B..D ${i}`);
  });

  test('Three properties, two negatively correlated (paper)', async () => {
    const w = weightsOf(corr(3, [0, 0, -0.35]), 3);
    expectFloat(w[0], 0.245, 1e-3, 'A');
    expectFloat(w[1], 0.377, 1e-3, 'B');
    expectFloat(w[2], 0.377, 1e-3, 'C');
  });

  test('Independent properties get equal weights', async () => {
    const d = 5;
    const w = weightsOf(corr(d, new Array(d * (d - 1) / 2).fill(0)), d);
    for (let i = 0; i < d; i++)
      expectFloat(w[i], 1 / d, 1e-9);
  });

  test('Negative raw weights are clipped and the rest renormalised', async () => {
    // C is almost A + B: its raw weight goes negative
    const m = 3;
    const inv = robustInverse(corr(m, [0.2, 0.7, 0.8]), m)!;
    const {raw, base, clipped} = campoWeights(inv.Rinv, m, Int32Array.from([0, 1, 2]), m);
    expect(clipped, true, 'clipped');
    expect(raw[2] < 0, true, 'C raw is negative');
    expect(base[2], 0, 'C is dropped');
    expectFloat(base[0] + base[1] + base[2], 1, 1e-12, 'sum');
  });

  test('pearsonMatrix matches the Correlation plot on complete rows', async () => {
    const n = 500;
    const cols = correlatedColumns(n, 4, 7);
    const R = pearsonMatrix(cols, n);
    const df = DG.DataFrame.fromColumns(cols.map((c, k) => DG.Column.fromFloat32Array(`c${k}`, c)));
    const plot = DG.Viewer.fromType(DG.VIEWER.CORR_PLOT, df) as DG.CorrelationPlot;
    for (let p = 0; p < 4; p++) {
      expect(R[p * 4 + p], 1, 'diagonal');
      for (let q = 0; q < 4; q++) {
        expect(R[p * 4 + q], R[q * 4 + p], 'symmetric');
        if (p !== q)
          expectFloat(R[p * 4 + q], plot.getCorrelation(df.columns.byIndex(p), df.columns.byIndex(q)), 1e-6);
      }
    }
  });

  test('Cholesky rejects a matrix that is not positive definite; robustInverse regularises', async () => {
    expect(cholesky(corr(2, [1.5]), 2), null, 'r > 1');
    const inv = robustInverse(corr(3, [1, 0.3, 0.3]), 3);
    expect(inv != null, true, 'two identical columns');
    expect(inv!.lambda > 0, true, 'regularised');
  });

  test('finalWeights: the pin example of the spec', async () => {
    const base = Float64Array.from([0.258, 0.038, 0, 0.023, 0, 0.103, 0.363, 0, 0.215]);
    const pinned = new Float64Array(9).fill(NaN);
    pinned[2] = 0.08;
    const {final, error} = finalWeights(base, new Float64Array(9).fill(1), pinned);
    expect(error, null);
    expectFloat(final![2], 0.08, 1e-12, 'HBA pinned');
    expectFloat(final![0], 0.258 * 0.92 / 1.0, 1e-3, 'TPSA');
    expectFloat(final![6], 0.334, 1e-3, 'cLogD');
    expectFloat(final!.reduce((a, b) => a + b, 0), 1, 1e-9, 'sum');
  });

  test('finalWeights: priorities, pins summing to 1 and to more than 1', async () => {
    const base = Float64Array.from([0.5, 0.25, 0.25]);
    const none = new Float64Array(3).fill(NaN);
    const doubled = finalWeights(base, Float64Array.from([2, 1, 1]), none).final!;
    expectFloat(doubled[0], 2 / 3, 1e-12, 'priority 2');
    expectFloat(doubled[1], 1 / 6, 1e-12);

    const full = finalWeights(base, Float64Array.from([1, 1, 1]), Float64Array.from([0.6, 0.4, NaN])).final!;
    expect(full[2], 0, 'pins take all the weight');

    const over = finalWeights(base, Float64Array.from([1, 1, 1]), Float64Array.from([0.7, 0.4, NaN]));
    expect(over.final, null);
    expect(over.error, 'Pinned weights sum to 1.100 (> 1)');

    expect(finalWeights(base, new Float64Array(3), none).error, 'All weighted priorities are 0');
    expect(finalWeights(Float64Array.from([1, 0, 0]), Float64Array.from([1, 1, 1]), Float64Array.from([0.5, NaN, NaN]))
      .error, 'Nothing to scale: all unpinned weights are 0');
  });

  test('Constant, unmapped and incomplete columns are handled without NaN', async () => {
    const n = 50;
    const [a, b] = correlatedColumns(n, 2, 3);
    const constant = new Float32Array(n).fill(1);
    const prep = prepareDesirabilities([a, null, constant, b], n);
    expect(Array.from(prep.active).join(','), '0,3', 'unmapped and constant columns are left out');
    expect(prep.rows, n);

    const holes = Float32Array.from(a);
    holes[0] = holes[5] = NaN;
    expect(prepareDesirabilities([holes, b], n).rows, n - 2, 'listwise deletion');

    const empty = prepareDesirabilities([Float32Array.from(a, () => NaN), b], n);
    expect(empty.rows, 0, 'no complete rows');
    expect(empty.active.length, 0, 'no active columns');

    const inv = robustInverse(pearsonMatrix(prep.cols, prep.rows), 2)!;
    const w = campoWeights(inv.Rinv, 2, prep.active, 4);
    expect(w.base.some((x) => Number.isNaN(x)), false, 'no NaN in base weights');
    expect(Number.isNaN(w.raw[1]) && Number.isNaN(w.raw[2]), true, 'raw is NaN outside R');
  });

  test('Average score with caMPO weights is the weighted sum on complete rows', async () => {
    const n = 200;
    const cols = correlatedColumns(n, 3, 11);
    const line = [[0, 0], [1, 1]];
    const R = pearsonMatrix(cols, n);
    const inv = robustInverse(R, 3)!;
    const {base} = campoWeights(inv.Rinv, 3, Int32Array.from([0, 1, 2]), 3);
    const dgCols = cols.map((c, k) => {
      const col = DG.Column.fromFloat32Array(`p${k}`, c);
      col.setTag('desirabilityTemplate', JSON.stringify(numerical(line, {weight: base[k]})));
      return col;
    });
    const score = mpo(DG.DataFrame.fromColumns(dgCols), dgCols, 'campo', 'Average').scoreColumn;
    for (const i of [0, 17, 199]) {
      const expected = cols.reduce((s, c, k) => s + base[k] * c[i], 0);
      expectFloat(score.get(i), expected, 1e-5, `row ${i}`);
    }
  });

  test('Stage A on 100 000 rows and 20 properties', async () => {
    const n = 100000;
    const d = 20;
    const cols = correlatedColumns(n, d, 5);
    const df = DG.DataFrame.fromColumns(cols.map((c, k) => DG.Column.fromFloat32Array(`p${k}`, c)));
    const tagged = applyDesirabilityTags(df, Object.fromEntries(cols.map((_, k) =>
      [`p${k}`, numerical([[0, 0], [1, 1]])])), true);
    const maps = new MpoCalculator().desirabilities(df, tagged);
    const prep = prepareDesirabilities(maps.map((m) => m.D), n);
    const inv = robustInverse(pearsonMatrix(prep.cols, prep.rows), prep.active.length);
    const {base} = campoWeights(inv!.Rinv, prep.active.length, prep.active, d);
    expectFloat(base.reduce((a, b) => a + b, 0), 1, 1e-9, 'weights sum to 1');
  }, {benchmark: true});
});

category('MPO: correlation-aware plumbing', () => {
  const profileOf = (weighting?: DesirabilityProfile['weighting']): DesirabilityProfile => ({
    type: DESIRABILITY_PROFILE_TYPE, name: 'campo', description: '', weighting,
    properties: {
      a: numerical([[0, 0], [1, 1]], {weight: 0.4, priority: 2}),
      b: numerical([[0, 1], [1, 0]], {weight: 0.6, pinnedWeight: 0.6}),
    },
  });
  const rowOf = (p: DesirabilityProfile, id: string) => ({...toRow(p), id, version: 1}) as ProfileRow;

  test('MpoCalculator.desirabilities reads the cache compute() fills', async () => {
    const values = [0.1, DG.FLOAT_NULL, 0.9, 0.5];
    const col = DG.Column.fromFloat32Array('x', Float32Array.from(values));
    const df = DG.DataFrame.fromColumns([col]);
    const calc = new MpoCalculator();
    const tag = (mv: NumericalDesirability['missingValues']) =>
      col.setTag('desirabilityTemplate', JSON.stringify(numerical([[0, 0], [1, 1]], {missingValues: mv})));

    tag({strategy: 'skip'});
    calc.compute(df, [col], 'calc', 'Average', true);
    const first = calc.desirabilities(df, [col])[0];
    expect(calc.desirabilities(df, [col])[0].D === first.D, true, 'same cached array');
    expect(Number.isNaN(first.D[1]), true, 'skip → NaN');

    tag({strategy: 'exclude'});
    expect(Number.isNaN(calc.desirabilities(df, [col])[0].D[1]), true, 'exclude → NaN');
    tag({strategy: 'default', score: 0.3});
    expectFloat(calc.desirabilities(df, [col])[0].D[1], 0.3, 1e-6, 'default → its score');
  });

  test('stageAKey changes with a curve or a new column, not with weights, priorities or pins', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromFloat32Array('a', Float32Array.from([0, 1]))]);
    const profile = profileOf('correlation');
    const key = stageAKey(df, profile, {});
    profile.properties.a.weight = 0.9;
    profile.properties.a.priority = 3;
    profile.properties.b.pinnedWeight = 0.1;
    expect(stageAKey(df, profile, {}), key, 'weights, priorities and pins');
    (profile.properties.a as NumericalDesirability).line = [[0, 0], [0.5, 1]];
    const curveKey = stageAKey(df, profile, {});
    expect(curveKey !== key, true, 'curve');
    df.columns.add(DG.Column.fromFloat32Array('b', Float32Array.from([1, 0])));
    expect(stageAKey(df, profile, {}) !== curveKey, true, 'column added');
    expect(desirabilityKey(profile.properties.a).includes('priority'), false);
  });

  test('Profile row round trip keeps weighting, priorities, pins and weights', async () => {
    const back = fromRow(rowOf(profileOf('correlation'), 'x'));
    expect(back.weighting, 'correlation');
    expect(back.properties.a.priority, 2);
    expect(back.properties.b.pinnedWeight, 0.6);
    expect(back.properties.b.weight, 0.6);

    expect(toRow(profileOf('manual')).weighting ?? null, null, 'manual is not written');
    expect(fromRow(rowOf(profileOf(), 'y')).weighting == null, true, 'absent = manual');
  });

  test('Renaming a property keeps its priority and pin', async () => {
    const editor = new MpoProfileEditor(undefined, true);
    editor.weightingControl = true;
    const profile = profileOf('correlation');
    editor.setProfile(profile);
    expect(editor.renameProperty('b', 'b2'), true);
    expect(profile.properties.b2.pinnedWeight, 0.6);
    expect(profile.properties.b2.priority == null, true);
    expect(editor.renameProperty('a', 'a2'), true);
    expect(profile.properties.a2.priority, 2);
  });

  test('setWeights fires a change only when a weight changes at display precision', async () => {
    const editor = new MpoProfileEditor(undefined, true);
    editor.weightingControl = true;
    editor.setProfile(profileOf('correlation'));
    let changes = 0;
    const sub = editor.onChanged.subscribe(() => changes++);
    try {
      expect(editor.setWeights({a: 0.4001, b: 0.6}, null), false, 'same at 0.000');
      expect(changes, 0);
      expect(editor.setWeights({a: 0.45, b: 0.55}, null), true, 'changed');
      expect(changes, 1);
    } finally {
      sub.unsubscribe();
    }
  });

  test('Hosts without the Weighting control show computed weights read-only', async () => {
    const editor = new MpoProfileEditor();
    editor.setProfile(profileOf('correlation'));
    const inputs = Array.from(editor.root.querySelectorAll('.statistics-mpo-weight-input input')) as HTMLInputElement[];
    expect(inputs.length > 0, true);
    expect(inputs.every((i) => i.disabled), true, 'disabled, the slider included');

    editor.setProfile(profileOf());
    const manual = Array.from(editor.root.querySelectorAll('.statistics-mpo-weight-input input')) as HTMLInputElement[];
    expect(manual.every((i) => !i.disabled), true, 'manual stays editable');
  });

  test('Switching the weighting moves the weights between modes', async () => {
    const profile: DesirabilityProfile = {
      type: DESIRABILITY_PROFILE_TYPE, name: 'campo-switch', description: '',
      properties: {
        a: numerical([[0, 0], [1, 1]], {weight: 0.4, weightColumn: 'wa'}),
        b: numerical([[0, 1], [1, 0]], {weight: 1}),
      },
    };
    const view = new MpoProfileCreateView(profile, false);
    try {
      const {a, b} = view.profile.properties;
      const input = view.editor.weightingInput;
      input.value = MPO_WEIGHTING_CAPTIONS.correlation;
      expect(view.editor.weighting, 'correlation', 'the input switches the mode');
      expect(a.priority, 0.4, 'a manual weight becomes the priority');
      expect(b.priority == null, true, 'weight 1 is no priority');
      expect(a.weightColumn == null, true, 'no weight column outside Manual');

      view.setWeighting('equal');
      expect(input.value, MPO_WEIGHTING_CAPTIONS.equal, 'the input follows the mode');
      expect(a.weight, 1);
      expect(b.weight, 1);

      view.setWeighting('manual');
      expect(a.weight, 0.4, 'manual weights come back');
      expect(a.weightColumn, 'wa');
      expect(b.weight, 1);
    } finally {
      view.activeView.close();
    }
  });
});
