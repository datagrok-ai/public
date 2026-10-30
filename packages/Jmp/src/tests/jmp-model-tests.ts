import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';
import {JslNode, findAll, parseJsl} from '../formats/jmp/jsl';
import {readJmp} from '../formats/jmp/jmp-parser';
import {modelingType} from '../formats/jmp/jmp-types';
import {ModelColumn, UnfittableModel, modelColumn, modelDesign, modelResponses} from '../formats/jmp/jmp-model';

type Call = Extract<JslNode, {type: 'call'}>;

const COLUMNS: ModelColumn[] = [
  {name: 'y', values: Float64Array.of(1, 2, 3, 4, NaN, 6)},
  {name: 'x', values: Float64Array.of(1, 2, 3, 4, 5, 9)},
  {name: 'z', values: Float64Array.of(0, 10, 0, 10, 0, 10), coding: [0, 10]},
  {name: 'g', levels: ['a', 'b', 'c'], codes: Int32Array.of(0, 1, 2, 0, 1, 2)},
  {name: 'site', levels: ['A', 'B'], codes: Int32Array.of(0, 0, 0, 1, 0, 0)},
  {name: 'v', values: Float64Array.of(0, 0, 1, 0, 0, 0)},
];

function design(text: string) {
  const call = parseJsl(text) as Call;
  const column = (n: string) => COLUMNS.find((c) => c.name === n) ?? (() => {
    throw new Error(`no ${n}`);
  })();
  return modelDesign(call, modelResponses(call)[0], 6, column);
}

function refused(text: string): string {
  try {
    design(text);
  } catch (e) {
    if (e instanceof UnfittableModel)
      return e.message;
    throw e;
  }
  return '';
}

category('JMP: models', () => {
  test('a design as JMP builds it: effect coding, coded and centered factors', async () => {
    const d = design('Fit Model(Y(:y), Effects(:x, :g, :z, :x * :z, :x * :x & RS))');
    expectArray(d.terms, ['Intercept', 'x', 'g[a]', 'g[b]', 'z', 'x*z', 'x*x']);
    const row = (r: number) => Array.from(d.x.subarray(r * 7, r * 7 + 7));
    // the mean of x over the rows fit (all but row 4, whose y is missing) is 19/5
    expectArray(row(0), [1, 1, 1, 0, -1, (1 - 3.8) * -1, (1 - 3.8) * (1 - 3.8)]);
    expectArray(row(2).slice(2, 4), [-1, -1]);
    expectArray(Array.from(d.fit), [1, 1, 1, 1, 0, 1]);
    const raw = design('Fit Model(Y(:y), Effects(:x, :x * :x), Center Polynomials(0), No Intercept(1))');
    expectArray(raw.terms, ['x', 'x*x']);
    expectArray(Array.from(raw.x.subarray(2, 4)), [2, 4]);
  });

  test('the rows fit: Where(), the training rows of a Validation(), a transformed response', async () => {
    expectArray(Array.from(design('Fit Model(Y(:y), Effects(:x), Where(:site == "A"))').fit), [1, 1, 1, 0, 0, 1]);
    expectArray(Array.from(design('Fit Model(Y(:y), Effects(:x), Validation(:v))').fit), [1, 1, 0, 1, 0, 1]);
    expect(design('Fit Model(Y(Log(:y)), Effects(:x))').y[2], Math.log(3));
  });

  test('models plain least squares can\'t fit are refused, with a reason', async () => {
    expect(refused('Fit Model(Y(:y), Effects(:x), Random Effects(:g))'),
      'Random effects need a mixed model, which isn\'t converted yet');
    expect(refused('Fit Model(Y(:y), Effects(:x, :site[:g] & Random))'),
      'Random effects need a mixed model, which isn\'t converted yet');
    expect(refused('Fit Model(Y(:y), Effects(:x, :site[:g]))'), 'Nested effects aren\'t converted yet');
    expect(refused('Fit Model(Y(:g), Effects(:x))'),
      'The response "g" is categorical, which least squares doesn\'t fit');
    expect(refused('Fit Model(Y(:y), Effects(:x), Freq(:v))'), 'The Freq role isn\'t applied in a refit');
  });

  test('JMP 18: a designed experiment\'s factors, coded as its Coding property says', async () => {
    const t = await readJmp(await _package.files.readAsBytes('public/mdpi-bioengineering-11-01089-wp1.jmp'));
    const pp1 = t.columns.find((c) => c.info.name === 'PP 1')!;
    expect(modelingType(pp1), 'continuous');
    const column = modelColumn(pp1);
    expectArray('coding' in column ? column.coding! : [], [-1, 0.38]);
    const fit = findAll(t.scripts.find((s) => s.name === 'Fit OLS')!.body!, 'Fit Model')[0];
    const byName = (n: string) => modelColumn(t.columns.find((c) => c.info.name === n)!);
    const d = modelDesign(fit, modelResponses(fit)[0], t.info.nrows, byName);
    const at = d.terms.indexOf('PP 1');
    const raw = 'values' in column ? column.values[0] : NaN;
    // JMP's saved formula writes the factor as (:PP 1 - (-0.31)) / 0.69
    expect(Math.abs(d.x[at] - (raw + 0.31) / 0.69) < 1e-12, true);
  });
}, {node: true});
